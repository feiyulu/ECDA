module eakf_oda_mod

  ! this module is produced by snz for parallelism of ensemble-based filtering
  ! data assimilation algorithm, starting from August 9, 2002. A linear
  ! regression eakf is used as a primary algorithm to design the parallelism.
  ! Communications are broadcasting the information at the observation location.

  ! FMS shared modules
  use fms_mod, only : open_namelist_file, check_nml_error, close_file
  use fms_mod, only : stdout, error_mesg, FATAL, WARNING
  use mpp_mod, only : mpp_sync_self, pe=>mpp_pe, npes=>mpp_npes
  use mpp_mod, only : mpp_root_pe
  use time_manager_mod, only : time_type, get_time
  use constants_mod, only : DEG_TO_RAD, RADIUS
  use mpp_domains_mod, only : mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain
  use mpp_domains_mod, only : domain2d, mpp_update_domains
  use gsw_mod_toolbox, only : gsw_pt_from_t

  ! ODA Modules
  use ocean_da_types_mod, only : ocean_profile_type, TEMP_ID, SALT_ID, missing_value
  use ocean_da_types_mod, only : ODA_PFL, ODA_XBT, ODA_MRB, ODA_OISST
  use ocean_da_types_mod, only : ocean_control_struct, grid_type
  use kdtree, only : kd_root, kd_search_radius, kd_init

  ! EAKF modules
  use assim_tools_mod, only : assim_tools_init, obs_increment, update_from_inc
  use loc_and_dist_mod, only : loc_type
  use cov_cutoff_mod, only : comp_cov_factor

  implicit none

  logical, save :: first_run_call = .true.
  logical, save :: module_initialized = .false.

  integer, save :: prof_count=0, outlier_count=0
  real, allocatable, dimension(:,:) :: ens
  real, allocatable, dimension(:) :: ens_mean, ens_inc
  real, allocatable, dimension(:) :: enso_temp, inc_temp
  real, allocatable, dimension(:) :: enso_salt, inc_salt
  real, allocatable, dimension(:) :: enso_theta
  integer, allocatable, dimension(:) :: lon1d, lat1d
  integer, allocatable, dimension(:) :: kd_ind
  integer, allocatable, dimension(:) :: dist_seq
  real, allocatable, dimension(:) :: kd_dist
  real, allocatable, dimension(:) :: dist_sorted

  public ensemble_filter

contains

  subroutine ensemble_filter(Prior, Posterior, Profiles, kdroot, Domain, oda_grid)
    type(ocean_control_struct), pointer, intent(in) :: Prior
    type(ocean_control_struct), pointer, intent(inout) :: Posterior
    type(ocean_profile_type), pointer, intent(in) :: Profiles
    type(kd_root), pointer, intent(inout) :: kdroot
    type(domain2d), pointer :: Domain
    type(grid_type), pointer, intent(in) :: oda_grid

    !---- namelist with default values
    real :: e_flder_oer = 2000.0
    real :: depth_eq = 200.0
    real :: lat_eq = 20.0
    logical :: debug_eakf = .false.
    logical :: outlier_qc = .true.
    logical :: get_obs_forecast = .true.
    logical :: get_obs_analysis = .false.
    real :: sst_ice_limit = -2.0
    real :: obs_ice_limit = -2.0
    real :: temp_limit = 5.0
    real :: salt_limit = 2.0

    namelist /eakf_nml/ e_flder_oer, depth_eq, lat_eq, debug_eakf, outlier_qc, &
            get_obs_forecast, get_obs_analysis, sst_ice_limit, obs_ice_limit, &
            temp_limit, salt_limit

    !--- module name and version number ----
    !character(len=*), parameter :: module_name = 'eakf'
    !character(len=*), parameter :: vers_num = 'x00.0'

    !=======================================================================
    integer :: assim_var = 1 ! 1 for self update, 2 for cross-variable update
    integer :: flag_hyb = 1 ! 0 for hyb scheme, 1 for eakf
    integer :: ni, nj, nk, ndim
    integer :: stdout_unit
    integer :: ens_size
    type(ocean_profile_type), pointer :: Prof

    !=======================================================================
    real, dimension(:), allocatable :: glon1d, glat1d
    integer :: kd_num

    !=======================================================================
    integer :: id_eakf_total
    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: isg, ieg, jsg, jeg
    integer :: halox, haloy

    real :: cov_factor, cov_factor_v, cov_factor_t, cov_factor_h

    logical :: assim_flag, interp_flag
    integer :: model_basin
    type(loc_type) :: model_loc, obs_loc

    integer :: ii_ens, jj_ens, kk_ens
    integer :: ind_temp_h, ind_temp_l, ind_salt_h, ind_salt_l
    integer :: i, j, k, k0, kk, num, blk, i_idx, salt_offset
    integer :: t_tau, s_tau
    integer :: kk0, kk1, kk2
    integer :: model_size
    integer :: unit, ierr, io, j_ens, i_h
    integer :: diff_days, diff_seconds, window_days, window_seconds
    real :: diff_hours, diff_k, window_hours

    !---------------------------------------------------------------------------
    real :: forecast_t, forecast_s, analysis_t, analysis_s
    real :: obs_value, obs_sigma, obs_var
    real :: dist, dist0
    real :: v2_h, v2_l
    real :: depth_bot, depth_kk
    real :: obs_dist = 5.0
    real :: outlier_limit = 5.0

    !---------------------------------------------------------------------------
    character(len=40) :: file_name
    character(len=40) :: diag_file

    stdout_unit = stdout()
    !---------------------------------------------------------------------------
    
    if((.not. associated(Prior)) .or. (.not. associated(Posterior)) ) &
            call error_mesg('eakf_oda_mod::ensemble_filter', 'Prior or Posterior not initialized', FATAL)
    if( Prior%ensemble_size .ne. Posterior%ensemble_size ) &
            call error_mesg('eakf_oda_mod::ensemble_filter', 'Ensemble size mismatch in Prior and Posterior', FATAL)

    ens_size = Prior%ensemble_size
    ndim = size(shape(Prior%T))
    nk = size(Prior%T, ndim-1)

    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain, isd, ied, jsd, jed)
    call mpp_get_global_domain(Domain, isg, ieg, jsg, jeg)
    isc = isc - isd + 1 ; iec = iec - isd + 1 ; jsc = jsc - jsd + 1 ; jec = jec - jsd + 1
    ied = ied - isd + 1 ; jed = jed - jsd + 1 ; isd = 1 ; jsd = 1
    halox = isc - isd
    haloy = jsc - jsd
    ni = ieg-isg+1; nj = jeg-jsg+1
    blk = (jed-jsd+1)*(ied-isd+1)
    salt_offset = nk*blk

    ! Read namelist for run time control
    unit = open_namelist_file()
    read(unit, nml = eakf_nml, iostat=io)
    call close_file(unit)

    if ( check_nml_error(io, 'eakf_nml') < 0 ) then
       if ( pe() == mpp_root_pe() ) then
          call error_mesg('eakf_mod::ensemble_filter', &
                  'EAKF_NML not found in input.nml.  Using defaults.', WARNING)
       end if
    end if

    model_size = blk * nk * 2
    if ( .not. module_initialized ) then
       call eakf_oda_init(model_size, ens_size, blk)
    
       allocate( glon1d(blk), glat1d(blk) )
       do i = isd, ied; do j = jsd, jed
         lon1d((j-jsd)*(ied-isd+1)+i-isd+1) = i
         lat1d((j-jsd)*(ied-isd+1)+i-isd+1) = j
         glon1d((j-jsd)*(ied-isd+1)+i-isd+1) = oda_grid%x(i,j)
         glat1d((j-jsd)*(ied-isd+1)+i-isd+1) = oda_grid%y(i,j)
       enddo; enddo

       if ( .not. associated(kdroot) ) then
         allocate(kdroot)
         call kd_init(kdroot, glon1d, glat1d)
       endif

       deallocate( glon1d, glat1d )
    end if

    ! Initialize assim tools module
    call assim_tools_init()

    ! print namelist
    if ( pe() == mpp_root_pe() .and. first_run_call ) then
       write (stdout_unit, *) 'model size is ', model_size, 'ensemble size is ', ens_size
       write (stdout_unit, *) 'e_flder_oer is', e_flder_oer
       write (stdout_unit, *) 'outlier qc switch', outlier_qc
       write (stdout_unit, *) 'outlier qc limit for T/S', temp_limit, salt_limit
       write (stdout_unit, *) 'Eq. lat and depth limit', lat_eq, depth_eq
    end if

    if ( debug_eakf ) then
       write (UNIT=stdout_unit, FMT='("PE ",I5,": finished eakf initialization")') pe()
    end if

    ! Form the ensemble state vector: ens(:, :)
    call ens_ics(Prior%T, Prior%S, isd, ied, jsd, jed, halox, haloy, nk, ens)

    if ( debug_eakf ) then
       write (UNIT=stdout_unit, FMT='("PE ",I5,": finished ens_ics")') pe()
    end if

    ! ###########################################################
    ! The assimilation main part starts here

    ! Compute the ensemble mean of the initial ensemble before assimilation
    ens_mean = sum(ens, dim=2) / ens_size
    call mpp_sync_self()

    !===== Eakf assim start =====================================

    call mpp_sync_self()

    if ( debug_eakf ) then
       write (UNIT=stdout_unit, FMT='("PE ",I5,": start profile loop")') pe()
    end if

    if ( get_obs_forecast ) then
       Prof => Profiles
       do while (associated(Prof))
          k0 = Prof%levels
          if(.not. associated(Prof%forecast)) allocate(Prof%forecast(k0))
          do kk = 1, k0
             Prof%forecast(kk) = missing_value
             depth_kk = Prof%depth(kk)
             interp_flag = .true.
             v2_h = 0.0
             v2_l = 0.0
             if ( Prof%variable == TEMP_ID ) then
                do i_idx=1, 4
                   ind_temp_h = Prof%obs_def(kk)%state_var_index(i_idx)
                   ind_temp_l = Prof%obs_def(kk)%state_var_index(i_idx+4)
                   v2_h = v2_h + ens_mean(ind_temp_h)*Prof%obs_def(kk)%coef(i_idx)
                   v2_l = v2_l + ens_mean(ind_temp_l)*Prof%obs_def(kk)%coef(i_idx)
                   if(ens_mean(ind_temp_h).eq.0.0 .or. ens_mean(ind_temp_l).eq.0.0) then
                      interp_flag = .false.
                   end if
                end do
                forecast_t = v2_h*Prof%obs_def(kk)%coef(5) + v2_l*Prof%obs_def(kk)%coef(6)
             end if
             v2_h = 0.0
             v2_l = 0.0
             if ( Prof%variable == SALT_ID .or. Prof%variable == TEMP_ID ) then
                do i_idx=1, 4
                   ind_salt_h = Prof%obs_def(kk)%state_var_index(i_idx)+salt_offset
                   ind_salt_l = Prof%obs_def(kk)%state_var_index(i_idx+4)+salt_offset
                   v2_h = v2_h + ens_mean(ind_salt_h)*Prof%obs_def(kk)%coef(i_idx)
                   v2_l = v2_l + ens_mean(ind_salt_l)*Prof%obs_def(kk)%coef(i_idx)
                   if(ens_mean(ind_salt_h).eq.0.0 .or. ens_mean(ind_salt_l).eq.0.0) then
                      interp_flag = .false.
                   end if
                end do
                forecast_s = v2_h*Prof%obs_def(kk)%coef(5) + v2_l*Prof%obs_def(kk)%coef(6)
             end if
             if(interp_flag) then
                if ( Prof%variable == TEMP_ID ) then
                   Prof%forecast(kk) = gsw_pt_from_t(forecast_s,forecast_t,0.0,depth_kk)
                   if (Prof%inst_type == ODA_OISST) then
                      if ( Prof%forecast(kk) < sst_ice_limit) Prof%flag(kk) = .false.
                      if ( Prof%data(kk) < obs_ice_limit) Prof%flag(kk) = .false.
                   endif
                elseif ( Prof%variable == SALT_ID ) then
                   Prof%forecast(kk) = forecast_s
                end if
             else
                Prof%flag(kk) = .false.
             endif
          end do
          Prof => Prof%cnext
       end do
    end if

    Prof => Profiles
    doloop_9: do while (associated(Prof)) ! (9)
       k0 = Prof%levels
       depth_bot = Prof%depth(k0)
       if ( Prof%variable == TEMP_ID ) outlier_limit=temp_limit
       if ( Prof%variable == SALT_ID ) outlier_limit=salt_limit

       call get_time(Prof%tdiff, diff_seconds, diff_days)
       diff_hours = diff_seconds/3600 + diff_days * 24
       call get_time(Prof%time_window, window_seconds, window_days)
       window_hours = window_seconds/3600 + window_days * 24
       obs_loc%lon = Prof%lon
       obs_loc%lat = Prof%lat

       if ( abs(obs_loc%lat) < 60.0 ) then
          dist0 = Prof%loc_dist*cos(obs_loc%lat*DEG_TO_RAD)
       else
          dist0 = Prof%loc_dist*cos(60.0*DEG_TO_RAD)
       end if
       call kd_search_radius(kdroot, obs_loc%lon, obs_loc%lat, &
               2*dist0, kd_ind, kd_dist, kd_num, .true.)

       dist_sorted(1:kd_num) = kd_dist(1:kd_num)
       dist_seq = 0
       call hpsort_eps_epw (kd_num, dist_sorted(1:kd_num), dist_seq(1:kd_num))

       doloop_4: do kk = 1, k0 ! (4)
          if ( Prof%flag(kk) ) then ! add each level flag
             depth_kk = Prof%depth(kk)
             obs_sigma = Prof%obs_error*exp(-Prof%depth(kk)/e_flder_oer)
             obs_value = Prof%data(kk)
             obs_var = obs_sigma * obs_sigma
             
             kk0 = FLOOR(Prof%k_index(kk))
             kk1 = kk0 - 2 * Prof%impact_levels +1
             kk2 = kk0 + 2 * Prof%impact_levels
             if(Prof%inst_type .eq. ODA_PFL .and. kk.eq.k0 .AND. depth_bot.gt.1500.0) kk2 = nk
             if(kk1 < 1) kk1 = 1
             if(kk2 > nk) kk2 = nk

             doloop_8: do k=1, kd_num ! (8)
                ii_ens = lon1d(kd_ind(dist_seq(k)))
                jj_ens = lat1d(kd_ind(dist_seq(k)))
                i_h = (jj_ens-jsd)*(ied-isd+1)+ii_ens-isd+1

                model_loc%lon = oda_grid%x(ii_ens, jj_ens)
                model_loc%lat = oda_grid%y(ii_ens, jj_ens)

                assim_flag = .false.
          
                model_basin = oda_grid%basin_mask(ii_ens, jj_ens)
                if(Prof%basin_mask .eq. 1) then ! Southern ocean
                   if(model_basin == 1 .or. model_basin == 2 .or. &
                           model_basin == 3 .or. model_basin == 5) assim_flag = .true.
                elseif(Prof%basin_mask .eq. 2) then ! Atlantic
                   if(model_basin == 2 .or. model_basin == 1 .or. model_basin == 4) assim_flag = .true.
                elseif(Prof%basin_mask .eq. 3) then ! Pacific
                   if(model_basin == 3 .or. model_basin == 1) assim_flag = .true.
                elseif(Prof%basin_mask .eq. 4) then ! Arctic
                   if(model_basin == 4 .or. model_basin == 2) assim_flag = .true.
                elseif(Prof%basin_mask .eq. 5) then ! Indian
                   if(model_basin == 5 .or. model_basin == 1) assim_flag = .true.
                elseif(Prof%basin_mask .eq. 6) then ! Mediterranean
                   if(model_basin == 6) assim_flag = .true.
                end if

                ifblock_6: if ( assim_flag ) then ! (6)
                   dist = kd_dist(dist_seq(k))
                   cov_factor_h = comp_cov_factor(dist, dist0)*cos((model_loc%lat-obs_loc%lat)*DEG_TO_RAD)
                   cov_factor_t = comp_cov_factor(diff_hours, window_hours)

                   !if (abs(obs_loc%lat)<lat_eq .and. depth_kk<depth_eq) then
                     !if ((Prof%inst_type.eq.ODA_PFL .or. Prof%inst_type.eq.ODA_XBT) .and. window_hours>24.0) then
                       !cov_factor_t = comp_cov_factor(diff_hours, 24.0)
                     !endif
                     !if (dist>200.0e3) then
                       !cov_factor_h = comp_cov_factor(200.0e3, dist0)*cos((model_loc%lat-obs_loc%lat)*DEG_TO_RAD)
                     !endif
                   !endif
                   
                   do j_ens=1, ens_size
                      v2_h = 0.0
                      v2_l = 0.0
                      if ( Prof%variable == TEMP_ID ) then
                         do i_idx=1, 4
                            ind_temp_h = Prof%obs_def(kk)%state_var_index(i_idx)
                            ind_temp_l = Prof%obs_def(kk)%state_var_index(i_idx+4)
                            v2_h = v2_h + ens(ind_temp_h, j_ens)*Prof%obs_def(kk)%coef(i_idx)
                            v2_l = v2_l + ens(ind_temp_l, j_ens)*Prof%obs_def(kk)%coef(i_idx)
                         end do
                         enso_theta(j_ens) = v2_h*Prof%obs_def(kk)%coef(5) &
                                 + v2_l*Prof%obs_def(kk)%coef(6)
                      end if
                      v2_h = 0.0
                      v2_l = 0.0
                      if ( Prof%variable == TEMP_ID .or. Prof%variable == SALT_ID ) then
                         do i_idx=1, 4
                            ind_salt_h = Prof%obs_def(kk)%state_var_index(i_idx)+salt_offset
                            ind_salt_l = Prof%obs_def(kk)%state_var_index(i_idx+4)+salt_offset
                            v2_h = v2_h + ens(ind_salt_h, j_ens)*Prof%obs_def(kk)%coef(i_idx)
                            v2_l = v2_l + ens(ind_salt_l, j_ens)*Prof%obs_def(kk)%coef(i_idx)
                         end do
                         enso_salt(j_ens) = v2_h*Prof%obs_def(kk)%coef(5) &
                                 + v2_l*Prof%obs_def(kk)%coef(6)
                      end if
                   end do

                   if ( Prof%variable == TEMP_ID ) then
                      do j_ens=1, ens_size
                         enso_temp(j_ens) = gsw_pt_from_t(enso_salt(j_ens),enso_theta(j_ens),0.0,depth_kk)
                      end do
                      inc_temp = 0.0
                      call obs_increment(enso_temp, ens_size, obs_value, obs_var, inc_temp, obs_dist)
                   elseif ( Prof%variable == SALT_ID ) then
                      inc_salt = 0.0
                      call obs_increment(enso_salt, ens_size, obs_value, obs_var, inc_salt, obs_dist)
                   end if

                   prof_count=prof_count+1
                   if (obs_dist>outlier_limit) outlier_count=outlier_count+1
                   if (outlier_qc .and. obs_dist>outlier_limit) then
                      !print *,"outlier obs, var ratio=",obs_dist
                   else
                      doloop_3: do kk_ens=kk1, kk2 ! (3)
                         t_tau   = (kk_ens-1)*blk + i_h
                         s_tau   = salt_offset + t_tau
                         if ( kk_ens <= kk0 ) then
                           diff_k = Prof%k_index(kk) - REAL(kk_ens)
                         else
                           diff_k = REAL(kk_ens) - Prof%k_index(kk)
                         end if
                         cov_factor_v = comp_cov_factor( diff_k, REAL(Prof%impact_levels) )
                         if(kk.eq.k0 .AND. depth_bot.gt.1500.0 .AND. kk_ens.gt.kk0) cov_factor_v  = 1.0
                         cov_factor=abs(cov_factor_v*cov_factor_t*cov_factor_h)

                         ifblock_1: if ( cov_factor > 0.0 ) then ! (1)
                            ifblock_0_60: if ( Prof%variable == TEMP_ID .and. ens_mean(t_tau) /= 0.0) then ! (0.60)
                               ! using temperature adjusts temperature and salinity
                               ! temperature itself
                               assim_var = 1
                               ens_inc = 0.0
                               call update_from_inc(enso_temp, inc_temp, ens(t_tau, :), ens_size, ens_inc, cov_factor, assim_var)
                               ens(t_tau, :) = ens(t_tau, :) + ens_inc(:)

                               ! limit the unreasonable values if applicable
                               doloop_0: do j_ens=1, ens_size ! (0)
                                  if ( ens(t_tau, j_ens) > 39.0 ) then
                                     write (UNIT=stdout_unit, FMT='("T(",3I5,") = ",F15.8)')&
                                          & ii_ens, jj_ens, kk_ens, ens(t_tau,j_ens)
                                     ens(t_tau, j_ens) = 39.0
                                  end if
                                  if ( ens(t_tau, j_ens) < -4.0 ) then
                                     write (UNIT=stdout_unit, FMT='("T(",3I5,") = ",F25.8)')&
                                          & ii_ens, jj_ens, kk_ens, ens(t_tau,j_ens)
                                     ens(t_tau, j_ens) = -4.0
                                  end if
                               end do doloop_0 ! handle the extremeties (0)

                               ! temperature impact salinity
                               ifblock_0_50: if ( Prof%temp_to_salt .and. ens_mean(s_tau) /= 0.0) then ! (0.5)
                                  assim_var = 2
                                  ens_inc = 0.0
                                  call update_from_inc(enso_temp, inc_temp, ens(s_tau, :), ens_size, ens_inc, cov_factor, assim_var)
                                  ens(s_tau, :) = ens(s_tau, :) + ens_inc(:)

                                  ! limit the unreasonable values if applicable
                                  doloop_00: do j_ens = 1, ens_size ! (0)
                                     if ( ens(s_tau, j_ens) > 44.0 ) then
                                        write (UNIT=stdout_unit, FMT='("S(",3I5,") = ",F15.8)') &
                                             & ii_ens, jj_ens, kk_ens, ens(s_tau,j_ens)
                                        ens(s_tau, j_ens) = 44.0
                                     end if

                                     if ( ens(s_tau, j_ens) < 0.0 ) then
                                        write (UNIT=stdout_unit, FMT='("S(",3I5,") = ",F15.8)')&
                                             & ii_ens, jj_ens, kk_ens, ens(s_tau,j_ens)
                                        ens(s_tau, j_ens) = 0.0
                                     end if
                                  end do doloop_00 ! handle the extremeties (0)
                               end if ifblock_0_50 ! impact salinity or not (0.5)
                            
                            elseif ( Prof%variable == SALT_ID .and. ens_mean(s_tau) /= 0.0) then ! (0.60)
                               assim_var = 1
                               ens_inc = 0.0
                               call update_from_inc(enso_salt, inc_salt, ens(s_tau, :), ens_size, ens_inc, cov_factor, assim_var)
                               ens(s_tau, :) = ens(s_tau, :) + ens_inc(:)

                               ! limit the unreasonable values if applicable
                               do j_ens = 1, ens_size ! (0)
                                  if ( ens(s_tau, j_ens) > 44.0 ) then
                                     write (UNIT=stdout_unit, FMT='("S(",3I5,") = ",F15.8)')&
                                          & ii_ens, jj_ens, kk_ens, ens(s_tau,j_ens)
                                     ens(s_tau, j_ens) = 44.0
                                  end if
                                  if ( ens(s_tau, j_ens) < 0.0 ) then
                                     write (UNIT=stdout_unit, FMT='("S(",3I5,") = ",F15.8)')&
                                          & ii_ens, jj_ens, kk_ens, ens(s_tau,j_ens)
                                     ens(s_tau, j_ens) = 0.0
                                  end if
                               end do ! handle the extremeties (0)

                               ! salinity impact temperature
                               ifblock_0_5: if ( Prof%salt_to_temp .and. ens_mean(t_tau) /= 0.0) then ! (0.5)
                                  assim_var = 2
                                  ens_inc = 0.0
                                  call update_from_inc(enso_salt, inc_salt, ens(t_tau, :), ens_size, ens_inc, cov_factor, assim_var)
                                  ens(t_tau, :) = ens(t_tau, :) + ens_inc(:)

                                  ! limit the unreasonable values if applicable
                                  do j_ens = 1, ens_size ! (0)
                                     if ( ens(t_tau, j_ens) > 39.0 ) then
                                        write (UNIT=stdout_unit, FMT='("T(",3I5,") = ",F15.8)')&
                                             & ii_ens, jj_ens, kk_ens, ens(t_tau,j_ens)
                                        ens(t_tau, j_ens) = 39.0
                                     end if
                                     if(ens(t_tau, j_ens) < -4.0) then
                                        write (UNIT=stdout_unit, FMT='("T(",3I5,") = ",F15.8)')&
                                             & ii_ens, jj_ens, kk_ens, ens(t_tau,j_ens)
                                        ens(t_tau, j_ens) = -4.0
                                     end if
                                  end do ! handle the extremeties (0)
                               end if ifblock_0_5 ! impact temperature or not (0.5)
                            end if ifblock_0_60 ! finish processing both variables
                         end if ifblock_1 ! only use obs which has cov_factor > 0.0 (1)
                      end do doloop_3 ! finish adjustments for related model levels (3)
                   end if ! obs_dist
                end if ifblock_6 ! get rid of land points (6)
             end do doloop_8 ! finish the adjustments for all related model columns (8)
          end if ! add each level flag
       end do doloop_4 ! go through one profile column (4)
       Prof => Prof%cnext
    end do doloop_9 ! finish all profiles (9)

    !===== Eakf assim finish =====================================

    ! Compute ensemble mean of current assimilated state
    !
    ! Redistribute the sub ensemble state vector ens(:, :) back to the model grids
    ! in the local-domain.
    call mpp_sync_self()
    call red_ens(Posterior%T, Posterior%S, isd, ied, jsd, jed, halox, haloy, nk, ens)

    do j_ens=1,ens_size
      call mpp_update_domains(Posterior%T(:,:,:,j_ens), Domain)
      call mpp_update_domains(Posterior%S(:,:,:,j_ens), Domain)
    enddo
   
    ens_mean = sum(ens, dim=2) / ens_size
    if ( get_obs_analysis ) then
       Prof => Profiles
       do while (associated(Prof))
          k0 = Prof%levels
          if(.not. associated(Prof%analysis)) allocate(Prof%analysis(k0))
          do kk = 1, k0
             Prof%analysis(kk) = missing_value
             depth_kk = Prof%depth(kk)
             interp_flag = .true.
             v2_h = 0.0
             v2_l = 0.0
             if ( Prof%variable == TEMP_ID ) then
                do i_idx=1, 4
                   ind_temp_h = Prof%obs_def(kk)%state_var_index(i_idx)
                   ind_temp_l = Prof%obs_def(kk)%state_var_index(i_idx+4)
                   v2_h = v2_h + ens_mean(ind_temp_h)*Prof%obs_def(kk)%coef(i_idx)
                   v2_l = v2_l + ens_mean(ind_temp_l)*Prof%obs_def(kk)%coef(i_idx)
                   if(ens_mean(ind_temp_h).eq.0.0 .or. ens_mean(ind_temp_l).eq.0.0) then
                      interp_flag = .false.
                   end if
                end do
                analysis_t = v2_h*Prof%obs_def(kk)%coef(5) + v2_l*Prof%obs_def(kk)%coef(6)
             end if
             v2_h = 0.0
             v2_l = 0.0
             if ( Prof%variable == SALT_ID .or. Prof%variable == TEMP_ID ) then
                do i_idx=1, 4
                   ind_salt_h = Prof%obs_def(kk)%state_var_index(i_idx)+salt_offset
                   ind_salt_l = Prof%obs_def(kk)%state_var_index(i_idx+4)+salt_offset
                   v2_h = v2_h + ens_mean(ind_salt_h)*Prof%obs_def(kk)%coef(i_idx)
                   v2_l = v2_l + ens_mean(ind_salt_l)*Prof%obs_def(kk)%coef(i_idx)
                   if(ens_mean(ind_salt_h).eq.0.0 .or. ens_mean(ind_salt_l).eq.0.0) then
                      interp_flag = .false.
                   end if
                end do
                analysis_s = v2_h*Prof%obs_def(kk)%coef(5) + v2_l*Prof%obs_def(kk)%coef(6)
             end if
             if(interp_flag) then
                if ( Prof%variable == TEMP_ID ) then
                   Prof%analysis(kk) = gsw_pt_from_t(analysis_s,analysis_t,0.0,depth_kk)
                elseif ( Prof%variable == SALT_ID ) then
                   Prof%analysis(kk) = analysis_s
                end if
             endif
          end do
          Prof => Prof%cnext
       end do
    end if


    if ( debug_eakf ) then
       write (UNIT=stdout_unit, FMT='("PE ",I5,": finished red_ens")') pe()
       if(prof_count>0 .and. real(outlier_count)/real(prof_count)>0.05) then
          print *, 'pe=', pe(), 'outliers=', outlier_count, 'ratio=', real(outlier_count)/real(prof_count)
       endif
    end if

    first_run_call = .false.

  end subroutine ensemble_filter


  ! Initialize arrays and other module variables
  subroutine eakf_oda_init(model_size, ens_size, blk)
    integer, intent(in) :: model_size, ens_size, blk

    integer :: istat, j_ens

    character(len=128) :: emsg_local

    if ( module_initialized ) then
       call error_mesg('eakf_oda_mod::eakf_oda_init', 'Module already initialized.', WARNING)
    else
       allocate(ens(model_size, ens_size), STAT=istat);         ens = 0.0
       allocate(ens_mean(model_size), STAT=istat);              ens_mean = 0.0
       allocate(enso_theta(ens_size), STAT=istat);              enso_theta = 0.0
       allocate(enso_temp(ens_size), STAT=istat);               enso_temp = 0.0
       allocate(inc_temp(ens_size), STAT=istat);                inc_temp= 0.0
       allocate(ens_inc(ens_size), STAT=istat);                 ens_inc = 0.0
       allocate(enso_salt(ens_size), STAT=istat);               enso_salt = 0.0
       allocate(inc_salt(ens_size), STAT=istat);                inc_salt = 0.0
       allocate(lon1d(blk), lat1d(blk) )
       allocate(kd_ind(blk) )
       allocate(dist_seq(blk))
       allocate(dist_sorted(blk))
       allocate(kd_dist(blk) )
    end if

    module_initialized = .true.
  end subroutine eakf_oda_init

  subroutine ens_ics(T_ensemble, S_ensemble, isd, ied, jsd, jed, halox, haloy, nk, as)
    real, pointer, intent(in), dimension(:,:,:,:) :: T_ensemble, S_ensemble
    integer, intent(in) :: isd, ied, jsd, jed, halox, haloy, nk
    real, intent(inout), dimension(:,:) :: as

    integer :: m, i, j, k, idx, blk_h

    !  Get initial state for ensemble members
    blk_h = (ied-isd+1)*(jed-jsd+1)
    as(:,:) = 0.0

    do m=1, size(as, 2)
       do k=1, nk
          do j=jsd, jed ! add halo grids 4 each pedomain in y
             do i=isd, ied ! add halo grids 4 each pedomain in x
                idx = (k-1)*blk_h + (j-jsd)*(ied-isd+1) + i-isd+1
                as(idx, m) = T_ensemble(i,j,k,m)
                idx = blk_h*nk + idx
                as(idx, m) = S_ensemble(i,j,k,m)
             end do
          end do
       end do
    end do
  end subroutine ens_ics

  subroutine red_ens(T_ensemble, S_ensemble, isd, ied, jsd, jed, halox, haloy, nk, as)
    real, pointer, intent(inout), dimension(:,:,:,:) :: T_ensemble, S_ensemble
    integer, intent(in) :: isd, ied, jsd, jed, halox, haloy, nk
    real, intent(in), dimension(:,:) :: as

    integer :: m, i, j, k, idx, idx_s, blk_h

    !  Get initial state for ensemble members
    blk_h = (ied-isd+1)*(jed-jsd+1)

    do m=1, size(as, 2)
       do k=1, nk
          do j=jsd+haloy, jed-haloy ! drop halo
             do i=isd+halox, ied-halox ! drop halo
                idx = (k-1)*blk_h + (j-jsd)*(ied-isd+1) + i-isd+1
                T_ensemble(i,j,k,m) = as(idx, m)
                idx_s = blk_h*nk + idx
                S_ensemble(i,j,k,m) = as(idx_s, m)
             end do
          end do
       end do
    end do
  end subroutine red_ens

  subroutine hpsort_eps_epw (n, ra, ind)
    !-input/output variables
    integer, intent(in)   :: n  
    integer, dimension(:) :: ind
    real, dimension(:) :: ra
    !-local variables
    integer :: i, ir, j, l, iind  
    real :: rra  
  !
    ! initialize index array
    IF (ind (1) .eq.0) then  
       DO i = 1, n  
          ind (i) = i  
       ENDDO
    ENDIF
    ! nothing to order
    IF (n.lt.2) return  
    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1  
  
    ir = n  
  
    sorting: do 
    
      ! still in hiring phase
      IF ( l .gt. 1 ) then  
         l    = l - 1  
         rra  = ra (l)  
         iind = ind (l)  
         ! in retirement-promotion phase.
      ELSE  
         ! clear a space at the end of the array
         rra  = ra (ir)  
         !
         iind = ind (ir)  
         ! retire the top of the heap into it
         ra (ir) = ra (1)  
         !
         ind (ir) = ind (1)  
         ! decrease the size of the corporation
         ir = ir - 1  
         ! done with the last promotion
         IF ( ir .eq. 1 ) then  
            ! the least competent worker at all !
            ra (1)  = rra  
            !
            ind (1) = iind  
            exit sorting  
         ENDIF
      ENDIF
      ! wheter in hiring or promotion phase, we
      i = l  
      ! set up to place rra in its proper level
      j = l + l  
      !
      DO while ( j .le. ir )  
         IF ( j .lt. ir ) then  
            ! compare to better underling
            IF ( hslt( ra (j),  ra (j + 1) ) ) then  
               j = j + 1  
            !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
               ! this means ra(j) == ra(j+1) within tolerance
             !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
            ENDIF
         ENDIF
         ! demote rra
         IF ( hslt( rra, ra (j) ) ) then  
            ra (i) = ra (j)  
            ind (i) = ind (j)  
            i = j  
            j = j + j  
         !else if ( .not. hslt ( ra(j) , rra ) ) then
            !this means rra == ra(j) within tolerance
            ! demote rra
           ! if (iind.lt.ind (j) ) then
           !    ra (i) = ra (j)
           !    ind (i) = ind (j)
           !    i = j
           !    j = j + j
           ! else
               ! set j to terminate do-while loop
           !    j = ir + 1
           ! endif
            ! this is the right place for rra
         ELSE
            ! set j to terminate do-while loop
            j = ir + 1  
         ENDIF
      ENDDO
      ra (i) = rra  
      ind (i) = iind  
  
    END DO sorting    
  end subroutine hpsort_eps_epw

  !  internal function 
  !  compare two real number and return the result

  logical function hslt( a, b )
    REAL :: a, b
    IF( abs(a-b) <  10e-4 ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
  end function hslt

end module eakf_oda_mod
