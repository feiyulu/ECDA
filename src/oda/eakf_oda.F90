module eakf_oda_mod

  ! this module is produced by snz for parallelism of ensemble-based filtering
  ! data assimilation algorithm, starting from August 9, 2002. A linear
  ! regression eakf is used as a primary algorithm to design the parallelism.
  ! Communications are broadcasting the information at the observation location.

  ! FMS shared modules
  use fms_mod, only : open_namelist_file, check_nml_error, close_file
  use fms_mod, only : stdout, error_mesg, FATAL, WARNING
  use mpp_mod, only : mpp_sync_self, pe=>mpp_pe, npes=>mpp_npes
  use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_root_pe
  use time_manager_mod, only : time_type, get_time
  use constants_mod, only : DEG_TO_RAD, RADIUS
  use mpp_domains_mod, only : domain2d, mpp_get_data_domain, mpp_get_compute_domain, mpp_get_global_domain

  ! ODA Modules
  use oda_types_mod, only : ocean_profile_type, TEMP_ID, SALT_ID, missing_value
  use oda_types_mod, only : ocean_control_struct, grid_type
  use kdtree, only : kd_root, kd_search_radius, kd_init

  ! EAKF modules
  use assim_tools_mod, only : assim_tools_init, obs_increment_prf_eta_hyb, update_from_obs_inc_prf_hyb
  use loc_and_dist_mod, only : loc_type
  use cov_cutoff_mod, only : comp_cov_factor

  implicit none

  logical, save :: first_run_call = .true.
  logical, save :: module_initialized = .false.

  real, allocatable, dimension(:,:) :: ens
  real, allocatable, dimension(:) :: ens_mean, ens_inc
  real, allocatable, dimension(:) :: enso_temp, obs_inc_eakf_temp, obs_inc_oi_temp
  real, allocatable, dimension(:) :: enso_salt, obs_inc_eakf_salt, obs_inc_oi_salt
  integer, allocatable, dimension(:) :: lon1d, lat1d
  integer, allocatable, dimension(:) :: kd_ind
  real, allocatable, dimension(:) :: kd_dist

  public ensemble_filter

contains

  subroutine ensemble_filter(Prior, Posterior, Profiles, kdroot, Domain, oda_grid)
    type(ocean_control_struct), pointer, intent(in) :: Prior
    type(ocean_control_struct), pointer, intent(inout) :: Posterior
    type(ocean_profile_type), pointer, intent(in) :: Profiles
    type(kd_root), pointer, intent(inout) :: kdroot
    type(domain2d), intent(in) :: Domain
    type(grid_type), pointer, intent(in) :: oda_grid

    !---- namelist with default values
    real :: cov_inflate = 1.0
    real :: mean_inflate = 1.0
    real :: sigma_o_t = 1.0
    real :: sigma_o_s = 1.0
    integer :: impact_levels = 3
    real :: update_window = 240.0
    real :: temp_dist = 1000.0e3
    real :: salt_dist = 500.0e3
    logical :: salt_to_temp = .false.
    logical :: temp_to_salt = .false.
    real :: depth_cut = 2000.0
    real :: e_flder_oer = 2000.0
    logical :: debug_eakf = .false.

    namelist /eakf_nml/sigma_o_t, sigma_o_s, impact_levels, &
            update_window, temp_dist, salt_dist, depth_cut, e_flder_oer, &
            salt_to_temp, temp_to_salt, debug_eakf

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

    real :: cor_oi = 0.0
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
    integer :: unit, ierr, io, j_ens, i_h, kk_bot
    integer :: diff_days, diff_seconds
    real :: diff_hours, diff_k

    !---------------------------------------------------------------------------
    real :: obs_value, obs_sigma_t, obs_sigma_s, obs_var_t, obs_var_s
    real :: obs_var_t_oi = 0.0
    real :: obs_var_s_oi = 0.0
    real :: std_bg = 0.0
    real :: std_oi_o = 0.0
    real :: std_oi_g = 0.0
    real :: dist, dist0
    real :: v2_h, v2_l
    real :: depth_bot

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

    !---------------------------------------------------------------------------

    id_eakf_total = mpp_clock_id('(ODA filter computation)')
    call mpp_clock_begin(id_eakf_total)

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
       write (stdout_unit, *) 'mean_inflate is ', mean_inflate, 'cov_inflate is ', cov_inflate
       write (stdout_unit, *) 'temp obs standard derivation is ', sigma_o_t
       write (stdout_unit, *) 'salt obs standard derivation is ', sigma_o_s
       write (stdout_unit, *) 'impact_levels is ', impact_levels
       write (stdout_unit, *) 'depth_cut is', depth_cut
       write (stdout_unit, *) 'e_flder_oer is', e_flder_oer
       write (stdout_unit, *) 'temp(salt)_to_salt(temp) is', temp_to_salt, salt_to_temp

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

    Prof => Profiles
    do while (associated(Prof))
       k0 = Prof%levels
       if(.not. associated(Prof%forecast)) allocate(Prof%forecast(k0))
       do kk = 1, k0
          Prof%forecast(kk) = missing_value
          if ( Prof%flag(kk) ) then ! add each level flag
             interp_flag = .true.
             v2_h = 0.0
             v2_l = 0.0
             if ( Prof%variable == TEMP_ID ) then
                do i_idx=1, 4
                   ind_temp_h = Prof%obs_def(kk)%state_var_index(i_idx)
                   ind_temp_l = Prof%obs_def(kk)%state_var_index(i_idx+4)
                   v2_h = v2_h + ens_mean(ind_temp_h)*Prof%obs_def(kk)%coef(i_idx)
                   v2_l = v2_l + ens_mean(ind_temp_l)*Prof%obs_def(kk)%coef(i_idx)
                   if(ens_mean(ind_temp_h) .eq. 0.0 .or. ens_mean(ind_temp_l) .eq. 0.0) &
                           interp_flag = .false.
                end do
             elseif ( Prof%variable == SALT_ID ) then
                do i_idx=1, 4
                   ind_salt_h = Prof%obs_def(kk)%state_var_index(i_idx)+salt_offset
                   ind_salt_l = Prof%obs_def(kk)%state_var_index(i_idx+4)+salt_offset
                   v2_h = v2_h + ens_mean(ind_salt_h)*Prof%obs_def(kk)%coef(i_idx)
                   v2_l = v2_l + ens_mean(ind_salt_l)*Prof%obs_def(kk)%coef(i_idx)
                   if(ens_mean(ind_salt_h) .eq. 0.0 .or. ens_mean(ind_salt_l) .eq. 0.0) &
                           interp_flag = .false.
                end do
             end if
             if(interp_flag) then
                Prof%forecast(kk) = v2_h*Prof%obs_def(kk)%coef(5) + v2_l*Prof%obs_def(kk)%coef(6)
             else
                Prof%flag(kk) = .false.
             endif
          end if
       end do
       Prof => Prof%cnext
    end do

    Prof => Profiles
    doloop_9: do while (associated(Prof)) ! (9)
       k0 = Prof%levels

       call get_time(Prof%tdiff, diff_seconds, diff_days)
       diff_hours = diff_seconds/3600 + diff_days * 24
       cov_factor_t = comp_cov_factor(diff_hours, update_window)

       if ( cov_factor_t > 0.0 ) then ! control 10d time window (5+ and 5-)
          obs_loc%lon = Prof%lon
          obs_loc%lat = Prof%lat

          call kd_search_radius(kdroot, obs_loc%lon, obs_loc%lat, &
                  2*temp_dist*cos(obs_loc%lat*DEG_TO_RAD), &
                  kd_ind, kd_dist, kd_num, .true.)

          doloop_8: do k=1, kd_num ! (8)
             ii_ens = lon1d(kd_ind(k)); jj_ens = lat1d(kd_ind(k))
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

             !! do not do assim over some islands
             !! East bound
             !if ( (model_loc%lat > 10.0 .and. model_loc%lat < 30.0) .and. &
                     !(model_loc%lon > 278.0 .and. model_loc%lon < 305.0) )  assim_flag = .false.
             !if ( (model_loc%lat >= 17.0 .and. model_loc%lat <= 25.0) .and. &
                     !(model_loc%lon > 262.0 .and. model_loc%lon <= 278.0) ) assim_flag = .false.
             !if ( (model_loc%lat > -27.0 .and. model_loc%lat < -18.0) .and. &
                     !(model_loc%lon > 300.0 .and. model_loc%lon < 330.0) ) assim_flag = .false.
             !! do not do assim over some islands
             !! West bound
             !if ( (model_loc%lat > -15.0 .and. model_loc%lat < 30.0) .and. &
                     !(model_loc%lon > 25.0 .and. model_loc%lon < 100.0) ) assim_flag = .false.

             ifblock_6: if ( assim_flag ) then ! (6)
                dist = kd_dist(k)
                cov_factor = cov_factor_t

                ifblock_4_4: if ( Prof%variable == TEMP_ID .or. &
                     & Prof%variable == SALT_ID ) then ! T,S -> T,S (4.4)
                   kk_bot = Prof%k_index(k0)
                   depth_bot = Prof%depth(k0)
                   doloop_4: do kk = 1, k0 ! (4)
                      if ( Prof%flag(kk) .and. &
                              Prof%depth(kk) <= depth_cut) then ! add each level flag
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
                               enso_temp(j_ens) = v2_h*Prof%obs_def(kk)%coef(5) &
                                       + v2_l*Prof%obs_def(kk)%coef(6)
                            end if
                            v2_h = 0.0
                            v2_l = 0.0
                            if ( Prof%variable == SALT_ID ) then
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
                            obs_sigma_t = sigma_o_t*exp(-Prof%depth(kk)/e_flder_oer)
                            obs_value = Prof%data(kk)
                            obs_var_t = obs_sigma_t * obs_sigma_t
                            call obs_increment_prf_eta_hyb(enso_temp, ens_size, obs_value, obs_var_t,&
                                 & obs_inc_eakf_temp, obs_var_t_oi, obs_inc_oi_temp, std_bg)
                         elseif ( Prof%variable == SALT_ID ) then
                            obs_sigma_s = sigma_o_s*exp(-Prof%depth(kk)/e_flder_oer)
                            obs_value = Prof%data(kk)
                            obs_var_s = obs_sigma_s * obs_sigma_s
                            call obs_increment_prf_eta_hyb(enso_salt, ens_size, obs_value, obs_var_s,&
                                 & obs_inc_eakf_salt, obs_var_s_oi, obs_inc_oi_salt, std_bg)
                         end if

                         kk0 = Prof%k_index(kk)
                         kk1 = kk0 - 2 * impact_levels +1
                         kk2 = kk0 + 2 * impact_levels
                         if(kk1 < 1) kk1 = 1
                         if(kk2 > nk) kk2 = nk

                         doloop_3: do kk_ens=kk1, kk2 ! (3)
                            t_tau   = (kk_ens-1)*blk + i_h
                            s_tau   = salt_offset + t_tau
                            if ( kk_ens <= kk0 ) then
                              diff_k = kk0 - kk_ens + 1
                            else
                              diff_k = kk_ens - kk0
                            end if
                            cov_factor_v = comp_cov_factor( diff_k, REAL(impact_levels) )
                            cov_factor = cov_factor_v * cov_factor_t

                            ifblock_2: if ( ens_mean(t_tau) /= 0.0 ) then ! (2) using ens_mean(t_tau) as mask
                               ifblock_1: if ( cov_factor > 0.0 ) then ! (1)
                                  ifblock_0_60: if ( Prof%variable == TEMP_ID ) then ! (0.60)
                                     ! using temperature adjusts temperature and salinity
                                     ! temperature itself
                                     if ( abs(obs_loc%lat) < 80.0 ) then
                                        dist0 = temp_dist*cos(obs_loc%lat*DEG_TO_RAD)
                                     else
                                        dist0 = temp_dist*cos(80.0*DEG_TO_RAD)
                                     end if

                                     cov_factor_h = comp_cov_factor(dist, dist0) * &
                                             cos((model_loc%lat-obs_loc%lat)*DEG_TO_RAD)
                                     cov_factor = cov_factor_h * cov_factor_v * cov_factor_t
                                     assim_var = 1
                                     ens_inc(:) = 0.0
                                     ifblock_0_60_1: if ( ens_mean(t_tau) /= 0.0 .and. cov_factor /= 0.0 ) then
                                        call update_from_obs_inc_prf_hyb(enso_temp, obs_inc_eakf_temp, &
                                                obs_inc_oi_temp, ens(t_tau, :), ens_size, ens_inc, cov_factor, &
                                                cor_oi, std_oi_o, std_oi_g, assim_var, flag_hyb)
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
                                     end if ifblock_0_60_1

                                     ! temperature impact salinity
                                     ifblock_0_50: if ( temp_to_salt ) then ! (0.5)
                                        if ( abs(obs_loc%lat) < 80.0 ) then
                                           dist0 = temp_dist*cos(obs_loc%lat*DEG_TO_RAD)
                                        else
                                           dist0 = temp_dist*cos(80.0*DEG_TO_RAD)
                                        end if

                                        cov_factor_h = comp_cov_factor(dist, dist0) * &
                                                cos((model_loc%lat-obs_loc%lat)*DEG_TO_RAD)
                                        cov_factor = cov_factor_h * cov_factor_v * cov_factor_t
                                        assim_var = 2
                                        ens_inc(:) = 0.0

                                        ifblock_0_50_1: if ( ens_mean(s_tau) /= 0.0 .and. cov_factor /= 0.0 ) then
                                           call update_from_obs_inc_prf_hyb(enso_temp, obs_inc_eakf_temp, &
                                                   obs_inc_oi_temp, ens(s_tau, :), ens_size, ens_inc, cov_factor, &
                                                   cor_oi, std_oi_o, std_oi_g, assim_var, flag_hyb)
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
                                        end if ifblock_0_50_1
                                     end if ifblock_0_50 ! impact salinity or not (0.5)
                                  
                                  elseif ( Prof%variable == SALT_ID ) then ! (0.60)
                                     if ( abs(obs_loc%lat) < 80.0 ) then
                                        dist0 = salt_dist*cos(obs_loc%lat*DEG_TO_RAD)
                                     else
                                        dist0 = salt_dist*cos(80.0*DEG_TO_RAD)
                                     end if

                                     cov_factor_h = comp_cov_factor(dist, dist0)*&
                                          & cos((model_loc%lat-obs_loc%lat)*DEG_TO_RAD)
                                     cov_factor = cov_factor_h * cov_factor_v * cov_factor_t
                                     assim_var = 1
                                     ens_inc(:) = 0.0
                                     ifblock_0_61_1: if ( ens_mean(s_tau) /= 0.0 .and. cov_factor /= 0.0 ) then
                                        call update_from_obs_inc_prf_hyb(enso_salt, obs_inc_eakf_salt, &
                                                obs_inc_oi_salt, ens(s_tau, :), ens_size, ens_inc, cov_factor, &
                                                cor_oi, std_oi_o, std_oi_g, assim_var, flag_hyb)
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
                                     end if ifblock_0_61_1

                                     ! salinity impact temperature
                                     ifblock_0_5: if ( salt_to_temp ) then ! (0.5)
                                        if ( abs(obs_loc%lat) < 80.0 ) then
                                           dist0 = salt_dist*cos(obs_loc%lat*DEG_TO_RAD)
                                        else
                                           dist0 = salt_dist*cos(80.0*DEG_TO_RAD)
                                        end if

                                        cov_factor_h = comp_cov_factor(dist, dist0)*&
                                             & cos((model_loc%lat-obs_loc%lat)*DEG_TO_RAD)
                                        cov_factor = cov_factor_h * cov_factor_v * cov_factor_t
                                        assim_var = 2
                                        ens_inc(:) = 0.0
                                        ifblock_0_5_1: if ( ens_mean(t_tau) /= 0.0 .and. cov_factor /= 0.0 ) then
                                           call update_from_obs_inc_prf_hyb(enso_salt, obs_inc_eakf_salt, &
                                                   obs_inc_oi_salt, ens(t_tau, :), ens_size, ens_inc, cov_factor, &
                                                   cor_oi, std_oi_o, std_oi_g, assim_var, flag_hyb)
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
                                        end if ifblock_0_5_1
                                     end if ifblock_0_5 ! impact temperature or not (0.5)
                                 end if ifblock_0_60 ! finish processing profiles (0.60)
                              end if ifblock_1 ! only use obs which has cov_factor > 0.0 (1)
                            end if ifblock_2 ! only do non-masked points (2)
                         end do doloop_3 ! finish adjustments for related model levels (3)
                      end if ! add each level flag
                   end do doloop_4 ! go through one profile column (4)
                end if ifblock_4_4 ! T,S -> T,S,U,V,tx,ty; ETA -> T,S,U,V (4.4)
             end if ifblock_6 ! get rid of land points (6)
          end do doloop_8 ! finish the adjustments for all related model columns (8)
       end if ! control 10d time window (5+ & 5-)
       Prof => Prof%cnext
    end do doloop_9 ! finish all profiles (9)

    !===== Eakf assim finish =====================================

    ! Compute ensemble mean of current assimilated state
    !
    ! Redistribute the sub ensemble state vector ens(:, :) back to the model grids
    ! in the local-domain.
    call mpp_sync_self()
    call red_ens(Posterior%T, Posterior%S, isd, ied, jsd, jed, halox, haloy, nk, ens)

    ens_mean = sum(ens, dim=2) / ens_size
    Prof => Profiles
    do while (associated(Prof))
       k0 = Prof%levels
       if(.not. associated(Prof%analysis)) allocate(Prof%analysis(k0))
       do kk = 1, k0
          Prof%analysis(kk) = missing_value
          if ( Prof%flag(kk) ) then ! add each level flag
             v2_h = 0.0
             v2_l = 0.0
             if ( Prof%variable == TEMP_ID ) then
                do i_idx=1, 4
                   ind_temp_h = Prof%obs_def(kk)%state_var_index(i_idx)
                   ind_temp_l = Prof%obs_def(kk)%state_var_index(i_idx+4)
                   v2_h = v2_h + ens_mean(ind_temp_h)*Prof%obs_def(kk)%coef(i_idx)
                   v2_l = v2_l + ens_mean(ind_temp_l)*Prof%obs_def(kk)%coef(i_idx)
                end do
             elseif ( Prof%variable == SALT_ID ) then
                do i_idx=1, 4
                   ind_salt_h = Prof%obs_def(kk)%state_var_index(i_idx)+salt_offset
                   ind_salt_l = Prof%obs_def(kk)%state_var_index(i_idx+4)+salt_offset
                   v2_h = v2_h + ens_mean(ind_salt_h)*Prof%obs_def(kk)%coef(i_idx)
                   v2_l = v2_l + ens_mean(ind_salt_l)*Prof%obs_def(kk)%coef(i_idx)
                end do
             end if
             Prof%analysis(kk) = v2_h*Prof%obs_def(kk)%coef(5) + v2_l*Prof%obs_def(kk)%coef(6)
         end if
       end do
       Prof => Prof%cnext
    end do

    if ( debug_eakf ) then
       write (UNIT=stdout_unit, FMT='("PE ",I5,": finished red_ens")') pe()
    end if

    first_run_call = .false.

    call mpp_clock_end(id_eakf_total)
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
       allocate(enso_temp(ens_size), STAT=istat);               enso_temp = 0.0
       allocate(obs_inc_eakf_temp(ens_size), STAT=istat);       obs_inc_eakf_temp= 0.0
       allocate(obs_inc_oi_temp(ens_size), STAT=istat);         obs_inc_oi_temp = 0.0
       allocate(ens_inc(ens_size), STAT=istat);                 ens_inc = 0.0
       allocate(enso_salt(ens_size), STAT=istat);               enso_salt = 0.0
       allocate(obs_inc_eakf_salt(ens_size), STAT=istat);       obs_inc_eakf_salt = 0.0
       allocate(obs_inc_oi_salt(ens_size), STAT=istat);         obs_inc_oi_salt = 0.0
       allocate(lon1d(blk), lat1d(blk) )
       allocate(kd_ind(blk) )
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

    integer :: m, i, j, k, idx, blk_h

    !  Get initial state for ensemble members
    blk_h = (ied-isd+1)*(jed-jsd+1)

    do m=1, size(as, 2)
       do k=1, nk
          do j=jsd+haloy, jed-haloy ! drop halo
             do i=isd+halox, ied-halox ! drop halo
                idx = (k-1)*blk_h + (j-jsd)*(ied-isd+1) + i-isd+1
                T_ensemble(i,j,k,m) = as(idx, m)
                idx = blk_h*nk + idx
                S_ensemble(i,j,k,m) = as(idx, m)
             end do
          end do
       end do
    end do
  end subroutine red_ens

end module eakf_oda_mod
