module eakf_tab_mod

  ! this module is produced by snz for initializing EAKF

  USE fms2_io_mod,    ONLY : ascii_read
  USE fms2_io_mod,    ONLY : FmsNetcdfFile_t
  USE fms2_io_mod,    ONLY : open_file, close_file, read_data
  USE fms2_io_mod,    ONLY : get_num_variables, get_num_dimensions, variable_exists

  use cov_cutoff_mod, only : comp_cov_factor
  use fms_mod,        only : check_nml_error
  use fms_mod,        only : write_version_number, error_mesg, WARNING, FATAL
  USE mpp_mod,        ONLY : input_nml_file
  use mpp_mod,        only : stdlog, mpp_pe, mpp_root_pe
  use constants_mod,  only : DEG_TO_RAD

  implicit none

  real, dimension(60) :: rrt
  real, dimension(200) :: rrv
  real, dimension(291) :: rcos, vcos
  real, dimension(600) :: rrr, rrh

  real, dimension(:,:), allocatable :: sgm_e, sgm_tx, sgm_ty
  real, dimension(:,:), allocatable :: cor_etx, cor_ety

  real, dimension(:,:,:), allocatable :: sgm_t, sgm_s, sgm_u, sgm_v
  real, dimension(:,:,:), allocatable :: cor_et, cor_es
  real, dimension(:,:,:), allocatable :: cor_eu, cor_ev

  real, dimension(:,:,:,:), allocatable :: cor_sta

  character(len=*), parameter :: MODULE_NAME = 'eakf_tab_mod'
  character(len=*), parameter :: VERS_NUM = '$Id$'

  public eakf_init, rrr, rrh, rrv, rrt, rcos, vcos
  public sgm_e, sgm_tx, sgm_ty, sgm_t, sgm_s, sgm_u, sgm_v, cor_etx, &
       & cor_ety, cor_et, cor_es, cor_eu, cor_ev, cor_sta

contains

  subroutine eakf_init()

    real, parameter :: CC = 1.25e-6
    real, parameter :: LL = 1.25e6
    real, parameter :: LV = 1./LL
    real, parameter :: LCV = 1./(LL*CC)

    real :: rr
    real, dimension(41,41,11) :: cor

    integer :: i, unit, j, k, k0, j0, j00, l0, l, blk_j0, blk_j1, blk_jk0
    integer :: ierr, istat
    integer :: stdlog_unit

    integer :: nlon, nlat, nlev
    integer :: ndim, nvar, natt, ntime

    character(len=356) :: emsg_local
    character(len=80) :: sgm_fname, cor_fname, cor_sta_fname
    character(len=32) :: axisname

    type(FmsNetcdfFile_t) :: fileobj

    !---- namelist with default values
    real :: cutoff = 2000.0e3
    real :: cutoff_v = 50.0
    real :: cutoff_t = 5.0

    namelist /eakf_tab_nml/ cutoff, cutoff_v, cutoff_t

    stdlog_unit = stdlog()

    ! Read namelist for run time control
    read (input_nml_file, nml=eakf_tab_nml, iostat=istat)
    ierr = check_nml_error(IOSTAT=istat,NML_NAME="EAKF_TAB_NML")

    if ( ierr < 0 ) then
       if ( mpp_pe() == mpp_root_pe() ) then
          call error_mesg('eakf_tab_mod::eakf_init', 'EAKF_TAB_NML not found in input.nml.  Using defaults.', WARNING)
       end if
    end if

    call write_version_number(VERS_NUM, MODULE_NAME)
    write (UNIT=stdlog_unit, NML=eakf_tab_nml)

    do i=1, 600
       rr = 5000.0 * (i-1)
       rrh(i) = comp_cov_factor(rr, cutoff)
       rrr(i) = (cos(CC*rr) + sin(CC*rr)*LCV) * exp(-rr*LV)
    end do
    do i=1, 200
       rr = 2.0 * (i-1)
       rrv(i) = comp_cov_factor(rr, cutoff_v)
    end do
    do i=1, 60
       rr = 12.0 * (i-1)
       rrt(i) = comp_cov_factor(rr, cutoff_t*24.)
    end do

    do i=2, 291
       if ( i <= 201 ) then
          rr = (i-1) * 0.05 * DEG_TO_RAD
       else if ( i <= 221 ) then
          rr = ((i-201)*0.5 + 10.0) * DEG_TO_RAD
       else
          rr = ((i-221) + 20.0) * DEG_TO_RAD
       end if

       rcos(i) = cos(rr)

       if ( rcos(i) > 1.e-8 ) then
          vcos(i) = 1./rcos(i)
       else
          vcos(i) = 1.e8
       end if
    end do

    rcos(1) = 1.0
    vcos(1) = 1.0

    ! read out oi sigma and cor.
    sgm_fname ='STA_COV/sgm.nc'
    if ( open_file(fileobj, trim(sgm_fname), "read", is_restart = .false.)) then 
       ndim= get_num_dimensions(fileobj)
       nvar= get_num_variables(fileobj)
       call get_dimension_size(fileobj,"XT_OCEAN",nlon)
       call get_dimension_size(fileobj,"YT_OCEAN",nlat)
       call get_dimension_size(fileobj,"ZT_OCEAN",nlev)
    endif

    if ( nlon /= 360 .or. nlat /= 200 .or. nlev /= 50 ) then
       write (emsg_local, '("nlon(",I4,") .ne. 360 or  nlat(",I4,") .ne. 200 or nlev(",I4,") .ne 50")' ) nlon, nlat, nlev
       call error_mesg('eakf_tab_mod::eakf_init', trim(emsg_local)//' (sgm.nc)', FATAL)
    end if

    ! Allocate arrays for sgm.nc data
    allocate(sgm_e(nlon,nlat), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate sgm_e.', FATAL)
    end if
    allocate(sgm_tx(nlon,nlat), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate sgm_tx.', FATAL)
    end if
    allocate(sgm_ty(nlon,nlat), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate sgm_ty.', FATAL)
    end if
    allocate(sgm_t(nlon,nlat,nlev), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate sgm_t.', FATAL)
    end if
    allocate(sgm_s(nlon,nlat,nlev), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate sgm_s.', FATAL)
    end if
    allocate(sgm_u(nlon,nlat,nlev), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate sgm_u.', FATAL)
    end if
    allocate(sgm_v(nlon,nlat,nlev), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate sgm_v.', FATAL)
    end if

    call read_data(fileobj, "SGM_E", sgm_e)
    call read_data(fileobj, "SGM_TX", sgm_tx)
    call read_data(fileobj, "SGM_TY", sgm_ty)
    call read_data(fileobj, "SGM_T", sgm_t)
    call read_data(fileobj, "SGM_S", sgm_s)
    call read_data(fileobj, "SGM_U", sgm_u)
    call read_data(fileobj, "SGM_v", sgm_v)
    call close_file(fileobj)
    ! finish sgm read

    cor_fname ='STA_COV/cor.nc'
    if ( open_file(fileobj, trim(cor_fname), "read", is_restart = .false.)) then 
       ndim= get_num_dimensions(fileobj)
       nvar= get_num_variables(fileobj)
       call get_dimension_size(fileobj,"XT_OCEAN",nlon)
       call get_dimension_size(fileobj,"YT_OCEAN",nlat)
       call get_dimension_size(fileobj,"ZT_OCEAN",nlev)
    endif

    if ( nlon /= 360 .or. nlat /= 200 .or. nlev /= 50 ) then
       write (emsg_local, '("nlon(",I4,") .ne. 360 or  nlat(",I4,") .ne. 200 or nlev(",I4,") .ne 50")' ) nlon, nlat, nlev
       call error_mesg('eakf_tab_mod::eakf_init', 'nlon .ne. 360 or nlat .ne. 200 or nlev .ne. 50. (cor.nc)', FATAL)
    end if

    ! Allocate cor variables
    allocate(cor_etx(nlon,nlat), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate cor_etx.', FATAL)
    end if
    allocate(cor_ety(nlon,nlat), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate cor_ety.', FATAL)
    end if
    allocate(cor_et(nlon,nlat,nlev), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate cor_et.', FATAL)
    end if
    allocate(cor_es(nlon,nlat,nlev), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate cor_es.', FATAL)
    end if
    allocate(cor_eu(nlon,nlat,nlev), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate cor_eu.', FATAL)
    end if
    allocate(cor_ev(nlon,nlat,nlev), STAT=istat)
    if ( istat .ne. 0 ) then
       call error_mesg('eakf_tab_mod::eakf_init', 'Unable to allocate cor_ev.', FATAL)
    end if

    call read_data(fileobj, "COR_SSH_TX", cor_etx)
    call read_data(fileobj, "COR_SSH_TY", cor_ety)
    call read_data(fileobj, "COR_SSH_T",  cor_et)
    call read_data(fileobj, "COR_SSH_S",  cor_es)
    call read_data(fileobj, "COR_SSH_U",  cor_eu)
    call read_data(fileobj, "COR_SSH_V",  cor_ev)
    call close_file(fileobj)
    ! finish cor read

    do j=1, 200
       do i=1, 360
          if ( sgm_e(i,j) < -9.9 ) sgm_e(i,j) = 0.0
          if ( sgm_tx(i,j) < -9.9 ) sgm_tx(i,j) = 0.0
          if ( sgm_ty(i,j) < -9.9 ) sgm_ty(i,j) = 0.0
          if ( cor_etx(i,j) < -9.9 ) cor_etx(i,j) = 0.0
          if ( cor_ety(i,j) < -9.9 ) cor_ety(i,j) = 0.0
       end do
    end do
    do k=1, 50
       do j=1, 200
          do i=1, 360
             if(sgm_t(i,j,k) < -9.9 ) sgm_t(i,j,k) = 0.0
             if(sgm_s(i,j,k) < -9.9 ) sgm_s(i,j,k) = 0.0
             if(sgm_u(i,j,k) < -9.9 ) sgm_u(i,j,k) = 0.0
             if(sgm_v(i,j,k) < -9.9 ) sgm_v(i,j,k) = 0.0
             if(cor_et(i,j,k) < -9.9 ) cor_et(i,j,k) = 0.0
             if(cor_es(i,j,k) < -9.9 ) cor_es(i,j,k) = 0.0
             if(cor_eu(i,j,k) < -9.9 ) cor_eu(i,j,k) = 0.0
             if(cor_ev(i,j,k) < -9.9 ) cor_ev(i,j,k) = 0.0
          end do
       end do
    end do

    ! start to read cor_sta for in situ obs here
    allocate(cor_sta(41,41,11,17400))
    blk_j0 = 14*15
    blk_j1 = 15*15
    blk_jk0 = 14*15*4
    do l0  = 1, 8
       select case(l0)
       case (1)
          cor_sta_fname = 'STA_COV/cor_TT_P'
          j00 = 14
       case (2)
          cor_sta_fname = 'STA_COV/cor_SS_P'
          j00 = 14
       case (3)
          cor_sta_fname = 'STA_COV/cor_TS_P'
          j00 = 14
       case (4)
          cor_sta_fname = 'STA_COV/cor_ST_P'
          j00 = 14
       case (5)
          cor_sta_fname = 'STA_COV/cor_TT_A'
          j00 = 15
       case (6)
          cor_sta_fname = 'STA_COV/cor_SS_A'
          j00 = 15
       case (7)
          cor_sta_fname = 'STA_COV/cor_TS_A'
          j00 = 15
       case (8)
          cor_sta_fname = 'STA_COV/cor_ST_A'
          j00 = 15
       end select

       do k0=1, 15
          do j0=1, j00
             ! Append to filename
             write (cor_sta_fname, '(A16,I2.2,A1,I2.2,A3)') cor_sta_fname, j0, 'D', k0, '.nc'

             if ( open_file(fileobj, trim(cor_sta_fname), "read", is_restart = .false.)) then 
                ndim= get_num_dimensions(fileobj)
                nvar= get_num_variables(fileobj)
                call get_dimension_size(fileobj,"XT_OCEAN",nlon)
                call get_dimension_size(fileobj,"YT_OCEAN",nlat)
                call get_dimension_size(fileobj,"ZT_OCEAN",nlev)
             endif
             if ( nlon /= 41 .or. nlat /= 41 .or. nlev /= 11 ) then
                write (emsg_local, '("nlon(",I4,") .ne. 41 or  nlat(",I4,") .ne. 41 or nlev(",I4,") .ne. 11")' ) nlon, nlat, nlev
                call error_mesg('eakf_tab_mod::eakf_init', trim(emsg_local), FATAL)
             end if

             call read_data(fileobj, "COR", cor)

             if ( l0 <= 4 ) then
                l = (l0-1)*blk_j0 + (k0-1)*j00 + j0
             else
                l = blk_jk0 + (l0-5)*blk_j1 + (k0-1)*j00 + j0
             end if

             do k=1, 11
                do j=1, 41
                   do i=1, 41
                      if ( cor(i,j,k) < -9.9 ) cor(i,j,k) = 0.0
                      cor_sta(i,j,k,l) = cor(i,j,k)
                   end do
                end do
             end do

             call close_file(fileobj)

          end do
       end do
    end do
  end subroutine eakf_init
end module eakf_tab_mod
