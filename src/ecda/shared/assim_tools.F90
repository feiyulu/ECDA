module assim_tools_mod
  ! A variety of operations required by assimilation.

  use fms_mod, only : file_exist, open_namelist_file, check_nml_error, write_version_number, close_file
  use fms_mod, only : error_mesg, WARNING, FATAL
  use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin, mpp_clock_end, stdlog, stdout

  logical :: first_run_call = .true.

  private
  public assim_tools_init, obs_increment, update_from_inc

  !---- namelist ASSIM_TOOLS_NML with default values

  real :: cor_cutoff = 0.0

  namelist /assim_tools_nml/ cor_cutoff

  !---- module name and version number
  character(len=*), parameter :: MODULE_NAME = 'assim_tools_mod'
  character(len=*), parameter :: VERS_NUM = '$Id$'

contains

  subroutine assim_tools_init()

    integer :: unit, istat, stdlog_unit, stdout_unit

    ! Read namelist for run time control
    if ( file_exist('input.nml') ) then
       unit = open_namelist_file()
       read(UNIT=unit, NML=assim_tools_nml, IOSTAT=istat)
       call close_file(unit)
    else
       ! Set istat to an arbitrary positive number if input.nml does not exist
       istat = 100
    end if

    if ( check_nml_error(istat, 'assim_tools_nml') < 0 ) then
       call error_mesg('assim_tools_mod::assim_tools_init', 'ASSIM_TOOLS_NML not found in input.nml,  Using defaults.', WARNING)
    end if


    if ( mpp_pe() == mpp_root_pe() .and. first_run_call ) then
       stdout_unit = stdout()
       stdlog_unit = stdlog()
       call write_version_number(VERS_NUM, MODULE_NAME)
       write (UNIT=stdlog_unit, NML=assim_tools_nml)
       !write (UNIT=stdout_unit, NML=assim_tools_nml)
    end if

    first_run_call = .false.
  end subroutine assim_tools_init

  subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc, obs_dist)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: ens
    real, intent(in) :: obs, obs_var
    real, intent(inout), dimension(ens_size) :: obs_inc
    real, intent(out) :: obs_dist

    real :: a, cov
    real :: mean, new_cov, new_mean
    integer :: ie

    mean = sum(ens) / ens_size
    cov = 0.0

    do ie=1, ens_size
       cov = cov + (ens(ie)-mean)*(ens(ie)-mean)
    end do
    cov = cov / (ens_size-1)

    if ( cov == 0.0 ) return
    obs_dist = abs(obs-mean)

    new_cov = (cov *obs_var)/(cov + obs_var)
    new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
    a = sqrt(new_cov/cov)
    obs_inc = a * (ens - mean) + new_mean - ens

  end subroutine obs_increment

  subroutine update_from_inc(obs, obs_inc, state, ens_size, new_state, cov_factor, assim_var)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: obs, obs_inc, state
    real, intent(inout), dimension(ens_size) :: new_state
    real, intent(in) :: cov_factor
    integer, intent(in) :: assim_var

    real :: mean_s, mean_o, cv_s, cv_o, cv, cr
    real :: std_o, std_g
    real :: cor_cut

    integer :: ie

    ! Compute statistics up to the second moment
    cor_cut = cor_cutoff

    mean_s = sum(state) / ens_size
    mean_o = sum(obs) / ens_size

    cv_s = 0.0
    cv_o = 0.0
    cv   = 0.0

    do ie=1, ens_size
       cv_s = cv_s + (state(ie) - mean_s)*(state(ie) - mean_s)
       cv_o = cv_o + (obs(ie) - mean_o)*(obs(ie) - mean_o)
       cv   = cv + (state(ie) - mean_s)*(obs(ie) - mean_o)
    end do

    cv_s = cv_s / (ens_size-1)
    cv_o = cv_o / (ens_size-1)
    cv   = cv / (ens_size-1)
    if ( cv_s /= 0.0 .and. cv_o /= 0.0 ) then
       std_o = sqrt(cv_o)
       std_g = sqrt(cv_s)
       cr   = cv /(std_g*std_o)
    else
       std_o = 0.0
       std_g = 0.0
       cr = 0.0
    end if

    new_state = 0.0

    if ( abs(cr) >= cor_cut ) then
       new_state = (cov_factor * cv / cv_o) * obs_inc
    else
       new_state = 0.0
    end if
  end subroutine update_from_inc

end module assim_tools_mod
