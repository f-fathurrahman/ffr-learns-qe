!================================================================================
!
! Modules:
!
! (1) version_m Originally By FHJ
!
! Returns information on GIT repository date and commit hash.
!
!================================================================================

module version_m
  implicit none
  private
  character(len=*), parameter :: version_base = "3.1"
  logical, parameter :: is_git = .false.
  character(len=*), parameter :: git_date = ""
  character(len=*), parameter :: git_commit = ""
  character(len=*), parameter :: git_author = ""
  public :: get_version_info
contains
  subroutine get_version_info(ainfo)
    character(len=256), intent(out) :: ainfo
    if (is_git) then
      ainfo = 'version ' // trim(version_base) // '+git ' // &
              trim(git_commit) // ' by ' // trim(git_author) // ' at ' // trim(git_date)
    else
      ainfo = 'version ' // trim(version_base)
    endif
  end subroutine get_version_info
end module version_m
