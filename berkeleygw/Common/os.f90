!>===================================================================
!!
!! Module os_m Originally by FHJ Last Modified Jan/2018 (FHJ)
!!
!! Also, see some of the corresponding C routines in os_c.c
!!
!! Routines:
!!
!! - system
!! - sleep
!! - c_to_f_string
!! - get_host_name
!!
!!===================================================================
!The following macro puts any point/array in the [-0.5, 0.5) range:
!The following macro puts any point/array in the [0, 1) range:
!Integer division of a/b rounded up*/
!Rounds a up to the smallest multiple of b*/
! disable Fortran OMP pragmas if not -DOMP*/
! note: C standard does not permit $ in identifiers, however this seems acceptable
! as an extension, for all versions of cpp I tried. --DAS
! truncate spaces in string
!#!define TRUNC(s) trim(adjustl(s))
! Oracle compiler has a length limit of 132 characters and won`t support these macros
! No checking for faster performance, if not in debug mode
! Use this instead of the intrinsic 'deallocate' for pointers
! Use this instead of the intrinsic 'deallocate' for arrays
!the TOSTRING macro converts a macro into a string
! deprecated identifiers
! Created Sept 2011 by DAS.
! Define characteristics of various compilers, via compiler symbols (e.g. -DGNU)
! to be used directly from the arch.mk files, and then defining what we need to do
! for that compiler via the symbols for various properties (e.g. NOSIZEOF).
! Ideally, to support a new compiler, one need only change this file, adding a
! new block to define what -DNEWCOMPILER would mean.
! NOTE: of course, Makefile-level issues still need to be handled in common-rules.mk
! very ancient version may require NOSIZEOF
! FHJ: Support for Open64 will be removed shortly in favor of OpenUH
! open64 is very similar to path, it is an open-sourced version of it
! omp_lib.f90 needed to do OpenMP, see common-rules.mk.
! cce 7.4.4 and before support sizeof for intrinsic types, but need NOSIZEOF_TYPE
! cce 8.0.0 and later do not allow sizeof for multidimensional arrays, requiring us
! to turn sizeof off everywhere. Why would Cray do this?
! It is considered a bug in OPEN64 that sizeof will not work in our code.
! on some platforms there is a different return value for sizeof if build is 64-bit
! Intrinsic module for OpenMP. Almost all compilers that support OpenMP provide
! a "omp_lib.mod" module, though the OpenMP standard allow them to only ship a
! "omp_lib.h" Fortran header.
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
module os_m
  use, intrinsic :: iso_c_binding
  use nrtype_m
  implicit none
  private
  public :: &
    system, &
    sleep, &
    c_to_f_string, &
    get_host_name
contains
  !> Call command `command` using call to POSIX `system` routine.
  subroutine system(command)
    character(len=*), intent(in) :: command
    character(len=1,kind=c_char) :: command_c(len(command)+1)
    integer(c_int) :: ret
    integer :: ii
    interface !int system(const char *command);
      integer(c_int) function system_c(command) bind(c, name='system')
        use, intrinsic :: iso_c_binding, only: c_int, c_char
        character(len=1,kind=c_char), intent(in) :: command(*)
      end function system_c
    end interface
    do ii = 1, len(command)
      command_c(ii) = command(ii:ii)
    end do
    command_c(len(command)+1:len(command)+1) = c_null_char
    ret = system_c(command_c)
  end subroutine system
  !> Sleep for `time` seconds. Internally calls POSIX `nanosecond`.
  subroutine sleep(time)
    real(DP), intent(in) :: time
    integer :: ret
    interface !int sleep_c(double s); // defined in os_c.c
      integer(c_int) function sleep_c(s) bind(c)
        use, intrinsic :: iso_c_binding, only: c_int, c_double
        real(c_double), value, intent(in) :: s ! pass by value, not by reference
      end function sleep_c
    end interface
    ret = sleep_c(time)
  end subroutine sleep
  !> Convert a C string to a Fotran string.
  function c_to_f_string(c_string) result(f_string)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    character(kind=c_char,len=1), intent(in) :: c_string(*)
    character(len=:), allocatable :: f_string
    integer(c_size_t) :: n
    interface !size_t strlen(const char *s);
      integer(c_size_t) function strlen(s) bind(c)
        use, intrinsic :: iso_c_binding, only: c_char, c_size_t
        character(len=1,kind=c_char), intent(in) :: s(*)
      end function strlen
    end interface
    n = strlen(c_string)
    allocate(character(len=n) :: f_string)
    f_string = transfer(c_string(1:n), f_string)
  end function c_to_f_string
  !> Get the host name for the current processor / MPI rank.
  function get_host_name() result(hostname)
    character(len=:), allocatable :: hostname
    integer, parameter :: max_hostname_len=1024
    character(len=max_hostname_len), target :: hostname_
    integer :: namelen, mpierr
    integer(c_int) :: c_status
    character(kind=c_char,len=1) :: c_hostname(max_hostname_len)
    integer(c_size_t) :: c_namelen
    interface !int gethostname(char *name, size_t namelen);
      integer(c_int) function gethostname(name, namelen) bind(c)
        use, intrinsic :: iso_c_binding, only: c_char, c_int, c_size_t
        integer(c_size_t), value, intent(in) :: namelen
        character(len=1,kind=c_char), dimension(namelen), intent(inout) :: name
      end function gethostname
    end interface
    ! FHJ: Otherwise use system call to gethostname
    c_namelen = int(max_hostname_len, kind=c_size_t)
    c_status = gethostname(c_hostname, c_namelen)
    hostname = c_to_f_string(c_hostname)
  end function get_host_name
end module os_m
