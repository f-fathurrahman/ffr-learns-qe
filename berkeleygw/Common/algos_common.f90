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
! Module to control selection of various accelerated algos in BerkeleyGW
! Any build-specific preprocessor directives should be as high-level
! as possible, i.e. no explicit framework-based calls
module algos_common_m
  use ISO_C_BINDING
  use message_m
  use peinfo_m, only: peinf
  use push_pop_m
  implicit none
  private
  ! Publicly exposed members
  public :: get_algo
  public :: get_algo_str
  public :: output_algos_common
  public :: set_algos_to_cpu_common
  public :: set_algos_to_best_available_gpu_common
  public :: verify_gpu_settings_common
  public :: initialize_gpu_common
  public :: algos_inread_common
  public :: die_algos
  public :: CPU_ALGO, OMP_TARGET_ALGO, OPENACC_ALGO
  enum, bind(c)
    enumerator :: CPU_ALGO
    enumerator :: OMP_TARGET_ALGO
    enumerator :: OPENACC_ALGO
  end enum
  ! General control variables
  logical, public :: use_gpu
contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! __ ____ _ ___ __ _ !
  ! / / / / /_(_) (_) /_(_)__ _____ !
  ! / / / / __/ / / / __/ / _ \/ ___/ !
  ! / /_/ / /_/ / / / /_/ / __(__ ) !
  ! \____/\__/_/_/_/\__/_/\___/____/ !
  ! !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function get_algo(algo_str) result(algo)
    character(len=*), intent(in) :: algo_str
    integer(kind(CPU_ALGO)) :: algo
    character(len=1024) :: errmsg
   
    select case (trim(algo_str))
    case('CPU_ALGO')
      algo = CPU_ALGO
    case('OMP_TARGET_ALGO')
      algo = OMP_TARGET_ALGO
    case('OPENACC_ALGO')
      algo = OPENACC_ALGO
    case default
      write(errmsg,'(3a)') 'Invalid version "', trim(algo_str), '" has been specified for one or more algorithms'
      call die(errmsg, only_root_writes=.true.)
    end select
   
  end function
  function get_algo_str(algo) result(algo_str)
    integer(kind(CPU_ALGO)), intent(in) :: algo
    character(len=:), allocatable :: algo_str
    character(len=1024) :: errmsg
   
    select case (algo)
    case(CPU_ALGO)
      algo_str = 'CPU_ALGO'
    case(OMP_TARGET_ALGO)
      algo_str = 'OMP_TARGET_ALGO'
    case(OPENACC_ALGO)
      algo_str = 'OPENACC_ALGO'
    case default
      write(errmsg,'(a)') 'Argument to get_algo_str() is invalid, this is an internal error'
      call die(errmsg, only_root_writes=.true.)
    end select
   
  end function
  subroutine output_algos_common()
   
    if (peinf%inode.eq.0) then
      if (use_gpu) then
        write(6,'(1x,a)') 'GPU acceleration is : ENABLED'
      else
        write(6,'(1x,a)') 'GPU acceleration is : DISABLED'
      end if
    endif
   
  end subroutine
  subroutine set_algos_to_cpu_common()
   
    use_gpu = .false.
   
  end subroutine
  subroutine set_algos_to_best_available_gpu_common()
   
    use_gpu = .true.
   
  end subroutine
  subroutine verify_gpu_settings_common()
   
   
  end subroutine
  subroutine initialize_gpu_common(algo_type)
    integer(kind(CPU_ALGO)), intent(in) :: algo_type
   
    select case (algo_type)
    case (OPENACC_ALGO)
      call die_algos("OpenACC")
    case (OMP_TARGET_ALGO)
      call die_algos("OpenMP Target")
    case default
      call die("Invald algorithm for initialize_gpu_common", &
               only_root_writes = .true.)
    end select
   
  end subroutine initialize_gpu_common
  subroutine algos_inread_common(keyword, line, found)
    character(len=*), intent(in) :: keyword, line
    logical, intent(out) :: found
   
    found = .false.
   
  end subroutine algos_inread_common
  subroutine die_algos(algo)
    character(len=*), intent(in) :: algo
    character*256 :: errmsg
   
    write(errmsg, '(5a)') "You have specified one or more ", trim(algo), &
               " algorithms, but ", trim(algo), &
               " is not compiled into this executable."
    call die(errmsg, only_root_writes=.true.)
   
  end subroutine die_algos
end module algos_common_m
