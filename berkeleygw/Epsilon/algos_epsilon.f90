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
! Module to control selection of various accelerated algos in Epsilon
! Any build-specific preprocessor directives should be as high-level
! as possible, i.e. no explicit framework-based calls
module algos_epsilon_m
  use ISO_C_BINDING
  use algos_common_m
  use message_m
  use peinfo_m, only: peinf
  use push_pop_m
  implicit none
  private
  ! Publicly exposed members
  public :: CPU_ALGO, OPENACC_ALGO, OMP_TARGET_ALGO
  public :: output_algos
  public :: set_algos_to_cpu
  public :: set_algos_to_best_available_gpu
  public :: verify_gpu_settings
  public :: initialize_gpu
  public :: algos_inread
  public :: use_gpu
  ! General control variables
  integer(kind(CPU_ALGO)), public :: mtxel_algo
  integer(kind(CPU_ALGO)), public :: chi_summation_algo
contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! __ ____ _ ___ __ _ !
  ! / / / / /_(_) (_) /_(_)__ _____ !
  ! / / / / __/ / / / __/ / _ \/ ___/ !
  ! / /_/ / /_/ / / / /_/ / __(__ ) !
  ! \____/\__/_/_/_/\__/_/\___/____/ !
  ! !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine output_algos()
   
    if (peinf%inode.eq.0) then
      call output_algos_common()
      write(6,*)
      write(6,'(1x,a)') 'Algorithms used:'
      write(6,'(1x,a,a)') '- mtxel         : ', trim(get_algo_str(mtxel_algo))
      write(6,'(1x,a,a)') '- chi_summation : ', trim(get_algo_str(chi_summation_algo))
      write(6,*)
    endif
   
  end subroutine
  subroutine set_algos_to_cpu()
   
    call set_algos_to_cpu_common()
    mtxel_algo = CPU_ALGO
    chi_summation_algo = CPU_ALGO
   
  end subroutine
  subroutine set_algos_to_best_available_gpu()
   
    call set_algos_to_best_available_gpu_common()
   
  end subroutine
  subroutine verify_gpu_settings()
   
    ! the main body of the code.
    call verify_gpu_settings_common()
    ! First check if the desired algorithm is compiled in
    if (mtxel_algo == OPENACC_ALGO .or. &
        chi_summation_algo == OPENACC_ALGO) then
      call die_algos("OpenACC")
    end if
    if (mtxel_algo == OMP_TARGET_ALGO .or. &
        chi_summation_algo == OMP_TARGET_ALGO) then
      call die_algos("OpenMP Target")
    end if
    ! Then check if we are using the GPU
    if (.not. use_gpu .and. (mtxel_algo /= CPU_ALGO .or. &
                             chi_summation_algo /= CPU_ALGO)) then
      call die("You have specified one or more GPU-accelerated algorithms, but the GPU is not enabled.", only_root_writes=.true.)
    end if
    ! Then check for unimplemented algorithms
    if (mtxel_algo == OPENACC_ALGO) then
      call die("You have specified the OpenACC version of mtxel(), but this has not been implemented.", only_root_writes=.true.)
    end if
    if (mtxel_algo == OMP_TARGET_ALGO) then
      call die("You have specified the OpenMP version of mtxel(), but this has not been implemented.", only_root_writes=.true.)
    end if
    if (chi_summation_algo == OMP_TARGET_ALGO) then
      call die("You have specified the OpenMP version of chi_summation(), but this has not been implemented.", only_root_writes=.true.)
    end if
    ! Finally, check for algorithms that don not play nicely with one another
    ! None so far...
   
  end subroutine
  subroutine initialize_gpu(algo_type)
    integer(kind(CPU_ALGO)), intent(in) :: algo_type
   
    call initialize_gpu_common(algo_type)
   
  end subroutine initialize_gpu
  subroutine algos_inread(keyword, line, found)
    character(len=*), intent(in) :: keyword, line
    logical, intent(out) :: found
    character*256 :: errmsg
   
    call algos_inread_common(keyword, line, found)
    if (.not. found) then
      found = .true.
      if (trim(keyword) .eq. 'use_gpu') then
        if (trim(line) .eq. '.true.') then
          write (errmsg,'(a)') 'GPU acceleration not compiled into this executable, use_gpu keyword is invalid.'
          call die(errmsg, only_root_writes = .true.)
        elseif(trim(line) .eq. '.false.') then
          call set_algos_to_cpu()
        else
          write(errmsg, '(3a)') 'Unexpected parameter "', trim(line), '" for keyword "use_gpu" was found in kernel.inp.'
          call die(errmsg, only_root_writes = .true.)
        end if
      else if (trim(keyword) .eq. 'mtxel_algo') then
        mtxel_algo = get_algo(trim(line))
      else if (trim(keyword) .eq. 'chi_summation_algo') then
        chi_summation_algo = get_algo(trim(line))
      else
        found = .false.
      end if
    end if
   
  end subroutine algos_inread
end module algos_epsilon_m
