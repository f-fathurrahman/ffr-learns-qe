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
module acc_linalg_m
  use ISO_C_BINDING
  use algos_common_m
  use blas_m
  use message_m, only: die
  use nrtype_m
  implicit none
  private
  public :: acc_dgemm
  public :: acc_xgemm
  public :: acc_zgemm
  ! WPH: Implementation of dgemm as a generic subroutine
  interface acc_xgemm
    module procedure acc_zgemm
    module procedure acc_dgemm
  end interface
contains
  subroutine acc_zgemm(transa, transb, &
                       m, n, k, &
                       alpha, a, lda, b, ldb, &
                       beta, c, ldc, &
                       algo)
    implicit none
    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double complex :: alpha, beta
    double complex :: a(lda, *), b(ldb, *), c(ldc, *)
    integer(kind(CPU_ALGO)) :: algo
    select case (algo)
    case (OPENACC_ALGO)
      call die("OpenACC version of acc_zgemm requested, but OpenACC not "&
            &"compiled into this executable", only_root_writes = .true.)
    case (CPU_ALGO)
      call zgemm(transa, transb, &
                 m, n, k, &
                 alpha, a, lda, b, ldb, &
                 beta, c, ldc)
    case default
      call die("Invald algorithm for acc_zgemm", only_root_writes = .true.)
    end select
  end subroutine acc_zgemm
  subroutine acc_dgemm(transa, transb, &
                       m, n, k, &
                       alpha, a, lda, b, ldb, &
                       beta, c, ldc, &
                       algo)
    implicit none
    character :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    double precision :: alpha, beta
    double precision :: a(lda, *), b(ldb, *), c(ldc, *)
    integer(kind(CPU_ALGO)) :: algo
    select case (algo)
    case (OPENACC_ALGO)
      call die("OpenACC version of acc_dgemm requested, but OpenACC not "&
            &"compiled into this executable", only_root_writes = .true.)
    case (CPU_ALGO)
      call dgemm(transa, transb, &
                 m, n, k, &
                 alpha, a, lda, b, ldb, &
                 beta, c, ldc)
    case default
      call die("Invald algorithm for acc_dgemm", only_root_writes = .true.)
    end select
  end subroutine acc_dgemm
end module acc_linalg_m
