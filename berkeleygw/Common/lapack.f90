!================================================================================
!
! Modules:
!
! (1) lapack_m Originally By DAS Last Modified 1/13/2011 (das)
!
! Interfaces for LAPACK functions, taken from http://www.netlib.org/lapack/double
! and http://www.netlib.org/lapack/complex16.
! Every LAPACK function used in the code should be listed here, and this
! module should be used in every routine containing LAPACK calls to ensure
! the argument types are correct.
!
! Note that if any array name from netlib.org is X, the interface will
! be interpreted as a preprocessor macro and cause a compilation failure,
! solved by changed to lower-case x.
!
!================================================================================
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
module lapack_m
  public ! only interfaces in this module
  interface
    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      implicit none
      INTEGER INFO, LDA, LDB, N, NRHS
      INTEGER IPIV( * )
      DOUBLE PRECISION A( LDA, * ), B( LDB, * )
    end SUBROUTINE DGESV
  end interface
  interface
    SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      implicit none
      INTEGER INFO, LDA, LDB, N, NRHS
      INTEGER IPIV( * )
      COMPLEX*16 A( LDA, * ), B( LDB, * )
    end SUBROUTINE ZGESV
  end interface
  interface
    SUBROUTINE ZHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, &
      ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO )
      implicit none
      CHARACTER JOBZ, RANGE, UPLO
      INTEGER IL, INFO, IU, LDZ, M, N
      DOUBLE PRECISION ABSTOL, VL, VU
      INTEGER IFAIL( * ), IWORK( * )
      DOUBLE PRECISION RWORK( * ), W( * )
      COMPLEX*16 AP( * ), WORK( * ), Z( LDZ, * )
    end SUBROUTINE ZHPEVX
  end interface
  interface
    SUBROUTINE ZHEEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
      ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL, INFO )
      implicit none
      CHARACTER JOBZ, RANGE, UPLO
      INTEGER IL, INFO, IU, LDA, LDZ, LWORK, M, N
      DOUBLE PRECISION ABSTOL, VL, VU
      INTEGER IFAIL( * ), IWORK( * )
      DOUBLE PRECISION RWORK( * ), W( * )
      COMPLEX*16 A( LDA, * ), WORK( * ), Z( LDZ, * )
    end SUBROUTINE ZHEEVX
  end interface
  interface
    SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
      implicit none
      INTEGER INFO, LDA, M, N
      INTEGER IPIV( * )
      DOUBLE PRECISION A( LDA, * )
    end SUBROUTINE DGETRF
  end interface
  interface
    SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      implicit none
      INTEGER INFO, LDA, LWORK, N
      INTEGER IPIV( * )
      DOUBLE PRECISION A( LDA, * ), WORK( * )
    end SUBROUTINE DGETRI
  end interface
  interface
    SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
      implicit none
      INTEGER INFO, LDA, M, N
      INTEGER IPIV( * )
      COMPLEX*16 A( LDA, * )
    end SUBROUTINE ZGETRF
  end interface
  interface
    SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      implicit none
      INTEGER INFO, LDA, LWORK, N
      INTEGER IPIV( * )
      COMPLEX*16 A( LDA, * ), WORK( * )
    end SUBROUTINE ZGETRI
  end interface
  interface
    SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
      implicit none
      CHARACTER JOBZ, UPLO
      INTEGER INFO, LDA, LWORK, N
      DOUBLE PRECISION RWORK( * ), W( * )
      COMPLEX*16 A( LDA, * ), WORK( * )
    end SUBROUTINE ZHEEV
  end interface
  interface
    SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      implicit none
      CHARACTER JOBZ, UPLO
      INTEGER INFO, LDA, LWORK, N
      DOUBLE PRECISION A( LDA, * ), W( * ), WORK( * )
    end SUBROUTINE DSYEV
  end interface
  interface
    SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      implicit none
      INTEGER INFO, LDA, LWORK, M, N
      DOUBLE PRECISION A( LDA, * ), TAU( * ), WORK( * )
    end SUBROUTINE DGEQRF
  end interface
  interface
    SUBROUTINE ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      implicit none
      INTEGER INFO, LDA, LWORK, M, N
      COMPLEX*16 A( LDA, * ), TAU( * ), WORK( * )
    end SUBROUTINE ZGEQRF
  end interface
  interface
    SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      implicit none
      INTEGER INFO, K, LDA, LWORK, M, N
      DOUBLE PRECISION A( LDA, * ), TAU( * ), WORK( * )
    end SUBROUTINE DORGQR
  end interface
  interface
    SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      implicit none
      INTEGER INFO, K, LDA, LWORK, M, N
      COMPLEX*16 A( LDA, * ), TAU( * ), WORK( * )
    end SUBROUTINE ZUNGQR
  end interface
  interface
    SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
      implicit none
      CHARACTER UPLO
      INTEGER INFO, LDA, N
      DOUBLE PRECISION A( LDA, * )
    end SUBROUTINE DPOTRF
  end interface
end module lapack_m
