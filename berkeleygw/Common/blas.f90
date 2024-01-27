!================================================================================
!
! Modules:
!
! (1) blas_m Originally By DAS Last Modified 1/13/2011 (das)
!
! Interfaces for BLAS functions, taken from http://www.netlib.org/blas/
! Every BLAS function used in the code should be listed here, and this
! module should be used in every routine containing BLAS calls to ensure
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
module blas_m
  public ! only interfaces in this module
  interface
    SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      implicit none
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
    end subroutine DGEMM
  end interface
  interface
    SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      implicit none
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      DOUBLE COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
    end subroutine ZGEMM
  end interface
  interface
    SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,x,INCX,BETA,Y,INCY)
      implicit none
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
      DOUBLE PRECISION A(LDA,*),x(*),Y(*)
    end SUBROUTINE DGEMV
  end interface
  interface
    SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,x,INCX,BETA,Y,INCY)
      implicit none
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
      DOUBLE COMPLEX A(LDA,*),x(*),Y(*)
    end SUBROUTINE ZGEMV
  end interface
  interface
    SUBROUTINE ZHEMV(UPLO,N,ALPHA,A,LDA,x,INCX,BETA,Y,INCY)
      implicit none
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
      DOUBLE COMPLEX A(LDA,*),x(*),Y(*)
    end SUBROUTINE ZHEMV
  end interface
  interface
    SUBROUTINE ZHPMV(UPLO,N,ALPHA,AP,x,INCX,BETA,Y,INCY)
      implicit none
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER INCX,INCY,N
      CHARACTER UPLO
      DOUBLE COMPLEX AP(*),x(*),Y(*)
    end SUBROUTINE ZHPMV
  end interface
  interface
    SUBROUTINE ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      implicit none
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER LDA,LDB,LDC,M,N
      CHARACTER SIDE,UPLO
      DOUBLE COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
    end SUBROUTINE ZHEMM
  end interface
  interface blas_dot
    DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
      implicit none
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DX(*),DY(*)
    end FUNCTION DDOT
    DOUBLE COMPLEX FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)
      implicit none
      INTEGER INCX,INCY,N
      DOUBLE COMPLEX ZX(*),ZY(*)
    end FUNCTION ZDOTC
  end interface blas_dot
  interface
    SUBROUTINE DSCAL(N,DA,DX,INCX)
      implicit none
      DOUBLE PRECISION DA
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
    end SUBROUTINE DSCAL
  end interface
  interface
    SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
      implicit none
      DOUBLE COMPLEX ZA
      INTEGER INCX,N
      DOUBLE COMPLEX ZX(*)
    end SUBROUTINE ZSCAL
  end interface
  interface blas_nrm2
    DOUBLE PRECISION FUNCTION DNRM2(N,x,INCX)
      implicit none
      INTEGER INCX,N
      DOUBLE PRECISION x(*)
    end FUNCTION DNRM2
    DOUBLE PRECISION FUNCTION DZNRM2(N,x,INCX)
      implicit none
      INTEGER INCX,N
      DOUBLE COMPLEX x(*)
    end FUNCTION DZNRM2
  end interface
end module blas_m
