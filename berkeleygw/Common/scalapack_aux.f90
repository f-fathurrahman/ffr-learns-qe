!================================================================================
!
! Modules:
!
! (1) scalapack_aux_m Originally by FHJ 07/24/2015
!
! Defines functions associated to matrix distributions. These functions
! were copied from netlib ScaLAPACK (BSD-licensed), and are included here
! as they can be used on builds with and without ScaLAPACK. None of these
! functions require ScaLAPACK.
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
module scalapack_aux_m
  implicit none
  private
  public :: &
    numroc, &
    indxl2g, &
    indxg2l, &
    indxg2p
contains
! -- ScaLAPACK tools routine (version 1.7) --
  INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
    integer, intent(in) :: N, NB, IPROC, ISRCPROC, NPROCS
    integer :: EXTRABLKS, MYDIST, NBLOCKS
    MYDIST = MOD( NPROCS+IPROC-ISRCPROC, NPROCS )
    NBLOCKS = N / NB
    NUMROC = (NBLOCKS/NPROCS) * NB
    EXTRABLKS = MOD( NBLOCKS, NPROCS )
    IF( MYDIST.LT.EXTRABLKS ) THEN
      NUMROC = NUMROC + NB
    ELSE IF( MYDIST.EQ.EXTRABLKS ) THEN
      NUMROC = NUMROC + MOD( N, NB )
    END IF
  END FUNCTION NUMROC
  INTEGER FUNCTION INDXL2G( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )
    integer, intent(in) :: INDXLOC, IPROC, ISRCPROC, NB, NPROCS
    INDXL2G = NPROCS*NB*((INDXLOC-1)/NB) + MOD(INDXLOC-1,NB) + &
      MOD(NPROCS+IPROC-ISRCPROC, NPROCS)*NB + 1
  END FUNCTION INDXL2G
  INTEGER FUNCTION INDXG2L( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
    integer, intent(in) :: INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
    INDXG2L = NB*((INDXGLOB-1)/(NB*NPROCS))+MOD(INDXGLOB-1,NB)+1
  END FUNCTION INDXG2L
  INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
    integer, intent(in) :: INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
    INDXG2P = MOD( ISRCPROC + (INDXGLOB - 1) / NB, NPROCS )
  END FUNCTION INDXG2P
end module scalapack_aux_m
