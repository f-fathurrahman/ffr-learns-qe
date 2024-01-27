!==========================================================================
!
! Routines:
!
! 1. subgrp() Originally By ? Last Modified 6/12/2008 (JRD)
!
! Determines a subgroup of the symmetry group that preserves a q-vector.
!
!===========================================================================
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
module subgrp_m
  use global_m
  use misc_m
  implicit none
  private
  public :: &
    subgrp
contains
subroutine subgrp(qq, syms)
  type (symmetry), intent(inout) :: syms
  real(DP), intent(in) :: qq(3)
  integer :: it, kg(3)
  real(DP) :: qk(3), dqk(3)
 
!---------------------
! Loop over transformations testing for r(q) = q + kg0
  syms%rq = qq
  syms%ntranq = 0
  do it = 1, syms%ntran
    qk(1:3) = matmul(syms%mtrx(1:3, 1:3, it), qq(1:3))
    dqk(1:3) = qk(1:3) - qq(1:3)
    call k_range(dqk(1:3), kg(1:3), TOL_Small)
    if (all(abs(dqk(1:3)) .lt. TOL_Small)) then
      !--------------------
      ! Store index of element of subgroup
      syms%ntranq = syms%ntranq + 1
      syms%indsub(syms%ntranq) = it
      syms%kgzero(1:3, syms%ntranq) = -kg(1:3)
    endif
  enddo
 
  return
end subroutine subgrp
end module subgrp_m
