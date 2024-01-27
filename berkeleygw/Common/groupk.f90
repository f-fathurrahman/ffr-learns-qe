!===============================================================================
!
! Routines:
!
! 1. groupk() Originally By gsm Last Modified 8/31/2010 (gsm)
!
! Generates symmetries of the k-point.
!
!===============================================================================
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
module groupk_m
  use global_m
  implicit none
  private
  public :: groupk
contains
subroutine groupk(kvec,nsyml,rotl,nsymk,rotk)
  real(DP), intent(in) :: kvec(:) !< (3)
  integer, intent(in) :: nsyml
  integer, intent(in) :: rotl(:,:,:) !< (3,3,48)
  integer, intent(out) :: nsymk
  integer, intent(out) :: rotk(:,:,:) !< (3,3,48)
  integer :: isyml
  real(DP) :: krot(3)
 
  nsymk=0
  do isyml=1,nsyml
    krot(1:3) = matmul(dble(rotl(1:3, 1:3, isyml)), kvec(1:3))
    if (all(abs(krot(1:3)-kvec(1:3)).lt.TOL_Small)) then
      nsymk=nsymk+1
      rotk(1:3, 1:3, nsymk) = rotl(1:3, 1:3, isyml)
    endif
  enddo
 
  return
end subroutine groupk
end module groupk_m
