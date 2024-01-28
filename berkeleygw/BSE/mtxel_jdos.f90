!===============================================================================
!
! Routines:
!
! (1) mtxel_jdos() Originally By MJ Last Modified 22/7/2010 MJ
!
! input: s1 and nmat
!
! output: s1 matrix element of the JDOS operator
!
! Calculates the JDOS operator
! exp(phi(ic,iv)) where phi(ic,iv) = random number between 0 and 2pi
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
module mtxel_jdos_m
  use global_m
  use random_m
  implicit none
  private
  public :: &
    mtxel_jdos
contains
subroutine mtxel_jdos(s1,nmat)
  integer, intent(in) :: nmat
  real(DP), intent(inout) :: s1(nmat)
  integer :: i
  real(DP) :: rnd,sq2
  complex(DPC) :: zi
  real(DP) :: random
 
  s1 = 0.0d0
  zi = (0.0d0,1.0d0)
  sq2 = sqrt(2.0d0)
!----------------------------------
! Calculate s1(iv) = exp(phi)
  do i = 1, nmat
    call genrand_real4(rnd)
    random = sq2*cos(2.0d0*PI_D*rnd)
    s1(i) = random
  enddo
 
  return
end subroutine mtxel_jdos
end module mtxel_jdos_m
