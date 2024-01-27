!==========================================================================
!
! Routines:
!
! (1) find_kpt_match() Originally by SIB
!
! Look for rkq in the list of kpoints rotated by the symmetries.
! If found, it puts its index in ikrkq, itqq is the index
! of the symmetry that worked, and kgqq is an "umklapp" vector (i.e.
! integer 3-vector) so that rkq = symm*kvec + kgqq.
!
!==========================================================================
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
module find_kpt_match_m
  use global_m
  implicit none
  private
  public :: find_kpt_match
contains
  subroutine find_kpt_match(kp_point, syms, rkq, ikrkq, itqq, kgqq)
    type(kpoints), intent(in) :: kp_point
    type(symmetry), intent(in) :: syms
    real(DP), intent(in) :: rkq(3)
    integer, intent(out) :: ikrkq
    integer, intent(out) :: itqq
    integer, intent(out) :: kgqq(3)
    integer :: ik, itq, ii
    real(DP) :: qk(3), del(3)
   
    ikrkq = 0
    ik_loop: do ik = 1, kp_point%nrk
      do itq = 1, syms%ntran
        qk(1:3) = matmul(syms%mtrx(1:3, 1:3, itq), kp_point%rk(1:3, ik))
        do ii = 1, 3
          del(ii) = rkq(ii) - qk(ii)
          if(del(ii).ge.0.0d0) kgqq(ii) = del(ii) + TOL_Small
          if(del(ii).lt.0.0d0) kgqq(ii) = del(ii) - TOL_Small
        enddo
        if(all(abs(del(1:3)-kgqq(1:3)) .lt. TOL_Small)) then
          ikrkq=ik
          itqq = itq
          exit ik_loop
        endif
      enddo
    enddo ik_loop
   
    return
  end subroutine find_kpt_match
end module find_kpt_match_m
