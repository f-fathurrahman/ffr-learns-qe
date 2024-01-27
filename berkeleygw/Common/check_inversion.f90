!===============================================================================
!
! Module:
!
! (1) check_inversion_m Originally By DAS Last Modified 10/14/2010
!
! Check whether our choice of real/complex version is appropriate given the
! presence or absence of inversion symmetry about the origin, and a guess
! about time-reversal symmetry depending on the number of spin-components.
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
module check_inversion_m
  use global_m
  implicit none
  private
  public :: check_inversion, check_inversion_type
contains
subroutine check_inversion(iflavor, ntran, mtrx, nspin, warn, real_need_inv, tnp)
  integer, intent(in) :: iflavor
  integer, intent(in) :: ntran
  integer, intent(in) :: mtrx(3, 3, 48) !< symmetry operations matrices
  integer, intent(in) :: nspin
  logical, intent(in) :: warn !< set to false to suppress warnings, for converters
  logical, intent(in) :: real_need_inv !< use for generating routines to block real without inversion
     !! this is not always true so that it is possible to run real without using symmetries
  real(DP), optional, intent(in) :: tnp(3, 48) !< fractional translations.
     !! optional only to avoid changing external interface for library.
  integer :: invflag, isym, ii, jj, itest, real_or_complex
  logical :: origin_inv
  character(len=7) :: sflavor
 
  if(iflavor .eq. 0) then
    real_or_complex = 1
  elseif(iflavor .eq. 1 .or. iflavor .eq. 2) then
    real_or_complex = iflavor
  else
    write(sflavor, '(i7)') iflavor
    call die("Illegal value iflavor = " // TRUNC(sflavor) // " passed to check_inversion: must be 0,1,2.", &
      only_root_writes=.true.)
  endif
  invflag = 0
  origin_inv = .false.
  do isym = 1, ntran
    itest = 0
    do ii = 1, 3
      do jj = 1, 3
        if(ii .eq. jj) then
          itest = itest + (mtrx(ii, jj, isym) + 1)**2
        else
          itest = itest + mtrx(ii, jj, isym)**2
        endif
      enddo
    enddo
    if(itest .eq. 0) then
      invflag = invflag + 1
      if(present(tnp)) then
        if(sum(abs(tnp(1:3, isym))) < TOL_Small) origin_inv = .true.
      else
        origin_inv = .true.
      endif
    endif
  enddo
  if(invflag > 0 .and. .not. origin_inv .and. peinf%inode==0) then
    write(0, '(a)') "WARNING: Inversion symmetry is present only with a fractional translation."
    write(0, '(a)') "Apply the translation so inversion is about the origin, to be able to use the real version."
  endif
  if(invflag .gt. 1 .and. peinf%inode==0) &
    write(0, '(a)') "WARNING: More than one inversion symmetry operation is present."
  if(invflag > 0 .and. .not. present(tnp) .and. peinf%inode==0) then
    write(0, '(a)') "WARNING: check_inversion did not receive fractional translations."
    write(0, '(a)') "Cannot confirm that inversion symmetry is about the origin for use of real version."
  endif
  if(real_or_complex .eq. 2) then
    if(origin_inv .and. warn .and. nspin == 1) then
      if(peinf%inode .eq. 0) &
        write(0, '(a)') "WARNING: Inversion symmetry about the origin is present. The real version would be faster."
    endif
  else
    if(.not. origin_inv) then
      if(real_need_inv) then
        call die("The real version cannot be used without inversion symmetry about the origin.", only_root_writes = .true.)
      endif
      if(peinf%inode .eq. 0) then
        write(0, '(a)') "WARNING: Inversion symmetry about the origin is absent in symmetries used to reduce k-grid."
        write(0, '(a)') "Be sure inversion about the origin is still a spatial symmetry, or you must use complex version instead."
      endif
    endif
    if(nspin > 1) then
      call die("Real version may only be used for spin-unpolarized calculations.", only_root_writes = .true.)
    endif
  endif
 
  return
end subroutine check_inversion
!=========================================================================
!> wrapper routine that uses typedefs types
subroutine check_inversion_type(iflavor, syms, nspin, warn, real_need_inv)
  integer, intent(in) :: iflavor
  type (symmetry), intent(in) :: syms
  integer, intent(in) :: nspin
  logical, intent(in) :: warn
  logical, intent(in) :: real_need_inv !< use for generating routines to block real without inversion
 
  call check_inversion(iflavor, syms%ntran, syms%mtrx, nspin, warn, real_need_inv, tnp = syms%tnp)
 
  return
end subroutine check_inversion_type
end module check_inversion_m
