!============================================================================
!
! Routines:
!
! (1) check_screening_trunc Originally by JRD Last Modified: 2/09/2009 (JRD)
!
! Die if screening, truncation, and q0vec are not set consistently.
!
!============================================================================
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
module check_screening_m
  use global_m
  implicit none
  private
  public :: check_screening_trunc
contains
subroutine check_screening_trunc(itruncflag,iscreen,q0vec,bdot)
  integer, intent(in) :: itruncflag
  integer, intent(in) :: iscreen
  real(DP), intent(in) :: q0vec(3)
  real(DP), intent(in) :: bdot(3,3)
  real(DP) :: q0len
 
  q0len = sqrt(DOT_PRODUCT(q0vec,MATMUL(bdot,q0vec)))
  if (iscreen==SCREEN_METAL .and. q0len<TOL_SMALL) then
    if(peinf%inode == 0) then
      write(0,*) ' '
      write(0,*) 'You want metallic screening but didn''t specify q0vec!!'
    endif
    call die('Inconsistent Screening', only_root_writes = .true.)
  endif
  if (iscreen==SCREEN_GRAPHENE .and. q0len<TOL_SMALL .and. itruncflag==TRUNC_NONE) then
    if(peinf%inode == 0) then
      write(0,*) ' '
      write(0,*) 'You want graphene screening with no truncation'
      write(0,*) 'but didn''t specify q0vec!!'
    endif
    call die('Inconsistent Screening', only_root_writes = .true.)
  endif
  if ((itruncflag==TRUNC_NONE .or. itruncflag==TRUNC_WIRE .or. &
    itruncflag==TRUNC_SLAB) .and. q0len<TOL_SMALL) then
    if(peinf%inode == 0) then
      write(0,*) ' '
      write(0,*) 'You have a divergent Coulomb interaction but didn''t specify q0vec!!'
    endif
    call die('Inconsistent Screening', only_root_writes = .true.)
  endif
  if ((itruncflag/=TRUNC_NONE .and. itruncflag/=TRUNC_WIRE .and. &
    itruncflag/=TRUNC_SLAB) .and. iscreen/=SCREEN_METAL .and. q0len>=TOL_SMALL) then
    if(peinf%inode == 0) then
      write(0,*) ''
      write(0,*) 'You want semiconductor or graphene screening with truncation'
      write(0,*) 'but specified nonzero q0vec!!'
    endif
    call die('Inconsistent Screening', only_root_writes = .true.)
  endif
 
  return
end subroutine check_screening_trunc
end module check_screening_m
