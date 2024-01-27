!===============================================================================
!
! Routines:
!
! 1. norm_wfng() Originally By gsm Last Modified 9/1/2010 (gsm)
!
! Normalizes wavefunctions in G-space.
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
module norm_m
  use blas_m
  use global_m
  implicit none
  private
  public :: &
    norm_wfng
contains
subroutine norm_wfng(ngk_l,nbstart,nbend,nb,ns,nk,wfn_d)
  integer, intent(in) :: ngk_l,nbstart,nbend,nb,ns,nk
  real(DP), intent(inout) :: wfn_d(ngk_l,nb,ns,nk) !< parallelized over G-vectors
  integer :: ib,is,ik,nbnorm
  character(len=16) :: s1,s2,s3
  real(DP), allocatable :: norm2band(:)
  real(DP), allocatable :: norm2dummy(:)
 
  nbnorm=nbend-nbstart+1
  allocate(norm2band (nbnorm))
  do ik=1,nk
    do is=1,ns
      do ib=nbstart,nbend
        norm2band(ib-nbstart+1) = blas_nrm2(ngk_l, wfn_d(:,ib,is,ik), 1)**2
      enddo ! ib
      do ib=nbstart,nbend
        if (norm2band(ib-nbstart+1) .gt. TOL_Zero) then
          call dscal(ngk_l, 1.0d0/sqrt(norm2band(ib-nbstart+1)), wfn_d(:,ib,is,ik), 1)
        else ! norm2.gt.TOL_Zero
          if (peinf%inode.eq.0) then
            write(s1,111)ib
            write(s2,111)is
            write(s3,111)ik
            write(0,211)TRUNC(s1),TRUNC(s2),TRUNC(s3)
          endif ! peinf%inode.eq.0
        endif ! norm2.gt.TOL_Zero
      enddo ! ib
    enddo ! is
  enddo ! ik
  if(allocated(norm2band))then;deallocate(norm2band);endif
 
  return
111 format(i16)
211 format(1x,"WARNING: zero norm for k-point =",1x,a,1x,"spin =", &
      1x,a,1x,"band =",1x,a,/)
end subroutine norm_wfng
end module norm_m
