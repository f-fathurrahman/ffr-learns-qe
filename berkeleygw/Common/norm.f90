!===============================================================================
!
! Routines:
!
! 1. norm_wfng() Originally By gsm Last Modified 9/1/2010 (gsm)
!
! Normalizes wavefunctions in G-space.
!
!===============================================================================

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
