module exposed

! to "expose" fields of objects of derived types 
! "Exposed variables" are the global variables defined in this module
! Call expose_* subroutines before accessing them. These subroutines
! simply copy the relevant field to global variables.


use kinds, only: dp
!
implicit none
!
real(dp), allocatable :: rho_of_r(:,:)
integer :: dfftp_nnr

!-----------------------------------------
  contains
!-----------------------------------------
! copy to global variable


!----------------------------
subroutine expose_dfftp_nnr()
!----------------------------
  use fft_base, only: dfftp
  implicit none
  dfftp_nnr = dfftp%nnr
  return
end subroutine




!---------------------------
subroutine expose_rho_of_r()
!---------------------------
  use lsda_mod, only : nspin
  use fft_base, only: dfftp
  use scf, only: scf_type, rho
  implicit none
  !
  if(allocated(rho_of_r)) deallocate(rho_of_r)
  allocate( rho_of_r(dfftp%nnr, nspin) )
  !
  rho_of_r(:,:) = rho%of_r(:,:) ! copy
  !
  return
end subroutine



end module