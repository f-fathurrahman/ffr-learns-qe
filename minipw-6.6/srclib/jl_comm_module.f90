!-------------
module jl_comm
!-------------

use kinds, only: dp
implicit none

integer :: Nrow_becp_k ! nkb
integer :: Ncol_becp_k ! nbnd
complex(dp), allocatable :: becp_k(:,:)
real(dp), allocatable :: rho_of_r(:,:)


contains

!-----------------------
subroutine copy_becp_k()
!-----------------------
  use becmod, only: becp
  implicit none
  integer :: m, n

  m = size(becp%k, 1)
  n = size(becp%k, 2)

  Nrow_becp_k = m
  Ncol_becp_k = n

  if(allocated(becp_k)) deallocate(becp_k)
  allocate(becp_k(m,n))

  becp_k(:,:) = becp%k(:,:)

  return

end subroutine


! TODO: make a more general subroutine
! e.g.: export_scf_rho
!--------------------------------------
subroutine copy_rho_of_r( rho_of_r_in )
!--------------------------------------
  implicit none
  ! Argument
  real(8), intent(in) :: rho_of_r_in(:,:)
  ! Local
  integer :: Nrows, Ncols

  Nrows = size(rho_of_r, 1)
  Ncols = size(rho_of_r, 2)

  if( allocated(rho_of_r) ) deallocate(rho_of_r)
  allocate(rho_of_r(Nrows, Ncols))

  rho_of_r(:,:) = rho_of_r_in(:,:)

  write(*,*) 'copied rho_of_r to jl_comm_module'

  return
end



end module