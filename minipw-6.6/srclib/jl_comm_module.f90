!-------------
module jl_comm
!-------------

use kinds, only: dp
implicit none

integer :: Nrow_becp_k ! nkb
integer :: Ncol_becp_k ! nbnd
complex(dp), allocatable :: becp_k(:,:)


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



!CALL calbec( n, vkb, psi, becp, m )


end module