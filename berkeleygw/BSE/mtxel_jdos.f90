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
