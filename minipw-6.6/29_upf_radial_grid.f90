include 'prepare_all.f90'

PROGRAM main
  IMPLICIT NONE 

  CALL prepare_all()
  call test_radial_grid()

END PROGRAM 



SUBROUTINE test_radial_grid()
  use uspp_param, only: upf
  implicit none

  write(*,*) 'upf%dx    = ', upf(1)%dx
  write(*,*) 'upf%xmin  = ', upf(1)%xmin
  write(*,*) 'upf%zmesh = ', upf(1)%zmesh
  write(*,*) 'upf%mesh  = ', upf(1)%mesh

END SUBROUTINE

