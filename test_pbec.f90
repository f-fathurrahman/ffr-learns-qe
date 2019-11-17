PROGRAM test_pbec
  IMPLICIT NONE 
  REAL(8) :: rho, grho
  REAL(8) :: sc, v1c, v2c
  INTEGER :: iflag

  rho = 1.2d0
  grho = 0.1d0
  iflag = 1
  CALL pbec( rho, grho, iflag, sc, v1c, v2c )

  WRITE(*,'(1x,F18.10)') sc  

END PROGRAM 
