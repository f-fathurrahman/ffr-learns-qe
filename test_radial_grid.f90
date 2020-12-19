PROGRAM test_radial_grid
  USE kinds, ONLY: DP
  USE radial_grids, ONLY: series
  REAL(DP):: f(4),r(4),r2(4),b(0:3)

  f = (/ 1.d0, 2.d0, 3.d0, 4.d0 /)
  r = (/ 1.1d0, 1.2d0, 1.3d0, 1.4d0 /)
  r2 = r**2
  b = (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
  CALL series(f,r,r2,b)

  WRITE(*,*) 'b = ', b

  WRITE(*,*) 'Pass here ...'
END PROGRAM 

