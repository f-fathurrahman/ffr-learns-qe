program test
  implicit none
  real(8) :: besr

  call sph_bes( 1, (/ 1.1d0 /), 2.d0, 0, besr )
  write(*,*) 'besr = ', besr

  call sph_bes( 1, (/ 1.1d0 /), 2.d0, 1, besr )
  write(*,*) 'besr = ', besr

  call sph_bes( 1, (/ 1.1d0 /), 2.d0, 2, besr )
  write(*,*) 'besr = ', besr

  call sph_bes( 1, (/ 1.1d0 /), 2.d0, 3, besr )
  write(*,*) 'besr = ', besr


end program