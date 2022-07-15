program main

  implicit none

  call prepare_all()

  call test_interpolate()

end program


!----------------------------
subroutine test_interpolate()
!----------------------------
  USE fft_base, ONLY : dffts, dfftp
  USE fft_interfaces,  ONLY : fft_interpolate
  implicit none
  real(8), allocatable :: vin(:), vout(:)

  allocate( vin(dfftp%nnr) )
  allocate( vout(dffts%nnr) )

  vin(:) = 1.d0
  vin(1:10) = 2.5d0
  vout(:) = 0.d0

  write(*,*) 'sum vin before fft_interpolate = ', sum(vin)
  write(*,*) 'sum vout before fft_interpolate = ', sum(vout)

  call fft_interpolate( dfftp, vin, dffts, vout )

  write(*,*) 'sum vin after fft_interpolate = ', sum(vin)
  write(*,*) 'sum vout after fft_interpolate = ', sum(vout)

  deallocate( vin )
  deallocate( vout )

end subroutine