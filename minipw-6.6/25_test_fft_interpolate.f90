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
  integer :: i

  allocate( vin(dfftp%nnr) )
  allocate( vout(dfftp%nnr) )

  vin(:) = 1.1d0
  vin(1:5) = 2.5d0
  vout(:) = 0.d0

  write(*,*) 'sum vin before fft_interpolate = ', sum(vin)
  write(*,*) 'sum vout before fft_interpolate = ', sum(vout)

  call fft_interpolate( dfftp, vin, dffts, vout )

  write(*,*) 'sum vin after fft_interpolate = ', sum(vin)
  write(*,*) 'sum vout after fft_interpolate = ', sum(vout(1:dffts%nnr))
  write(*,*) 'sum vout after fft_interpolate dfftp = ', sum(vout(1:dfftp%nnr))

  write(*,*) 'Some vin and vout'
  do i = 1,10
    write(*,'(1x,I5,2F18.10)') i, vin(i), vout(i)
  enddo

  deallocate( vin )
  deallocate( vout )

end subroutine