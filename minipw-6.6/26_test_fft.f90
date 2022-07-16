program main

  implicit none

  call prepare_all()

  call test_fft()

end program


!----------------------------
subroutine test_fft()
!----------------------------
  USE fft_base, ONLY : dffts, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  implicit none
  complex(8), allocatable :: vin(:), vout(:)
  integer :: i

  allocate( vin(dfftp%nnr) )
  allocate( vout(dfftp%nnr) )

  vin(:) = 1.d0
  vin(1:10) = 2.5d0

  write(*,*) 'sum vin = ', sum(vin)
  ! To G-space
  vout(:) = vin(:)
  !CALL fwfft('Rho', vout, dfftp)
  CALL invfft('Rho', vout, dfftp)
  write(*,*) 'sum vout = ', sum(vout)

  write(*,*) 'Some vout'
  do i = 1,10
    write(*,'(1x,I5,2F18.10)') i, vout(i)
  enddo

  deallocate( vin )
  deallocate( vout )

end subroutine