!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!----------------------------------------------------------------
subroutine fft_interpolate_real( dfft_in, v_in, dfft_out, v_out )
!----------------------------------------------------------------
  !
  !   This subroutine interpolates an array  v_in   defined on fft grid  dfft_in
  !                           to   an array  v_out  defined on fft grid  dfft_out
  !   v_in and v_out are assumed to be real arrays and may concide
  !
  USE fft_param,      ONLY : DP
  USE fft_types,      ONLY : fft_type_descriptor
  USE fft_interfaces, ONLY : fwfft, invfft
  ! I/O variables
  TYPE(fft_type_descriptor), INTENT(IN) :: dfft_in, dfft_out
  REAL(DP), INTENT(IN) :: v_in(:) !dfft_in%nnr)
  REAL(DP), INTENT(OUT) :: v_out(:) !dfft_out%nnr)
  ! local variables
  integer :: i
  INTEGER :: ngm
  COMPLEX(DP), ALLOCATABLE :: aux_in (:), aux_out (:)

  call start_clock ('interpolate')

  write(*,*)
  write(*,*) 'fft_interpolate_real: shape v_in = ', shape(v_in)
  write(*,*) 'fft_interpolate_real: shape v_out = ', shape(v_out)
  write(*,*)
  write(*,*) 'dfft_in%nnr  = ', dfft_in%nnr
  write(*,*) 'dfft_out%nnr = ', dfft_out%nnr
  write(*,*)


  IF( dfft_out%grid_id == dfft_in%grid_id ) THEN
     
    write(*,*) 'fft_interpolate_real: simply copy' 
    v_out (1:dfft_in%nnr) = v_in (1:dfft_in%nnr)

  ELSE
     
    write(*,*) 'fft_interpolate_real: two different grids'
    write(*,*) 'fft_interpolate_real: dfft_in%nnr = ', dfft_in%nnr
    write(*,*) 'fft_interpolate_real: dfft_out%nnr = ', dfft_out%nnr

    if (dfft_in%lgamma .neqv. dfft_out%lgamma) &
       call fftx_error__ ('fft_interpolate_real','two grids with inconsistent lgamma values', 1)

    ALLOCATE( aux_in(dfft_in%nnr) )
    ALLOCATE( aux_out(dfft_out%nnr) )

    ! Copy all data to aux_in
    aux_in(1:dfft_in%nnr) = v_in(1:dfft_in%nnr)

    ! To G-space
    CALL fwfft('Rho', aux_in, dfft_in)

    write(*,*) 'Some aux_in after fwfft'
    do i = 1,10
      write(*,'(1x,I5,2F18.10)') i, aux_in(i)
    enddo


    ! Zero out aux_out
    aux_out(1:dfft_out%nnr) = (0.d0, 0.d0)

    ! find the minimum of ngm for both FFT grids
    ngm = min(dfft_in%ngm, dfft_out%ngm)
    write(*,*) 'fft_interpolate_real: ngm = ', ngm

    ! The indices dfft_in%nl(1:ngm)
    ! Copy aux_in to aux_out
    aux_out( dfft_out%nl(1:ngm) ) = aux_in( dfft_in%nl(1:ngm) )
    
    write(*,*) 'Some aux_out before invfft'
    do i = 1,10
      write(*,'(1x,I5,2F18.10)') i, aux_out(i)
    enddo


    ! additional work for gamma-trick
    IF(dfft_in%lgamma) aux_out( dfft_out%nlm(1:ngm) ) = aux_in( dfft_in%nlm(1:ngm) )

    write(*,*) 'Before invfft: sum(aux_out) = ', sum(aux_out)
    
    ! Back to R-space
    CALL invfft('Rho', aux_out, dfft_out)

    write(*,*) 'Some aux_out after invfft'
    do i = 1,10
      write(*,'(1x,I5,2F18.10)') i, aux_out(i)
    enddo
    
    !
    write(*,*) 'After invfft: sum(aux_out) = ', sum(aux_out)

    v_out(1:dfft_out%nnr) = aux_out(1:dfft_out%nnr)

    DEALLOCATE(aux_in, aux_out)

  ENDIF

  call stop_clock ('interpolate')

  return

end subroutine fft_interpolate_real


!-------------------------------------------------------------------
subroutine fft_interpolate_complex (dfft_in, v_in, dfft_out, v_out )
!-------------------------------------------------------------------
  !
  !   This subroutine interpolates an array  v_in   defined on fft grid  dfft_in
  !                           to   an array  v_out  defined on fft grid  dfft_out
  !   v_in and v_out are assumed to be complex arrays and may concide
  !
  USE fft_param,      ONLY : DP
  USE fft_types,      ONLY : fft_type_descriptor
  USE fft_interfaces, ONLY : fwfft, invfft
  ! I/O variables
  TYPE(fft_type_descriptor), INTENT(IN) :: dfft_in, dfft_out
  COMPLEX(DP),INTENT(IN)  :: v_in (:) !dfft_in%nnr)
  COMPLEX(DP),INTENT(OUT) :: v_out (:) !dfft_out%nnr)
  ! local variables
  INTEGER :: ngm
  COMPLEX(DP), ALLOCATABLE :: aux_in (:)

  if (dfft_out%lgamma.OR.dfft_in%lgamma) call fftx_error__('fft_interpolate_complex','lgamma not allowed', 1)

  call start_clock ('interpolate')

  IF (dfft_out%grid_id == dfft_in%grid_id) THEN

     v_out (1:dfft_in%nnr) = v_in (1:dfft_in%nnr)

  ELSE

     ALLOCATE (aux_in( dfft_in%nnr))

     aux_in (1:dfft_in%nnr) = v_in(1:dfft_in%nnr)

     CALL fwfft ('Rho', aux_in, dfft_in)

     v_out(1:dfft_out%nnr) = (0.d0, 0.d0)

     ngm = min(dfft_in%ngm, dfft_out%ngm)

     v_out (dfft_out%nl (1:ngm) ) = aux_in (dfft_in%nl (1:ngm) )

     CALL invfft ('Rho', v_out, dfft_out)

     DEALLOCATE (aux_in)

  END IF

  call stop_clock ('interpolate')

  return

end subroutine fft_interpolate_complex
!
