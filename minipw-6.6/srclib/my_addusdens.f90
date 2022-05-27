!----------------------------------------------------------------------
SUBROUTINE my_addusdens( rho )
  !----------------------------------------------------------------------
  !! Add US contribution to the charge density to \(\text{rho}(G)\).
  !
  USE control_flags,        ONLY : tqr
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_base,             ONLY : dfftp
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: rho(dfftp%ngm,nspin_mag)
  !! Charge density in G space
  !
  IF ( tqr ) THEN
     stop 'stop tqr == .true. is not supported here'
  ELSE
     CALL my_addusdens_g( rho )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE
!
