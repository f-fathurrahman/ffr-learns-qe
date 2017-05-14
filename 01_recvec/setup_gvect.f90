!------------------------------------------------------------------------------
SUBROUTINE setup_gvect( gamma_only )
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE recvec_subs, ONLY : ggen
  USE cell_base, ONLY : at, bg
  IMPLICIT NONE
  !
  LOGICAL :: gamma_only

  CALL data_structure( gamma_only )
  CALL ggen( gamma_only, at, bg )
  CALL gshells( .false. )
END SUBROUTINE

