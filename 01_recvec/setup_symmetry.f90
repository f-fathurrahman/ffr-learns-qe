!------------------------------------------------------------------------------
SUBROUTINE setup_symmetry()
!------------------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE io_global, ONLY : stdout, ionode
  USE symm_base, ONLY : set_sym_bl, find_sym
  USE ions_base, ONLY : nat, tau, ityp
  USE fft_base, ONLY : dfftp
  USE symm_base, ONLY : nrot, nsym
  IMPLICIT NONE
  ! Local
  LOGICAL :: magnetic_sym
  REAL(DP) :: m_loc(3,nat)

  ! We ignore the magnetization
  magnetic_sym = .FALSE.
  m_loc(:,:)   = 0_DP

  CALL set_sym_bl()
  CALL find_sym( nat, tau, ityp, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
    magnetic_sym, m_loc )

  IF(ionode) THEN
    WRITE(stdout,*) 'nat, nsym, nrot = ', nat, nsym, nrot
  ENDIF

END SUBROUTINE
