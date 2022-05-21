SUBROUTINE my_electrons()
!----------------------------------------------------------------------------
!! General self-consistency loop, also for hybrid functionals.  
!! For non-hybrid functionals it just calls \(\texttt{electron_scf}\).

  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE funct,                ONLY : dft_is_hybrid
  USE control_flags,        ONLY : restart
  !
  IMPLICIT NONE
  integer :: printout
  REAL(DP) :: exxen

  IF( restart ) THEN
   ! unsupported
   stop 'restart is not supported'
  ENDIF

  IF( dft_is_hybrid() ) THEN
    stop 'Hybrid DFT is not supported'
  ENDIF

  exxen = 0.0d0
  printout = 2
  CALL my_electrons_scf( printout, exxen )
  FLUSH( stdout )

  RETURN
END SUBROUTINE









