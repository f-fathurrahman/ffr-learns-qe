SUBROUTINE plot_vltot()
  USE io_global, ONLY : stdout, ionode
  USE cell_base, ONLY : alat, at
  USE ions_base, ONLY : nat, tau, atm, ityp
  USE scf, ONLY : vltot
  USE fft_base, ONLY : dfftp
  IMPLICIT NONE
  INTEGER :: iuxsf

  iuxsf = 111
  IF(ionode) THEN
    OPEN(unit=iuxsf, file='STRUCT_vltot.xsf', action='write')
    CALL xsf_struct(alat, at, nat, tau, atm, ityp, iuxsf)
    ! convert vltot to Hartree
    CALL xsf_fast_datagrid_3d(0.5d0*vltot, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
      dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, at, alat, iuxsf)
    CLOSE(iuxsf)
  ENDIF



END SUBROUTINE


!------------------------------------------------------------------------------
PROGRAM pwt_plot_vltot
!------------------------------------------------------------------------------
  
  IMPLICIT NONE 
  CALL prepare_all()

  CALL plot_vltot()

END PROGRAM   




