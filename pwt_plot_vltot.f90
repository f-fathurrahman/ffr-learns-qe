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
    CALL xsf_fast_datagrid_3d(vltot, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
      dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, at, alat, iuxsf)
    CLOSE(iuxsf)
  ENDIF



END SUBROUTINE


!------------------------------------------------------------------------------
PROGRAM pwt_plot_vltot
!------------------------------------------------------------------------------
  USE environment, ONLY: environment_start
  USE mp_global, ONLY: mp_startup
  USE read_input, ONLY: read_input_file
  USE command_line_options, ONLY: input_file_
  !
  USE check_stop, ONLY: check_stop_init
  IMPLICIT NONE

  CALL mp_startup()
  CALL environment_start( 'PWSCF' )
  CALL read_input_file( 'PW', input_file_ )
  
  ! Details of run_pwscf
  CALL iosys()
  CALL check_stop_init()
  CALL setup()
  CALL init_run()

  CALL plot_vltot()

  CALL stop_run( 0 )
  CALL do_stop( 0 )
  STOP

END PROGRAM   




