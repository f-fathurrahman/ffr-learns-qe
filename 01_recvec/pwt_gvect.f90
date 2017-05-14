PROGRAM pwt_gvect
  USE mp_global, ONLY: mp_startup
  USE environment, ONLY: environment_start
  USE read_input, ONLY: read_input_file
  USE check_stop, ONLY: check_stop_init
  USE parameters, ONLY: ntypx, npk, lmaxx
  ! for v-5.1 and later
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  ! Local
  INTEGER :: exit_status

  CALL mp_startup( )

  CALL environment_start('PWSCF')

  CALL read_input_file('PW', input_file_)

  CALL iosys()

  ! If we don't include the following two subroutine calls, there will
  ! be error before the program stops regarding find_unit.
  CALL setup()
  CALL init_run()

  CALL info_gvect()
  CALL info_gvecs()

  !CALL t_import_gvect()

  CALL info_structure()

  CALL stop_run( exit_status )
  CALL do_stop( exit_status )

END PROGRAM

