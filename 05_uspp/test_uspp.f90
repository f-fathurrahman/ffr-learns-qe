PROGRAM test
  
  USE mp_global, ONLY : mp_startup, mp_global_end
  USE environment, ONLY : environment_start
  USE read_input, ONLY : read_input_file
  USE command_line_options, ONLY : input_file_

  IMPLICIT NONE 
  !INTEGER :: exit_status

  CALL mp_startup( )
  CALL environment_start( 'PWSCF' )
  WRITE(*,*) 'input_file_ = ', trim(input_file_)
  CALL read_input_file( 'PW', input_file_ )
  CALL iosys()
  CALL setup()
  CALL init_run()

  CALL test_uspp()

  CALL mp_global_end()

END PROGRAM 



!---------------------
SUBROUTINE test_uspp()

  USE uspp, ONLY : nkb, nhtol, nhtolm, indv

  IMPLICIT NONE 

  WRITE(*,*) 'nkb = ', nkb

END SUBROUTINE 

