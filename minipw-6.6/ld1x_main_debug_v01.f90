!---------------------------------------------------------------
PROGRAM main_ld1x_debug_v01
!---------------------------------------------------------------
  !
  USE mp_global,         ONLY : mp_startup, mp_global_end
  USE environment,       ONLY : environment_start
  USE ld1inc,            ONLY : iswitch, grid
  USE radial_grids,      ONLY : deallocate_radial_grid
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  CHARACTER(LEN=9) :: code = 'LD1'
  !
  !   write initialization information
  !
  CALL mp_startup()
  CALL environment_start(code)
  !
  !    read input, possible pseudopotential and set the main variables
  !
  CALL ld1_readin(input_file_)
  CALL ld1_setup()


  call ld1x_print_variables()


  ! Errors on these cases
  IF( iswitch == 1 ) THEN
    
    WRITE(*,*) 'DEBUG STARTS HERE'
    CALL ld1x_debug_v01()

  ELSEIF( iswitch == 2) THEN
    call errore('debug_ld1x_v01', 'iswitch = 2 is not supported', 1)
  ELSEIF( iswitch == 3 ) THEN
     !
     !  pseudopotential generation and test
     !
     CALL all_electron(.FALSE., 1) ! do not compute log-deriv
     CALL my_gener_pseudo()
     !if(.not. lgipaw_reconstruction) 
     CALL run_test()
     CALL ld1_writeout()
  ELSEIF( iswitch == 4 ) THEN
    call errore('debug_ld1x_v01', 'iswitch = 4 is not supported', 1)
  ELSE 
    CALL errore('ld1', 'iswitch not implemented',1)
  ENDIF 


  ! XXX: why call this here?
  CALL deallocate_radial_grid( grid )

  CALL mp_global_end()

  write(*,*)
  write(*,*) 'Program ended normally'

END PROGRAM
