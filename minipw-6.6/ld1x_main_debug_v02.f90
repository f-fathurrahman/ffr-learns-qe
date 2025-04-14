!---------------------------------------------------------------
PROGRAM main_ld1x_debug_v02
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

  CALL mp_startup()
  CALL environment_start(code)

  CALL ld1_readin(input_file_)
  
  CALL ld1x_debug_setup_gen_and_test()
  ! only until this setup

  CALL mp_global_end()

  write(*,*)
  write(*,*) 'Program ended normally'

END PROGRAM
