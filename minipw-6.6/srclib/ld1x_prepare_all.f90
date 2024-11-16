subroutine ld1x_prepare_all()
  USE mp_global, ONLY : mp_startup
  USE environment, ONLY : environment_start
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  CHARACTER(LEN=9) :: code = 'LD1'
  
  write(*,*)
  write(*,*) '************ ENTER ld1x_prepare_all ***********'
  write(*,*)

  !
  ! write initialization information
  !
  CALL mp_startup()
  CALL environment_start(code)
  !
  ! read input, possible pseudopotential and set the main variables
  !
  CALL ld1_readin(input_file_)
  CALL ld1_setup()

  write(*,*)
  write(*,*) '************ EXIT ld1x_prepare_all ***********'
  write(*,*)

end subroutine

