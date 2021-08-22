PROGRAM simple_bse

  USE start_end
  USE input_simple_exc
  
  IMPLICIT NONE 

  TYPE(input_options) :: sinp

  !setup MPI environment
  CALL startup()

  CALL read_input_simple_exc( sinp )  
  SELECT CASE(sinp%task)
  CASE(0) !solve eigen-problem
    CALL simple_eigen(sinp)
  CASE(1) !find spectrum lanczos
    CALL lanczos(sinp)
  END SELECT

  CALL stop_run()
  STOP

END PROGRAM simple_bse
