! ffr, 22 Feb 2016

PROGRAM t_mp
  USE io_global, ONLY : stdout, ionode
  USE mp_global, ONLY : mp_startup, mp_global_end
  USE mp_world, ONLY : nproc, world_comm
  USE mp, ONLY : mp_size
  USE command_line_options, ONLY : nimage_, npool_, ndiag_, nband_, ntg_
  IMPLICIT NONE
  INTEGER :: exit_status
  INTEGER :: ncpu
  
  CALL mp_startup()
 
  IF(ionode) THEN
    WRITE(stdout,*) 'Hello World!'
    WRITE(stdout,*) 'nproc = ', nproc
    WRITE(stdout,*) 'nimage_ = ', nimage_
    WRITE(stdout,*) 'npool_  = ', npool_
    WRITE(stdout,*) 'ndiag_  = ', ndiag_
    WRITE(stdout,*) 'nband_  = ', nband_
    WRITE(stdout,*) 'ntg_    = ', ntg_
  ENDIF

  CALL mp_global_end()  ! only this call that is required from stop_run()
  exit_status = 0
  CALL do_stop( exit_status )

END PROGRAM

