SUBROUTINE prepare_all()

  USE mp_global, ONLY : mp_startup
  USE mp_world, ONLY : world_comm
  USE mp_pools, ONLY : intra_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE command_line_options, ONLY : ndiag_, input_file_
  USE environment, ONLY : environment_start
  USE read_input, ONLY : read_input_file
  USE check_stop, ONLY : check_stop_init

  IMPLICIT NONE 

  include 'laxlib.fh'
  
  CALL mp_startup(start_images=.true.)

  CALL laxlib_start( ndiag_, world_comm, intra_pool_comm, &
                     do_distr_diag_inside_bgrp_ = .FALSE. )

  CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, inter_bgrp_comm )

  CALL environment_start( 'PWSCF' )

  CALL read_input_file( 'PW', input_file_ )

  CALL iosys()

  CALL check_stop_init()  ! required in c_bands

  CALL setup()
  CALL my_init_run()

  !CALL unset_mpi_comm_4_solvers()
  !CALL stop_run( exit_status )
  !CALL do_stop( exit_status )

  WRITE(*,*) 'prepare_all ended normally'

END SUBROUTINE 
