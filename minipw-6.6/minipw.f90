program minipw
  use mp_global, only : mp_startup
  use mp_world, only : world_comm
  use mp_pools, only : intra_pool_comm
  use mp_bands, only : intra_bgrp_comm, inter_bgrp_comm
  use command_line_options, only : ndiag_, input_file_
  use environment, only : environment_start
  use read_input, only : read_input_file
  use check_stop, only : check_stop_init

  IMPLICIT NONE 

  include 'laxlib.fh'
  
  call mp_startup(start_images=.true.)

  call laxlib_start( ndiag_, world_comm, intra_pool_comm, &
                     do_distr_diag_inside_bgrp_ = .FALSE. )

  call set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, inter_bgrp_comm )

  call environment_start( 'PWSCF' )

  call read_input_file( 'PW', input_file_ )

  write(*,*) '================================= Enter iosys'
  call iosys()
  write(*,*) '================================= After iosys'

  call check_stop_init()  ! required in c_bands

  write(*,*) '================================= Enter setup'
  call setup()
  write(*,*) '================================= After setup'


  write(*,*) '================================= Enter init_run'
  call init_run()
  write(*,*) '================================= After init_run'

  call my_electrons()

end program
