# Common usage

One way to investigate the program is by letting it read the input file
and pseudopotentials, initialize internal data structure, prepare for
SCF and then proceed to print out the variables we are interested in
or call modified subroutines that we want to study.

For example, I want to study the subroutine symvector:

```fortran
include 'prepare_all.f90'
include "my_symme.f90"

PROGRAM main
  USE kinds, ONLY: DP
  IMPLICIT NONE 

  CALL prepare_all()
  CALL test_symmetry()
END PROGRAM 

SUBROUTINE test_symmetry()
  ! .... call my_symvector here with appropriate input
END SUBROUTINE


SUBROUTINE my_symvector()
  ! ... put modified definition of symvector here
END SUBROUTINE

```


Content of prepare_all:
```fortran
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
  CALL init_run()

END SUBROUTINE 

```


API of Quantum Espresso changes very often, so the subroutine might be modified in
the recent version.