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

  CALL test_gth()

  CALL mp_global_end()

END PROGRAM 



!---------------------
SUBROUTINE test_gth()

  USE m_gth, ONLY : gth_p, gth_parameters
  USE uspp_param, ONLY : upf
  USE pseudo_types, ONLY : pseudo_upf

  IMPLICIT NONE 
  INTEGER :: N_GTH, N_UPF
  INTEGER :: igth, iupf
  TYPE(gth_parameters) :: cgth
  TYPE(pseudo_upf) :: cupf

  N_GTH = size(gth_p)
  WRITE(*,*) 'N_GTH = ', N_GTH

  DO igth = 1, N_GTH

    WRITE(*,*)
    WRITE(*,*) 'GTH parameters ', igth
    
    cgth = gth_p(igth)  ! current GTH in this iteration

    WRITE(*,*) 'itype, lloc, lmax = ', cgth%itype, cgth%lloc, cgth%lmax
    WRITE(*,*)
  ENDDO

  N_UPF = size(upf)
  WRITE(*,*) 'N_UPF = ', N_UPF
  DO iupf = 1, N_UPF
    cupf = upf(iupf)
    WRITE(*,*) cupf%nbeta
  ENDDO 

END SUBROUTINE 

