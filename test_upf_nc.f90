PROGRAM test_upf_nc
  
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

  CALL test_upf_uspp()

  CALL mp_global_end()

END PROGRAM 


SUBROUTINE test_upf_uspp()

  USE ions_base, ONLY : ntyp => nsp
  USE uspp_param, ONLY : upf, lmaxkb, nh
  USE uspp, ONLY : indv, nkb, nhtolm
  USE pseudo_types, ONLY : pseudo_upf
  
  IMPLICIT NONE 
  TYPE(pseudo_upf) :: cupf
  INTEGER :: NumUPF, iupf
  INTEGER :: ityp, ih, ibeta, lm

  WRITE(*,*)
  WRITE(*,*) 'Calling test_upf_uspp'
  WRITE(*,*) '---------------------'
  WRITE(*,*)

  WRITE(*,'(1x,A,I5)') 'lmaxkb = ', lmaxkb
  WRITE(*,*) 

  NumUPF = size(upf)
  WRITE(*,'(1x,A,I5)') 'NumUPF = ', NumUPF

  DO iupf = 1, NumUPF
    cupf = upf(iupf)
    WRITE(*,'(1x,A,I5)') 'iupf = ', iupf
    WRITE(*,'(1x,A,I5)') 'nbeta = ', cupf%nbeta
  ENDDO 

  WRITE(*,'(1x,A,I5)') 'Total number of beta functions: nkb = ', nkb
  DO ityp = 1, ntyp
    WRITE(*,*)
    WRITE(*,*) 'ityp = ', ityp
    WRITE(*,*) 'Loop over ibeta = 1 until ', upf(ityp)%nbeta
    DO ibeta = 1, upf(ityp)%nbeta
      WRITE(*,*)
      WRITE(*,*) 'ibeta = ', ibeta
      WRITE(*,*) 'Loop over ih = 1 until ', nh(ityp)
      DO ih = 1, nh(ityp)
        WRITE(*,'(1x,4I5)') ityp, ibeta, ih, indv(ih,ityp)
        IF ( ibeta == indv(ih,ityp) ) THEN 
          lm = nhtolm(ih,ityp)
          WRITE(*,'(1x,A,I5,A,I5)') 'same ibeta = ', ibeta, ' lm = ', lm
        ENDIF 
      ENDDO 
    ENDDO 
  ENDDO 


END SUBROUTINE 

