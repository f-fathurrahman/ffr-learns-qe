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
  !CALL test_upf_uspp()

  CALL test_deeq()

  CALL mp_global_end()

END PROGRAM 


SUBROUTINE test_deeq()
  USE uspp, ONLY : deeq, indv_ijkb0, nkb
  USE constants, ONLY : e2
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE uspp_param, ONLY : nh
  IMPLICIT NONE 
  INTEGER :: isp, ia, ii, jj

  WRITE(*,*) 'nkb = ', nkb
  DO isp = 1,ntyp
    WRITE(*,*) isp, nh(isp)
  ENDDO 
  WRITE(*,*) 'shape(deeq) = ', shape(deeq)
  
  !WRITE(*,*) 'deeq = ', deeq*e2  ! convert to Ha unit

  DO isp = 1,ntyp
    IF( nh(isp) == 0 ) CYCLE 
    DO ia = 1,nat
      IF( ityp(ia) == isp ) THEN 
        WRITE(*,*) indv_ijkb0(ia) + 1
        DO ii = 1,nh(isp)
        DO jj = 1,nh(isp)
          IF( abs(deeq(ii,jj,ia,1)) > 0.d0 ) THEN 
            WRITE(*,'(4I3,F18.10)') isp, ia, ii, jj, deeq(ii,jj,ia,1)*e2
          ENDIF 
        ENDDO 
        ENDDO 
      ENDIF 
    ENDDO 
  ENDDO 

  WRITE(*,*) 'shape(indv_ijkb0) = ', shape(indv_ijkb0)
END SUBROUTINE 


!---------------------
SUBROUTINE test_gth()

  USE m_gth, ONLY : gth_p, gth_parameters

  IMPLICIT NONE 
  INTEGER :: NumGTH, igth
  TYPE(gth_parameters) :: cgth
  INTEGER :: nbeta, npr, i

  WRITE(*,*)
  WRITE(*,*) 'Calling test_gth'
  WRITE(*,*) '----------------'
  WRITE(*,*)

  NumGTH = size(gth_p)
  WRITE(*,*) 'NumGTH = ', NumGTH

  DO igth = 1, NumGTH

    WRITE(*,*)
    WRITE(*,*) 'GTH parameters ', igth
    
    cgth = gth_p(igth)  ! current GTH in this iteration
    nbeta = size(cgth%lll)
    npr = size(cgth%ipr)

    WRITE(*,'(1x,A,3I5)') 'itype, lloc, lmax = ', cgth%itype, cgth%lloc, cgth%lmax
    WRITE(*,'(1x,A,2I5)') 'nbeta, npr = ', nbeta, npr
    
    WRITE(*,*) 'List of beta functions:'
    DO i = 1, nbeta
      WRITE(*,*) i, cgth%lll(i), cgth%ipr(i)
    ENDDO 

  ENDDO

END SUBROUTINE 


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
    WRITE(*,'(1x,A,I5)') 'nh    = ', nh(iupf)   ! Is this OK for NumUPF /= ntyp ????
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
        !WRITE(*,'(1x,4I5)') ityp, ibeta, ih, indv(ih,ityp)
        IF ( ibeta == indv(ih,ityp) ) THEN 
          lm = nhtolm(ih,ityp)
          WRITE(*,'(1x,A,I5,A,I5)') 'same ibeta = ', ibeta, ' lm = ', lm
        ENDIF 
      ENDDO 
    ENDDO 
  ENDDO 


END SUBROUTINE 

