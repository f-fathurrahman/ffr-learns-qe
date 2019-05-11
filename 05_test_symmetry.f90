PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL test_symmetry()
END PROGRAM 


SUBROUTINE test_symmetry()
  USE symm_base, ONLY : nrot, fft_fact, ft, ftau, s, nsym
  USE symme, ONLY : no_rho_sym
  IMPLICIT NONE 
  INTEGER :: irot

  WRITE(*,*)
  WRITE(*,*) 'nrot = ', nrot

  WRITE(*,*)
  WRITE(*,*) 'fft_fact(1) = ', fft_fact(1)
  WRITE(*,*) 'fft_fact(2) = ', fft_fact(2)
  WRITE(*,*) 'fft_fact(3) = ', fft_fact(3)

  DO irot = 1,nrot
    WRITE(*,'(1x,A,I2,A,3F18.10)') 'ft(:,', irot, ') = ', ft(:,irot)
  ENDDO 

  DO irot = 1,nrot
    WRITE(*,'(1x,A,I2,A,3I5)') 'ftau(:,', irot, ') = ', ftau(:,irot)
  ENDDO 

  DO irot = 1,nrot
    WRITE(*,*) 'irot = ', irot
    WRITE(*,*) s(1,:,irot)
    WRITE(*,*) s(2,:,irot)
    WRITE(*,*) s(3,:,irot)
  ENDDO 

  WRITE(*,*) 'nsym = ', nsym
  WRITE(*,*) 'no_rho_sym = ', no_rho_sym


END SUBROUTINE 
