PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL test_symmetry()
END PROGRAM 


SUBROUTINE test_symmetry()
  USE symm_base, ONLY : nrot, fft_fact, ft, ftau, s, nsym
  IMPLICIT NONE 
  INTEGER :: isym

  WRITE(*,*)
  WRITE(*,*) 'nrot = ', nrot

  WRITE(*,*)
  WRITE(*,*) 'fft_fact(1) = ', fft_fact(1)
  WRITE(*,*) 'fft_fact(2) = ', fft_fact(2)
  WRITE(*,*) 'fft_fact(3) = ', fft_fact(3)

  DO isym = 1,nsym
    WRITE(*,'(1x,A,I2,A,3F18.10)') 'ft(:,', isym, ') = ', ft(:,isym)
  ENDDO 

  DO isym = 1,nsym
    WRITE(*,'(1x,A,I2,A,3I5)') 'ftau(:,', isym, ') = ', ftau(:,isym)
  ENDDO 

  DO isym = 1,nsym
    WRITE(*,*) 'irot = ', isym
    WRITE(*,*) s(1,:,isym)
    WRITE(*,*) s(2,:,isym)
    WRITE(*,*) s(3,:,isym)
  ENDDO 

  WRITE(*,*) 'nsym = ', nsym


END SUBROUTINE 
