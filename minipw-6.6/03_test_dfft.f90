include 'prepare_all.f90'

PROGRAM test_dfft
  IMPLICIT NONE 

  CALL prepare_all()

  CALL info_dffts()
END PROGRAM 


SUBROUTINE info_dffts()
  
  USE fft_base, ONLY : dffts
  USE klist, ONLY : igk_k
  USE gvect, ONLY : gg
  USE cell_base, ONLY : tpiba2

  IMPLICIT NONE 

  ! This is normally taken from wvfct module
  ! Here, we set it manually
  INTEGER :: current_k
  INTEGER :: i

  current_k = 1
  
  WRITE(*,*)
  WRITE(*,*) 'nnr = ', dffts%nnr
  WRITE(*,*) 'shape of dffts%nl = ', shape(dffts%nl)

  WRITE(*,*) 'Some value of dffts%nl = '
  DO i = 1,20
    WRITE(*,*) i, dffts%nl(i)
  ENDDO 

  WRITE(*,*) 'Some value of gg = '
  DO i = 1,20
    WRITE(*,*) i, gg(i)*tpiba2
  ENDDO 

  WRITE(*,*) 'shape of igk_k = ', shape(igk_k)

  WRITE(*,*) 'Some value of igk = '
  DO i = 1,757
    WRITE(*,*) i, igk_k(i, current_k), dffts%nl(igk_k(i, current_k))
  ENDDO 


END SUBROUTINE 
