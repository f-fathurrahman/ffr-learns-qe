include 'prepare_all.f90'

PROGRAM test_dfft
  IMPLICIT NONE 

  CALL prepare_all()

  CALL info_scf()
END PROGRAM 


SUBROUTINE info_scf()
  
  USE scf, ONLY : v_of_0

  IMPLICIT NONE 
  
  WRITE(*,*)
  WRITE(*,*) 'v_of_0 = ', v_of_0
  WRITE(*,*)

END SUBROUTINE 
