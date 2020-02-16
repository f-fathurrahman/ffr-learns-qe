include "my_symm_base.f90"


SUBROUTINE test_my_symm_base()
  
  USE my_symm_base, ONLY : set_sym_bl
  IMPLICIT NONE 
  
  CALL set_sym_bl()

END SUBROUTINE 


PROGRAM main

  IMPLICIT NONE 
  
  CALL prepare_all()

  CALL test_my_symm_base()

END PROGRAM 



