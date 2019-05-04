PROGRAM main

  IMPLICIT NONE 

  CALL prepare_all()

  CALL my_scf()

END PROGRAM 


SUBROUTINE my_scf()

  IMPLICIT NONE 

  WRITE(*,*)
  WRITE(*,*) 'my_scf is starting'
  WRITE(*,*)



  WRITE(*,*)
  WRITE(*,*) 'my_scf is finished'
  WRITE(*,*)

END SUBROUTINE 

