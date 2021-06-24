include 'prepare_all.f90'

PROGRAM new_gvectors

  IMPLICIT NONE 

  CALL prepare_all()

  CALL info_gvectors()

END PROGRAM 

SUBROUTINE info_gvectors()

  USE cell_base, ONLY : tpiba, alat
  USE gvect, ONLY : ngm, ecutrho, gcutm

  WRITE(*,*)
  WRITE(*,*) '-------------'
  WRITE(*,*) 'GVectors info'
  WRITE(*,*) '-------------'
  WRITE(*,*)
  WRITE(*,*) 'Number gvectors:', Ngm
  WRITE(*,*)
  WRITE(*,*) 'ecutrho = ', ecutrho
  WRITE(*,*) 'gcutm = ', gcutm
  WRITE(*,*)
  WRITE(*,*) 'alat = ', alat
  WRITE(*,*) 'tpiba = ', tpiba


END SUBROUTINE 


