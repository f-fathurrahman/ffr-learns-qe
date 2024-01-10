include 'prepare_all.f90'

PROGRAM new_gvectors

  IMPLICIT NONE 

  CALL prepare_all()

  CALL info_gvectors()

END PROGRAM 

SUBROUTINE info_gvectors()

  USE cell_base, ONLY : tpiba, alat, tpiba2
  USE gvect, ONLY : ngm, ecutrho, gcutm
  USE gvecs, ONLY: dual

  WRITE(*,*)
  WRITE(*,*) '-------------'
  WRITE(*,*) 'GVectors info'
  WRITE(*,*) '-------------'
  WRITE(*,*)
  WRITE(*,*) 'Number gvectors:', Ngm
  WRITE(*,*)
  WRITE(*,*) 'ecutrho = ', ecutrho
  WRITE(*,*) 'gcutm = ', gcutm
  WRITE(*,*) 'dual = ', dual
  WRITE(*,*) 'gcutm*tpiba2 = ', gcutm*tpiba2
  WRITE(*,*)
  WRITE(*,*) 'alat = ', alat
  WRITE(*,*) 'tpiba = ', tpiba


END SUBROUTINE 


