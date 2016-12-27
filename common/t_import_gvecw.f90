SUBROUTINE t_import_gvecw()
  USE gvecw, ONLY : ngw, ngw_g, ngwx, ecutwfc, gcutw, ekcut, gkcut
  IMPLICIT NONE

  WRITE(*,*)
  WRITE(*,*) '-----------------------------------------------------------------'
  WRITE(*,*) 't_import_gvecw(): START'
  WRITE(*,*) '-----------------------------------------------------------------'
  WRITE(*,*) 
  WRITE(*,fmt=9) 'ngw     = ', ngw
  WRITE(*,fmt=9) 'ngw_g   = ', ngw_g
  WRITE(*,fmt=9) 'ngwx    = ', ngwx
  WRITE(*,fmt=99) 'ecutwfc = ', ecutwfc
  WRITE(*,fmt=99) 'gcutw   = ', gcutw
  WRITE(*,fmt=99) 'ekcut   = ', ekcut
  WRITE(*,fmt=99) 'gkcut   = ', gkcut
  WRITE(*,*)
  WRITE(*,*) 't_import_gvecw(): PASSED'
  WRITE(*,*) '-----------------------------------------------------------------'

9  FORMAT(1x,A,I10)
99 FORMAT(1x,A,F18.10)
END SUBROUTINE

