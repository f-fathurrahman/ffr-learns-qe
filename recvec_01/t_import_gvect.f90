! ffr, 24 June 2016

SUBROUTINE t_import_gvect()
  USE gvect, ONLY : ngm, ngm_g, ngl, ngmx, ecutrho, gcutm
  USE gvect, ONLY : nl, nlm
  USE gvect, ONLY : gstart
  USE gvect, ONLY : gg
  USE gvect, ONLY : gl, igtongl
  USE gvect, ONLY : g
  USE gvect, ONLY : mill
  USE gvect, ONLY : ig_l2g
  USE gvect, ONLY : sortedig_l2g
  USE gvect, ONLY : mill_g
  USE gvect, ONLY : eigts1, eigts2, eigts3


  WRITE(*,*) '-----------------------------------------------------------------'
  WRITE(*,*) 't_import_gvect(): START'
  WRITE(*,*) '-----------------------------------------------------------------'
  WRITE(*,*)
  WRITE(*,fmt=9) 'ngm     = ', ngm
  WRITE(*,fmt=9) 'ngm_g   = ', ngm_g
  WRITE(*,fmt=9) 'ngl     = ', ngl
  WRITE(*,fmt=9) 'ngmx    = ', ngmx
  WRITE(*,fmt=99) 'ecutrho = ', ecutrho
  WRITE(*,fmt=99) 'gcutm   = ', gcutm
  WRITE(*,*)
  WRITE(*,*) 'shape(nl)  = ', shape(nl)
  WRITE(*,*) 'shape(nlm) = ', shape(nlm)
  WRITE(*,*)
  WRITE(*,fmt=9) 'gstart = ', gstart
  WRITE(*,*)
  WRITE(*,*) 'shape(gg) = ', shape(gg)
  WRITE(*,*)
  WRITE(*,*) 'shape(gl)      = ', shape(gl)
  WRITE(*,*) 'shape(igtongl) = ', shape(igtongl)
  WRITE(*,*)
  WRITE(*,*) 'shape(g)    = ', shape(g)
  WRITE(*,*) 'shape(mill) = ', shape(mill)
  WRITE(*,*)
  WRITE(*,*) 'shape(ig_l2g)     = ', shape(ig_l2g)
  WRITE(*,*) 'shape(sorted_l2g) = ', shape(sortedig_l2g)
  WRITE(*,*)
  WRITE(*,*) 'shape(mill_g) = ', shape(mill_g)
  WRITE(*,*)
  WRITE(*,*) 'shape(eigts1) = ', shape(eigts1)
  WRITE(*,*) 'shape(eigts2) = ', shape(eigts2)
  WRITE(*,*) 'shape(eigts3) = ', shape(eigts3)
  WRITE(*,*)
  WRITE(*,*) 't_import_gvect(): PASSED'
  WRITE(*,*) '-----------------------------------------------------------------'

9  FORMAT(1x,A,I10)
99 FORMAT(1x,A,F18.10)
END SUBROUTINE

