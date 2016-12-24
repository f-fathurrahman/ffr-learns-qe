SUBROUTINE t_import_gvecs()
  USE gvecs, ONLY : ngms, ngms_g, ngsx
  USE gvecs, ONLY : nls, nlsm
  USE gvecs, ONLY : ecuts, gcutms
  USE gvecs, ONLY : dual, doublegrid

  WRITE(*,*)
  WRITE(*,*) '-----------------------------------------------------------------'
  WRITE(*,*) 't_import_gvecs(): START'
  WRITE(*,*) '-----------------------------------------------------------------'
  WRITE(*,*)
  WRITE(*,fmt=9) 'ngms   = ', ngms
  WRITE(*,fmt=9) 'ngms_g = ', ngms_g
  WRITE(*,fmt=9) 'ngsx   = ', ngsx
  WRITE(*,*)
  WRITE(*,*) 'shape(nls)  = ', shape(nls)
  WRITE(*,*) 'shape(nlsm) = ', shape(nlsm)
  WRITE(*,*)
  WRITE(*,fmt=99) 'ecuts  = ', ecuts
  WRITE(*,fmt=99) 'gcutms = ', gcutms
  WRITE(*,*)
  WRITE(*,fmt=99) 'dual       = ', dual
  WRITE(*,*) 'doublegrid = ', doublegrid
  WRITE(*,*)
  WRITE(*,*) 't_import_gvecs(): PASSED'
  WRITE(*,*) '-----------------------------------------------------------------'

9  FORMAT(1x,A,I10)
99 FORMAT(1x,A,F18.10)

END SUBROUTINE

