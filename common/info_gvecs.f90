SUBROUTINE info_gvecs()
  USE gvecs, ONLY : ngms, ngms_g, ngsx
  IMPLICIT NONE

  WRITE(*,*)
  WRITE(*,*) 'INFO_GVECS'
  WRITE(*,*) '----------'
  WRITE(*,*)
  WRITE(*,*) 'ngms   = ', ngms
  WRITE(*,*) 'ngms_g = ', ngms_g
  WRITE(*,*) 'ngsx   = ', ngsx

END SUBROUTINE 

