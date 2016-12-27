SUBROUTINE info_gvect()
  USE gvect, ONLY : g, gg, ngm
  IMPLICIT NONE
  INTEGER :: ig

  WRITE(*,*)
  WRITE(*,*) 'INFO_GVECT'
  WRITE(*,*) '----------'
  WRITE(*,*)
  WRITE(*,*) 'ngm       = ', ngm
  WRITE(*,*) 'shape(g)  = ', shape(g)
  WRITE(*,*) 'shape(gg) = ', shape(gg)

  WRITE(*,*)
  WRITE(*,*) 'Some g-vectors: (10 first)'
  DO ig = 1, 10
    WRITE(*,'(1x,I8,3F10.5,F18.10)') ig, g(1:3,ig), gg(ig)
  ENDDO

  WRITE(*,*)
  WRITE(*,*) 'Some g-vectors: (10 last)'
  DO ig = ngm-10, ngm
    WRITE(*,'(1x,I8,3F10.5,F18.10)') ig, g(1:3,ig), gg(ig)
  ENDDO

END SUBROUTINE 

