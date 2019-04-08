PROGRAM new_vlocal

  IMPLICIT NONE 

  CALL prepare_all()

  CALL test_vlocal()

END PROGRAM


SUBROUTINE test_vlocal()

  USE vlocal, ONLY : vloc, strf
  USE gvect, ONLY : ngl, ngm
  USE ions_base, ONLY : ntyp => nsp
  USE fft_base,  ONLY : dfftp
  USE gvect, ONLY : ngm
  USE gvect, ONLY : igtongl, gcutm
  USE scf, ONLY : v_of_0
  USE fft_interfaces, ONLY : invfft
  USE cell_base, ONLY : tpiba2

  IMPLICIT NONE 

  REAL(8), PARAMETER :: PI=4.d0*atan(1.d0)
  INTEGER :: nt, ng
  
  COMPLEX(8), ALLOCATABLE :: aux(:)
  REAL(8), ALLOCATABLE :: vloc_r(:)

  ALLOCATE( aux(dfftp%nnr) )
  ALLOCATE( vloc_r(dfftp%nnr) )

  aux(:) = (0.d0,0.d0)
  vloc_r(:) = 0.d0


  WRITE(*,*) 'dfftp%nnr = ', dfftp%nnr

  WRITE(*,*) 'ngl = ', ngl  ! no. of g-shells
  WRITE(*,*) 'ngm = ', ngm
  WRITE(*,*) 'size vloc = ', size(vloc)

  DO nt = 1, ntyp
      DO ng = 1, ngm
          aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vloc(igtongl(ng),nt)*strf(ng, nt)
      END DO
  END DO

  WRITE(*,*) 'v_of_0 = ', v_of_0

  CALL invfft('Rho', aux, dfftp)
  vloc_r = dble(aux)

  WRITE(*,*) 'sum vloc_r in Ha', sum(vloc_r)/2.d0
  WRITE(*,*) 'vloc_r', vloc_r(1:5)

  WRITE(*,*) 'gcutm = ', gcutm*tpiba2
  WRITE(*,*) 'PI    = ', PI
  WRITE(*,*) 4.d0*15.d0/2.d0/PI**2

END SUBROUTINE 
