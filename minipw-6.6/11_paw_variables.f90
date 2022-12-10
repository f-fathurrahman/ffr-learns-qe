PROGRAM main
  USE kinds, ONLY: DP
  IMPLICIT NONE 

  CALL prepare_all()
  CALL test_paw_variables()
  CALL test_uspp_paw()
  
END PROGRAM main


SUBROUTINE test_uspp_paw()
  USE uspp_param, ONLY: upf
  USE us, ONLY: spline_ps, qrad
  ! local variables
  INTEGER :: Nupf, isp, i, Ni, j, Nj

  WRITE(*,*)
  WRITE(*,*) '----------------------------------------'
  WRITE(*,*) 'Some fields of UPF (per atomic species):'
  WRITE(*,*) '----------------------------------------'
  
  Nupf = size(upf,1)
  WRITE(*,*) 'Nupf = ', Nupf

  DO isp = 1,Nupf
    WRITE(*,*)
    WRITE(*,*) 'isp = ', isp
    WRITE(*,*) 'typ = ', upf(isp)%typ
    WRITE(*,*) 'ecutwfc = ', upf(isp)%ecutwfc
    WRITE(*,*) 'ecutrho = ', upf(isp)%ecutrho
    WRITE(*,*) 'tvanp = ', upf(isp)%tvanp
    WRITE(*,*) 'is_multiproj = ', upf(isp)%is_multiproj
    WRITE(*,*) 'nbeta = ', upf(isp)%nbeta
    WRITE(*,*) 'lll = ', upf(isp)%lll
    WRITE(*,*) 'shape dion = ', SHAPE(upf(isp)%dion)
    Ni = SIZE(upf(isp)%dion,1)
    Nj = SIZE(upf(isp)%dion,2)
    WRITE(*,*) 'dion = '
    DO i = 1,Ni
      DO j = 1,Nj
        !write(*,'(1x,F18.5)',advance='no') upf(isp)%dion(i,j)
        WRITE(*,'(1x,ES18.10)',advance='no') upf(isp)%dion(i,j)
      ENDDO
      WRITE(*,*)
    ENDDO  
  ENDDO

  WRITE(*,*) 'spline_ps = ', spline_ps
  WRITE(*,*) 'shape qrad: ', SHAPE(qrad)
  WRITE(*,*) 'qrad: ', qrad(1,1,1,1)

END SUBROUTINE test_uspp_paw


SUBROUTINE test_paw_variables()
  USE paw_variables, ONLY: total_core_energy, okpaw
  USE uspp, ONLY: okvan

  WRITE(*,*)
  WRITE(*,*) '---------------------------------------'
  WRITE(*,*) 'Some variables in paw_variables module:'
  WRITE(*,*) '---------------------------------------'
  
  WRITE(*,*) 'okpaw = ', okpaw
  WRITE(*,*) 'okvan = ', okvan
  WRITE(*,*) 'total_core_energy = ', total_core_energy

END SUBROUTINE test_paw_variables

