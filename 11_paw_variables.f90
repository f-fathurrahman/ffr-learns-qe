PROGRAM main
  USE kinds, ONLY: DP
  IMPLICIT NONE 

  CALL prepare_all()
  !CALL test_paw_variables()
  call test_uspp_paw()

END PROGRAM 



SUBROUTINE test_uspp_paw()
  USE uspp_param, ONLY: upf
  ! local variables
  INTEGER :: Nupf, isp

  Nupf = size(upf,1)
  WRITE(*,*) 'Nupf = ', Nupf

  do isp = 1,Nupf
    write(*,*)
    write(*,*) 'isp = ', isp
    write(*,*) 'typ = ', upf(isp)%typ
    write(*,*) 'ecutwfc = ', upf(isp)%ecutwfc
    write(*,*) 'ecutrho = ', upf(isp)%ecutrho
    write(*,*) 'tvanp = ', upf(isp)%tvanp
    write(*,*) 'is_multiproj = ', upf(isp)%is_multiproj
    write(*,*) 'nbeta = ', upf(isp)%nbeta
    write(*,*) 'lll = ', upf(isp)%lll
    write(*,*) 'shape dion = ', size(upf(isp)%dion)
    write(*,*) 'dion = ', upf(isp)%dion
  enddo

END SUBROUTINE



SUBROUTINE test_paw_variables()
  USE paw_variables, ONLY: total_core_energy, okpaw
  USE uspp, ONLY: okvan

  WRITE(*,*) 'okpaw = ', okpaw
  WRITE(*,*) 'okvan = ', okvan
  WRITE(*,*) 'total_core_energy = ', total_core_energy

END SUBROUTINE