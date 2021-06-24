include 'prepare_all.f90'

PROGRAM main
  IMPLICIT NONE 

  CALL prepare_all()
  call test_uspp()

END PROGRAM 



SUBROUTINE test_uspp()
  use uspp_param, only: upf
  use ener, only: etxcc
  use us, only: spline_ps, qrad
  ! local variables
  INTEGER :: Nupf, isp, i, Ni, j, Nj

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
    write(*,*) 'shape dion = ', shape(upf(isp)%dion)
    Ni = size(upf(isp)%dion,1)
    Nj = size(upf(isp)%dion,2)
    write(*,*) 'dion = '
    do i = 1,Ni
      do j = 1,Nj
        !write(*,'(1x,F18.5)',advance='no') upf(isp)%dion(i,j)
        write(*,'(1x,ES18.10)',advance='no') upf(isp)%dion(i,j)
      enddo
      write(*,*)
    enddo
    write(*,*) 'nlcc = ', upf(isp)%nlcc
    write(*,*) 'kkbeta = ', upf(isp)%kkbeta
    write(*,*) 'kbeta = ', upf(isp)%kbeta
  enddo

  write(*,*) 'spline_ps = ', spline_ps

  write(*,*) 'shape qrad: ', shape(qrad)
  write(*,*) 'qrad: ', qrad(1,1,1,1)

  write(*,*) 'etxcc = ', etxcc

END SUBROUTINE

