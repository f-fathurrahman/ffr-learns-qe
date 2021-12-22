include 'prepare_all.f90'

PROGRAM main
  IMPLICIT NONE 

  CALL prepare_all()
  call test_uspp()

END PROGRAM 



SUBROUTINE test_uspp()
  use uspp_param
  use us, only: spline_ps, qrad
  use uspp
  implicit none
  ! local variables
  INTEGER :: Nupf, isp, i, Ni, j, Nj

  Nupf = size(upf,1)
  WRITE(*,*) 'Nupf = ', Nupf ! should be the same as ntyp

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


  write(*,*)
  write(*,*) 'nh = ', nh(1:Nupf)   ! number of beta functions per atomic type
  write(*,*) 'nh = ', nhm     ! max number of different beta functions per atom
  write(*,*) 'nbetam = ', nbetam  ! max number of beta functions


  write(*,*) 'lmaxkb = ', lmaxkb ! max angular momentum
  write(*,*) 'lmaxq  = ', lmaxq  ! max angular momentum + 1 for Q functions
  write(*,*) 'nvb    = ', nvb    ! number of species with Vanderbilt PPs (CPV)
  write(*,*) 'ish    = ', ish(1:Nupf)
  ! for each specie the index of the first beta 
  ! function: ish(1)=1, ish(i)=1+SUM(nh(1:i-1))


  write(*,*)
  write(*,*) 'shape(dvan) = ', shape(dvan)
  write(*,*) 'shape(deeq) = ', shape(deeq)
  write(*,*) 'shape(qq_nt) = ', shape(qq_nt)
  write(*,*) 'shape(qq_at) = ', shape(qq_at)


END SUBROUTINE

