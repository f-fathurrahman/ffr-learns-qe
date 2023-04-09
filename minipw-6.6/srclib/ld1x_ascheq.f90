!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
SUBROUTINE ascheq(nn,lam,e,mesh,grid,vpot,ze2,thresh0,y,nstop)
!---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation for
  !  bound states in a local potential.
  !  thresh determines the absolute accuracy for the eigenvalue
  !
  USE kinds, ONLY : DP
  USE radial_grids, ONLY: radial_grid_type, series
  !
  IMPLICIT NONE
  !
  TYPE(radial_grid_type), INTENT(in) :: grid
  INTEGER :: mesh,lam, ierr
  INTEGER:: nn,nstop,maxter,iter,l1,i,ik,ncross,n, &
       nstart,ns,n2,nst2,ndcr
  REAL(DP) :: ze2,ddx12,eup,elw,b0e,ymx,rap,rstart,di,expn,  &
       c1,c2,fe,a0,a1,a2,sum0,f2, ss, sqlhf,f0,f1,dfe,de,eps, &
       yln,xp,sum1,xl1,x4l6
  REAL(DP):: vpot(mesh), y(mesh)
  REAL(DP), ALLOCATABLE:: c(:), el(:), f(:)
  REAL(DP):: b(0:3), e, thresh0, thresh
  DATA maxter/50/
  

  WRITE(*,*)
  WRITE(*,*) 'Solving radial Schroedinger equation'
  WRITE(*,*) 'ze2       = ', ze2
  WRITE(*,*) 'nn        = ', nn
  WRITE(*,*) 'll        = ', lam
  WRITE(*,*) 'Trial Enl = ', e
  WRITE(*,*) 'sum vpot  = ', sum(vpot)
  WRITE(*,*)


  !
  ! set up constants and initialize
  !
  ALLOCATE( c(mesh), stat=ierr )
  ALLOCATE( f(mesh), stat=ierr )
  ALLOCATE( el(mesh), stat=ierr )

  thresh = thresh0
  IF( e < -5.e+2 ) thresh = thresh0*10.0_DP
  iter = 0
  ddx12 = grid%dx * grid%dx / 12.0_dp
  l1 = lam + 1
  sqlhf = (DBLE(lam) + 0.5_dp)**2
  
  ndcr = nn - lam - 1
  ! number of node crossings
  
  !
  ! set initial lower and upper bounds to the eigenvalue
  !
  eup = vpot(mesh) + sqlhf/grid%r2(mesh)
  elw = eup
  DO i = 1,mesh
    elw = min(elw, vpot(i) + sqlhf/grid%r2(i))
  ENDDO

  nstop = 200

  IF( eup == elw ) GOTO 900
  
  IF( e > eup ) e = 0.9_DP*eup + 0.1_DP*elw

  IF( e < elw ) e = 0.9_DP*elw + 0.1_DP*eup

  WRITE(*,*)
  WRITE(*,*) 'Initial lower and upper bounds to the eigenvalue (in Ha)'
  WRITE(*,*)
  WRITE(*,'(1x,A,F18.10)') 'e   = ', e/2
  WRITE(*,'(1x,A,F18.10)') 'eup = ', eup/2
  WRITE(*,'(1x,A,F18.10)') 'elw = ', elw/2

  !
  ! series developement of the potential near the origin
  !
  DO i = 1,4
    y(i) = vpot(i) - ze2/grid%r(i)
  ENDDO
  
  CALL series(y, grid%r, grid%r2, b)

  WRITE(*,*)
  WRITE(*,*) 'Near origin: r, vpot, yi, b (in Ha)'
  WRITE(*,*)
  DO i = 1,4
    WRITE(*,'(1x,4ES18.10)') grid%r(i), vpot(i)/2.d0, y(i)/2.d0, b(i-1)/2.d0
    ! Note that b has shape of (0:3)
  ENDDO
  
  !stop 'ffr ascheq line 66'

  !
  ! Loop starts here
  !

300 CONTINUE

  iter = iter + 1

  
  WRITE(*,*) '***** ======================='
  WRITE(*,*) '***** Starting iterSch = ', iter
  WRITE(*,*) '***** ======================='

  nstop = 300
  IF( iter > maxter ) GOTO 900

  !
  ! set up the f-function and determine the position of its last
  ! change of sign
  ! f < 0 (approximatively) means classically allowed   region
  ! f > 0         "           "        "      forbidden   "
  f(1) = ddx12*( grid%r2(1)*(vpot(1) - e) + sqlhf )

  DO i = 2,mesh
    f(i) = ddx12*( grid%r2(i)*( vpot(i) - e ) + sqlhf)
    IF( f(i) /= SIGN(f(i), f(i-1)) ) THEN
      ik = i
    ENDIF
  ENDDO
  nstop = 302

  IF( ik >= mesh-2 ) GOTO 900

  DO i = 1,mesh
    f(i) = 1.0_dp - f(i)
  ENDDO

  !DO i = 1,4
  !  WRITE(*,'(1x,A,F18.10)') 'f = ', f(i)
  !ENDDO
  !WRITE(*,'(1x,A,F18.10)') 'ddx12 = ', ddx12
 
  !
  y(:) = 0.0_dp
  !
  ! determination of the wave-function in the first two points by
  ! series developement
  ! 
  xl1 = lam + 1.0_DP
  x4l6 = 4.0_dp*lam + 6.0_dp
  b0e = b(0) - e
  c1 = 0.5_dp*ze2/xl1
  c2 = (c1*ze2 + b0e)/x4l6

  !WRITE(*,*)
  !WRITE(*,'(1x,A,F18.10)') 'xl1 = ', xl1
  !WRITE(*,'(1x,A,F18.10)') 'x4l6 = ', x4l6
  !WRITE(*,'(1x,A,F18.10)') 'b0e = ', b0e/2 ! Ha
  !WRITE(*,'(1x,A,F18.10)') 'c1 = ', c1
  !WRITE(*,'(1x,A,F18.10)') 'c2 = ', c2/2 ! Ha

  !WRITE(*,*)
  !WRITE(*,*) 'before start_scheq: e = ', e*0.5
  !
  CALL start_scheq( lam, e, b, grid, ze2, y )

  write(*,*)
  WRITE(*,*) 'Two initial points from series development:'
  WRITE(*,'(1x,A,F18.10)') 'y(1) = ', y(1) !*0.5d0
  WRITE(*,'(1x,A,F18.10)') 'y(2) = ', y(2) !*0.5d0
  !stop 'ffr ascheq line 110'

  !
  !  start outward integration and count number of crossings
  !
  ncross = 0
  ymx = 0.0_dp
  WRITE(*,*)
  WRITE(*,*) 'Begin outward integration'
  WRITE(*,*) 'until ik = ', ik
  !WRITE(*,*) 'f = ', f(1:4)
  DO n = 2,ik-1
    y(n+1) = ( (12.0_dp - 10.0_dp*f(n) )*y(n) - f(n-1)*y(n-1))/f(n+1)
    IF( y(n) /= SIGN(y(n),y(n+1)) ) THEN
      ncross = ncross + 1
    ENDIF
    ymx = MAX(ymx, ABS(y(n+1)))
  ENDDO
  
  WRITE(*,'(1x,A,F18.10)') 'ymx = ', ymx
  WRITE(*,*) 'ncross = ', ncross
  !stop 'ffr ascheq line 130'


  !
  ! matching radius has been reached going out. if ncross is not
  ! equal to ndcr, modify the trial eigenvalue.
  !
  IF( ndcr < ncross ) THEN
    ! too many crossings. e is an upper bound to the true eigen-
    ! value. increase abs(e)
    eup = e
    rap = (DBLE(ncross+l1)/DBLE(nn))**2
    e = (e - vpot(mesh))*rap + vpot(mesh)
    IF(e < elw) THEN
      e = 0.9_dp*elw + 0.1_dp*eup
    ENDIF
    WRITE(*,*) 'Too many crossings. e is increased to ', e
    GOTO 300

  ELSEIF( ndcr > ncross ) THEN
    !
    ! too few crossings. e is a lower bound to the true eigenvalue. decrease abs(e)
    !
    elw = e
    rap = ( DBLE(ncross + l1)/DBLE(nn) )**2
    e = (e - vpot(mesh))*rap + vpot(mesh)
    IF(e > eup) THEN
      e = 0.9_dp*eup + 0.1_dp*elw
    ENDIF
    WRITE(*,*) 'Too few crossings. e is decreased to ', e
    GOTO 300
  ENDIF

  !stop 'ffr ascheq 179'

  !
  ! prepare inward integration
  ! charlotte froese can j phys 41,1895(1963)
  !
  ! start at  min( rmax, 10*rmatch )
  !
  nstart = mesh
  ns = 10
  rstart = ns*grid%r(ik)
  IF(rstart < grid%r(mesh)) THEN
    DO i=ik,mesh
      nstart = i
      IF(grid%r(i) >= rstart) GOTO 403
    ENDDO
403 nstart = nstart/2
    nstart = 2*nstart + 1
  ENDIF

  WRITE(*,*)
  WRITE(*,*) 'Inward integration from nstart = ', nstart
  !stop 'ffr ascheq 201'

  !
  ! set up a, l, and c vectors
  !
  n = ik + 1
  el(n) = 10.0_dp*f(n) - 12.0_dp
  c(n) = -f(ik)*y(ik)
  n2 = ik + 2
  DO n = n2,nstart
    di = 10.0_dp*f(n) - 12.0_dp
    el(n) = di-f(n)*f(n-1)/el(n-1)
    c(n) = -c(n-1)*f(n-1)/el(n-1)
  ENDDO

  WRITE(*,*) 'el(nstart) = ', el(nstart)
  WRITE(*,*) 'c(nstart) = ', c(nstart)
  !stop 'ffr ascheq 201'


  !
  ! start inward integration by the froese's tail procedure
  !
  expn = exp(-sqrt(12.0_dp*abs(1.0_dp - f(nstart-1))))
  y(nstart-1) = c(nstart-1)/( el(nstart-1) + f(nstart)*expn)
  y(nstart) = expn*y(nstart-1)

  DO n = nstart-2,ik+1,-1
    y(n) = ( c(n) - f(n+1)*y(n+1) )/el(n)
  ENDDO
  
  WRITE(*,'(1x,A,F18.10)') 'y = ', y(ik+1)
  !stop 'ffr ascheq 231'


  !
  ! if necessary, improve the trial eigenvalue by the cooley's
  ! procedure. jw cooley math of comp 15,363(1961)
  !
  fe = ( 12.0_dp - 10.0_dp*f(ik) ) * y(ik) - f(ik-1)*y(ik-1) - f(ik+1)*y(ik+1)
  WRITE(*,'(1x,A,F18.10)') 'fe = ', fe
  !stop 'ffr ascheq 240'
  

  !
  ! calculate the normalization
  !
  IF( ymx >= 1.0e10_dp ) THEN
    DO i = 1,mesh
      y(i) = y(i)/ymx
    ENDDO  
  ENDIF
  

  WRITE(*,'(1x,A,F18.10)') 'ymx = ', ymx
  WRITE(*,'(1x,A,F18.10)') 'y[1] = ', y(1)
  WRITE(*,'(1x,A,F18.10)') 'y[Nrmesh] = ', y(mesh)
  !stop 'ffr ascheq 252'


  a0 = 1.0_dp/DBLE(2*lam + 3)
  a1 = c1/DBLE(lam + 2)
  a2 = (c1*c1 + c2 + c2)/DBLE(2*lam+5)
  sum0 = (a0 + grid%r(1)*(a1 + grid%r(1)*a2))*grid%r(1)**(2*lam + 3)
  nst2 = nstart - 2
  f2 = grid%r2(1)*y(1)*y(1)
  ss = grid%r(1)*f2/DBLE(2*l1 + 1)
  DO n = 1,nst2,2
    f0 = f2
    f1 = grid%r2(n+1)*y(n+1)*y(n+1)
    f2 = grid%r2(n+2)*y(n+2)*y(n+2)
    ss = ss + f0 + f2 + 4.0_DP*f1
  ENDDO
   
  write(*,'(1x,A,F18.10)') 'ss = ', ss
  !stop 'ffr ascheq 269'

  ss = sum0 + grid%dx*ss/3.0_dp
  dfe = -y(ik)*f(ik)/grid%dx/ss
  de = -fe*dfe
  eps = abs(de/e)

  WRITE(*,'(1x,A,F18.10)') 'ss = ', ss
  WRITE(*,'(1x,A,F18.10)') 'de = ', de*0.5
  !if(iter == 4) then
  !  stop 'ffr ascheq 277'
  !ENDIF


  write(*,*) 'iter = ', iter, ' e = ', e/2, ' de = ', abs(de)/2
  if(abs(de).lt.thresh) then
    write(*,*) 'GOTO 600 here ......'
    go to 600
  ENDIF
  
  if(eps .gt. 0.25_dp) then
    write(*,*) 'de is updated here 295 ...'
    de=0.25_dp*de/eps
  ENDIF
  
  if(de.gt.0.0_dp) then
    write(*,*) 'elw is updated here 300 ...'
    elw=e
  ENDIF
  
  if(de.lt.0.0_dp) then
    write(*,*) 'eup is updated here 305 ...'
    eup=e
  ENDIF

  e = e + de
  write(*,*) 'de = ', de/2
  write(*,*) 'e = ', e/2

  if(e .gt. eup) then
    write(*,*) 'e is updated here 312 ...'
    e = 0.9_dp*eup+0.1_dp*elw
  ENDIF
  
  IF(e < elw) THEN
    WRITE(*,*) 'e is updated here 317 ...'
    e = 0.9_dp*elw + 0.1_dp*eup
  ENDIF

  write(*,*) 'New data'
  write(*,'(1x,A,F18.10)') 'e = ', e/2
  write(*,'(1x,A,F18.10)') 'eup = ', eup/2

  WRITE(*,'(1x,A,F18.10)') 'elw = ', elw/2
  
  IF( iter < maxter) GOTO 300
  
  nstop=50

600 CONTINUE

  write(*,*)
  write(*,*) 'exit loop'
  write(*,*) 'iter = ', iter
  write(*,*) 'e = ', e/2

  !stop 'ffr ascheq 281'

  !
  !  normalize the eigenfunction and exit
  !
  DO n=nstart,mesh-1
    y(n+1) = 0.0_dp
    ! FIXME: is this OK?
    IF( y(n) == 0.0_dp ) GOTO 601
    yln = LOG(ABS(y(n)))
    xp = -SQRT(12.0_dp*ABS(1.0_dp-f(n)))
    expn = yln + xp
    IF(expn < -80.0_dp) GOTO 601
    y(n+1) = SIGN(EXP(expn),y(n))
601 CONTINUE
  ENDDO
   
  sum1 = 0.0_dp
  DO n = nstart,mesh-2,2
    f0 = f2
    f1 = grid%r2(n+1) * y(n+1) * y(n+1)
    f2 = grid%r2(n+2) * y(n+2) * y(n+2)
    sum1 = sum1 + f0 + f2 + 4.0_dp*f1
  ENDDO
   

  ss = ss + grid%dx*sum1/3.0_dp
  ss = SQRT(ss)
  WRITE(*,*) 'ss for normalization = ', ss

  DO n=1,mesh
    y(n) = grid%sqr(n)*y(n)/ss
  ENDDO
  
  IF(nstop .LT. 100) GOTO 900

  nstop = 0

  DEALLOCATE(el)
  DEALLOCATE(f)
  DEALLOCATE(c)

  RETURN

  !
  !  error exit
  !
  ! 900  write(6,9000) nstop,nn,lam,elw,eup
  ! 9000 format(5x,'error in ascheq: nstop =',i4,'. n l =',2i3,/ &
  !     & 5x,'elw =',f15.10,' eup =',f15.10)

900 CONTINUE

  DEALLOCATE(el)
  DEALLOCATE(f)
  DEALLOCATE(c)
  RETURN

END SUBROUTINE ascheq

