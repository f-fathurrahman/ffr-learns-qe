SUBROUTINE my_lschps( mode, z, eps, grid, nin, n, l, e, v, u, nstop )
  !
  ! integrates radial pauli-type scalar-relativistic equation
  ! on a logarithmic grid
  ! modified routine to be used in finding norm-conserving
  ! pseudopotential
  !
  ! on input:
  !   mode = 1 find energy and wavefunction of bound states,
  !            scalar-relativistic (all-electron)
  !   mode = 2 find energy and wavefunction of bound state,
  !            nonrelativistic (pseudopotentials)
  !   mode = 3 fixed-energy calculation, for logarithmic derivatives
  !   mode = 4 find energy which produces a specified logarithmic
  !            derivative (nonrelativistic, pseudopotentials)
  !   mode = 5 is for pseudopotential to produce wavefunction beyond
  !            radius used for pseudopotential construction
  !   z    = atomic number
  !   eps  = convergence factor: eiganvalue is considered converged if
  !          the correction to eigenvalue is smaller in magnitude than
  !          eps times the magnitude of the current guess
  !   grid = structure containing radial grid information
  !   l, n = main and angular quantum numbers
  !   e    = starting estimate of the energy (mode=1,2)
  !          fixed energy at which the wavefctn is calculated (mode=3,4)
  !   v(i) = self-consistent potential
  !   nin  = integration up to r(nin) (mode=3,4,5)
  !
  ! on output:
  !   e    = final energy (mode=1,2)
  !   u(i) = radial wavefunction (defined as the radial part of the wavefct
  !          multiplied by r)
  !   nstop= 0 if regular termination, 1 otherwise
  !   nin  = last grid point for which the wavefct is calculated (mode=1,2)
  !
  USE kinds, ONLY : DP
  USE radial_grids, ONLY: radial_grid_type
  USE ld1inc, ONLY : cau_fact
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER, INTENT (in) :: mode, n, l
  REAL(DP), INTENT(in) :: z, eps
  TYPE (radial_grid_type), INTENT(in) :: grid
  REAL(DP), INTENT(in) :: v(grid%mesh)
  INTEGER, INTENT(inout) :: nin
  REAL(DP), INTENT(inout) :: e
  INTEGER, INTENT(out) :: nstop
  REAL(DP), INTENT(out) :: u(grid%mesh)
  !
  ! local variables
  !
  INTEGER, PARAMETER :: maxter=60
  REAL(DP), EXTERNAL:: my_lschps_aei, my_lschps_aeo, my_lschps_aii, my_lschps_aio
  ! arrays  used as work space
  REAL(DP),ALLOCATABLE :: up(:),upp(:),cf(:),dv(:),fr(:),frp(:)
  REAL(DP):: al, als, cn
  real(DP):: de, emax, emin
  real(DP):: fss, gamma, ro, sc
  real(DP):: sls, sn, uld, uout,  upin, upout
  real(DP):: xkap
  INTEGER:: i, it, mmax, n_it, node, mch, ierr
  !
  !
  nstop=0
  al   = grid%dx
  mmax = grid%mesh

  ALLOCATE(up(mmax), stat=ierr)
  ALLOCATE(upp(mmax), stat=ierr)
  ALLOCATE(cf(mmax), stat=ierr)
  ALLOCATE(dv(mmax), stat=ierr)
  ALLOCATE(fr(mmax), stat=ierr)
  ALLOCATE(frp(mmax), stat=ierr)

  uld=0.0_dp
  !
  !
  IF(mode == 1 .or. mode == 3) THEN
     !     relativistic calculation
     !     fss=(1.0_dp/137.036_dp)**2
     fss=(1.0_dp/cau_fact)**2
     IF(l == 0) THEN
        gamma=sqrt(1.0_dp-fss*z**2)
     ELSE
        gamma=(l*sqrt(l**2-fss*z**2) + &
             (l+1)*sqrt((l+1)**2-fss*z**2))/(2*l+1)
     ENDIF
  ELSE
     !     non-relativistic calculation
     fss=1.0e-20_dp
     gamma=l+1
  ENDIF
  !
  sls=l*(l+1)
  !
  ! emin, emax = estimated bounds for e
  !
  IF(mode == 1 .or. mode == 2) THEN
     emax=v(mmax)+sls/grid%r(mmax)**2
     emin=0.0_dp
     DO i=1,mmax
        emin=min(emin,v(i)+sls/grid%r(i)**2)
     ENDDO
     IF(e > emax) e=1.25_dp*emax
     IF(e < emin) e=0.75_dp*emin
     IF(e > emax) e=0.5_dp*(emax+emin)
  ELSEIF(mode == 4) THEN
     emax=e + 10.0_dp
     emin=e - 10.0_dp
  ENDIF
  !
  DO i=1,4
     u(i)=0.0_dp
     up(i)=0.0_dp
     upp(i)=0.0_dp
  ENDDO
  als=al**2
  !
  ! calculate dv/dr for darwin correction
  !
  CALL my_lschps_derv(mmax, al, grid%r, v, dv )
  !
  !     starting of loop on energy for bound state
  !
  DO n_it = 1, maxter
     !
     ! coefficient array for u in differential eq.
     DO i=1,mmax
        cf(i)=als*(sls + (v(i)-e)*grid%r(i)**2)
     ENDDO
     !
     ! find classical turning point for matching
     !
     IF(mode == 1 .or. mode == 2) THEN
        DO i=mmax,2,-1
           IF(cf(i-1) <= 0.0_dp .and. cf(i) > 0.0_dp) THEN
              mch=i
              GOTO 40
           ENDIF
        ENDDO
        !PRINT '('' warning: wfc '',2i2,'' no turning point'')', n, l
        e=0.0_dp
        DO i=1,mmax
           u (i)=0.0_dp
        ENDDO
        nstop=1
        GOTO 999
     ELSE
        mch=nin
     ENDIF
40   CONTINUE

     !  relativistic coefficient arrays for u (fr) and up (frp).
     DO i=1,mmax
        fr(i)=als*(grid%r(i)**2)*0.25_dp*(-fss*(v(i)-e)**2 + &
             fss*dv(i)/ (grid%r(i)*(1.0_dp+0.25_dp*fss*(e-v(i)))))
        frp(i)=-al*grid%r(i)*0.25_dp*fss*dv(i)/(1.0_dp+0.25_dp*fss*(e-v(i)))
     ENDDO
     !
     ! start wavefunction with series
     !
     DO i=1,4
        u(i)=grid%r(i)**gamma
        up(i)=al*gamma*grid%r(i)**gamma
        upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
     ENDDO
     !
     ! outward integration using predictor once, corrector
     ! twice
     node=0
     !
     DO i=4,mch-1
        u(i+1) = u(i) + my_lschps_aeo(up,i)
        up(i+1) = up(i) + my_lschps_aeo(upp,i)
        DO it=1,2
           upp(i+1) = (al + frp(i+1))*up(i+1)+(cf(i+1)+fr(i+1))*u(i+1)
           up(i+1) = up(i) + my_lschps_aio(upp,i)
           u(i+1) = u(i) + my_lschps_aio(up,i)
        ENDDO
        IF(u(i+1)*u(i) <= 0.0_dp) node=node+1
     ENDDO
     !
     uout=u(mch)
     upout=up(mch)
     !
     IF(node-n+l+1 == 0 .or. mode == 3 .or. mode == 5) THEN
        !
        IF(mode == 1 .or. mode == 2) THEN
           !
           ! start inward integration at 10*classical turning
           ! point with simple exponential
           nin=mch+2.3_dp/al
           IF(nin+4 > mmax) nin=mmax-4
           xkap=sqrt(sls/grid%r(nin)**2 + 2.0_dp*(v(nin)-e))
           !
           DO i=nin,nin+4
              u(i)=exp(-xkap*(grid%r(i)-grid%r(nin)))
              up(i)=-grid%r(i)*al*xkap*u(i)
              upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
           ENDDO
           !
           ! integrate inward
           !
           DO i=nin,mch+1,-1
              u(i-1) = u(i) + my_lschps_aei(up,i)
              up(i-1) = up(i) + my_lschps_aei(upp,i)
              DO it = 1,2
                 upp(i-1) = (al + frp(i-1))*up(i-1) + (cf(i-1)+fr(i-1))*u(i-1)
                 up(i-1) = up(i) + my_lschps_aii(upp,i)
                 u(i-1) = u(i) + my_lschps_aii(up,i)
              ENDDO
           ENDDO
           !
           ! scale outside wf for continuity
           sc=uout/u(mch)
           !
           DO i=mch,nin
              up(i)=sc*up(i)
              u (i)=sc*u (i)
           ENDDO
           !
           upin=up(mch)
           !
        ELSE
           !
           upin=uld*uout
           !
        ENDIF
        !
        ! perform normalization sum
        !
        ro=grid%r(1)*exp(-0.5_dp*grid%dx)
        sn=ro**(2.0_dp*gamma+1.0_dp)/(2.0_dp*gamma+1.0_dp)
        !
        DO i=1,nin-3
           sn=sn+al*grid%r(i)*u(i)**2
        ENDDO
        !
        sn=sn + al*(23.0_dp*grid%r(nin-2)*u(nin-2)**2 &
             + 28.0_dp*grid%r(nin-1)*u(nin-1)**2 &
             +  9.0_dp*grid%r(nin  )*u(nin  )**2)/24.0_dp
        !
        ! normalize u
        cn=1.0_dp/sqrt(sn)
        uout=cn*uout
        upout=cn*upout
        upin=cn*upin
        !
        DO i=1,nin
           up(i)=cn*up(i)
           u(i)=cn*u(i)
        ENDDO
        DO i=nin+1,mmax
           u(i)=0.0_dp
        ENDDO
        !
        ! exit for fixed-energy calculation
        !
        IF(mode == 3 .or. mode == 5) GOTO 999

        ! perturbation theory for energy shift
        de=uout*(upout-upin)/(al*grid%r(mch))
        !
        ! convergence test and possible exit
        !
        IF ( abs(de) < max(abs(e),0.2_dp)*eps) GOTO 999
        !
        IF(de > 0.0_dp) THEN
           emin=e
        ELSE
           emax=e
        ENDIF
        e=e+de
        IF(e > emax .or. e < emin) e=0.5_dp*(emax+emin)
        !
     ELSEIF(node-n+l+1 < 0) THEN
        ! too few nodes
        emin=e
        e=0.5_dp*(emin+emax)

     ELSE
        ! too many nodes
        emax=e
        e=0.5_dp*(emin+emax)
     ENDIF
  ENDDO

  !PRINT '('' warning: wfc '',2i2,'' not converged'')', n, l
  u=0.0_dp
  nstop=1
  !
  ! deallocate arrays and exit
  !
999 CONTINUE
  DEALLOCATE(frp)
  DEALLOCATE(fr)
  DEALLOCATE(dv)
  DEALLOCATE(cf)
  DEALLOCATE(upp)
  DEALLOCATE(up)
  RETURN

END SUBROUTINE my_lschps


!---------------------------------------------------------------
FUNCTION my_lschps_d2u(al,cf,fr,frp,u,up) result(res)
!---------------------------------------------------------------
  ! second derivative from radial KS equation
  USE kinds, ONLY : DP
  IMPLICIT NONE
  real(dp):: res
  real(dp), INTENT(in):: al,cf,fr,frp,u,up
  !
  res = (al+frp)*up + (cf+fr)*u
  RETURN
END FUNCTION


FUNCTION my_lschps_aei(y,j) result(res)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER j
  real(DP):: y(j+3), res
  !
  res = -(4.16666666667e-2_dp)*(55.0_dp*y(j)-59.0_dp*y(j+1) + 37.0_dp*y(j+2)-9.0_dp*y(j+3))
  RETURN
END FUNCTION


!
! adams extrapolation and interpolation formulas for
! outward and inward integration, abramowitz and stegun, p. 896
FUNCTION my_lschps_aeo(y,j) result(res)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER:: j
  real(DP):: y(j), res
  !
  res = (4.16666666667e-2_dp)*(55.0_dp*y(j)-59.0_dp*y(j-1) + 37.0_dp*y(j-2)-9.0_dp*y(j-3))
  RETURN
END FUNCTION


FUNCTION my_lschps_aii(y,j) result(res)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER:: j
  real(DP) :: y(j+2), res
  !
  res = -(4.16666666667e-2_dp)*(9.0_dp*y(j-1)+19.0_dp*y(j) - 5.0_dp*y(j+1)+y(j+2))
  RETURN
END FUNCTION




FUNCTION my_lschps_aio(y,j) result(res)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: j
  real(DP):: y(j+1), res
  !
  res = (4.16666666667e-2_dp)*(9.0_dp*y(j+1)+19.0_dp*y(j) -5.0_dp*y(j-1)+y(j-2))
  RETURN
END FUNCTION


SUBROUTINE my_lschps_derV(mmax,al,r,v,dv)
  ! dv = dv/dr
  USE kinds, ONLY : dp
  IMPLICIT NONE
  INTEGER, INTENT(in)  :: mmax
  REAL(dp), INTENT(in) :: al, r(mmax), v(mmax)
  REAL(dp), INTENT(out):: dv(mmax)
  !
  INTEGER :: i
  !
  dv(1)=(-50.0_dp*v(1)+96.0_dp*v(2)-72.0_dp*v(3)+32.0_dp*v(4) &
       -6.0_dp*v(5))/(24.0_dp*al*r(1))
  dv(2)=(-6.0_dp*v(1)-20.0_dp*v(2)+36.0_dp*v(3)-12.0_dp*v(4) &
       +2.0_dp*v(5))/(24.0_dp*al*r(2))
  !
  DO i=3,mmax-2
     dv(i)=(2.0_dp*v(i-2)-16.0_dp*v(i-1)+16.0_dp*v(i+1) &
          -2.0_dp*v(i+2))/(24.0_dp*al*r(i))
  ENDDO
  !
  dv(mmax-1)=( 3.0_dp*v(mmax)+10.0_dp*v(mmax-1)-18.0_dp*v(mmax-2)+ &
       6.0_dp*v(mmax-3)-v(mmax-4))/(12.0_dp*al*r(mmax-1))
  dv(mmax)=( 25.0_dp*v(mmax)-48.0_dp*v(mmax-1)+36.0_dp*v(mmax-2)-&
       16.0_dp*v(mmax-3)+3.0_dp*v(mmax-4))/(12.0_dp*al*r(mmax))
  !
  RETURN
  !
END SUBROUTINE



