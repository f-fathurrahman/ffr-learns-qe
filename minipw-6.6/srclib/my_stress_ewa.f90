!-----------------------------------------------------------------------
SUBROUTINE my_stress_ewa( alat, nat, ntyp, ityp, zv, at, bg, tau, &
                      omega, g, gg, ngm, gstart, gamma_only,  &
                      gcutm, sigmaewa )
  !---------------------------------------------------------------------
  !! Ewald contribution. Both real- and reciprocal-space terms are 
  !! present.
  !
  USE kinds
  USE constants,    ONLY : tpi, e2, eps6
  USE Coul_cut_2D,  ONLY : do_cutoff_2D, cutoff_stres_sigmaewa
  !
  IMPLICIT NONE
  !
  INTEGER :: nat !! input: number of atoms in the unit cell
  INTEGER :: ntyp !! input: number of different types of atoms
  INTEGER :: ityp(nat) !! input: the type of each atom
  INTEGER :: ngm !! input: number of plane waves for G sum
  INTEGER :: gstart !! input: first nonzero g vector
  LOGICAL, INTENT(IN) :: gamma_only !! gamma point only
  REAL(DP), INTENT(IN) :: tau(3,nat) !! input: the positions of the atoms in the cell
  REAL(DP), INTENT(IN) :: g(3,ngm) !! input: the coordinates of G vectors
  REAL(DP), INTENT(IN) :: gg(ngm) !! input: the square moduli of G vectors
  REAL(DP), INTENT(IN) :: zv(ntyp) !! input: the charge of each type of atoms
  REAL(DP), INTENT(IN) :: at(3,3) !! input: the direct lattice vectors
  REAL(DP), INTENT(IN) :: bg(3,3) !! input: the reciprocal lattice vectors
  REAL(DP), INTENT(IN) :: omega !! input: the volume of the unit cell
  REAL(DP), INTENT(IN) :: alat !! input: measure of length
  REAL(DP), INTENT(IN) :: gcutm !! input: cut-off of g vectors
  REAL(DP), INTENT(OUT) :: sigmaewa(3,3) ! output: the ewald stress
  !
  ! ... local variables
  !
  INTEGER, PARAMETER :: mxr = 50
  ! the maximum number of R vectors included in r sum
  integer :: ng,  nr, na, nb, l, m, nrm
  ! counter over reciprocal G vectors
  ! counter over direct vectors
  ! counter on atoms
  ! counter on atoms
  ! counter on atoms
  ! number of R vectors included in r sum
  !
  REAL(DP) :: charge, arg, tpiba2, dtau(3), alpha, r(3,mxr),     &
              r2(mxr), rmax, rr, upperbound, fact, fac, g2, g2a, &
              sdewald, sewald
  ! total ionic charge in the cell
  ! the argument of the phase
  ! length in reciprocal space
  ! the difference tau_s - tau_s'
  ! alpha term in ewald sum
  ! input of the rgen routine ( not used here )
  ! the square modulus of R_j-tau_s-tau_s'
  ! the maximum radius to consider real space sum
  ! buffer variable
  ! used to optimize alpha
  ! auxiliary variables
  ! diagonal term
  ! nondiagonal term
  COMPLEX(DP) :: rhostar
  REAL(DP), EXTERNAL :: qe_erfc ! the erfc function

  write(*,*)
  write(*,*) '*** ENTER my_stres_ewa'
  write(*,*)

  tpiba2 = (tpi / alat)**2
  sigmaewa(:,:) = 0.d0
  charge = 0.d0
  !
  DO na = 1, nat
     charge = charge + zv(ityp(na))
  ENDDO
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error ON THE ENERGY
  !
  alpha = 2.9d0
12 alpha = alpha - 0.1d0
  !
  IF (alpha==0.0) CALL errore( 'stres_ew', 'optimal alpha not found', 1 )
  upperbound = e2 * charge**2 * SQRT(2 * alpha / tpi) * &
               qe_erfc ( SQRT(tpiba2 * gcutm / 4.0d0 / alpha) )
  !write(*,*) 'alpha = ', alpha, 'upperbound = ', upperbound
  !
  IF (upperbound > 1d-7) GOTO 12

  write(*,*) 'Using alpha = ', alpha

  !
  ! G-space sum here
  !
  ! Determine if this processor contains G=0 and set the constant term
  !
  IF (gstart == 2) THEN
    sdewald = tpi * e2 / 4.d0 / alpha * (charge / omega)**2
  ELSE
    sdewald = 0.d0
  ENDIF
  write(*,*) 'e2 = ', e2
  write(*,*) 'sdewald before subtracted  (in Ha/bohr**2) = ', 0.5d0*sdewald
  !
  ! sdewald is the diagonal term
  IF( gamma_only ) THEN
    fact = 2.d0
  ELSE
    fact = 1.d0
  ENDIF
  !
  IF( do_cutoff_2D ) THEN 
    CALL cutoff_stres_sigmaewa( alpha, sdewald, sigmaewa )
  ELSE
    DO ng = gstart, ngm
      g2 = gg(ng) * tpiba2
      g2a = g2 / 4.d0 / alpha
      rhostar = (0.d0, 0.d0)
      DO na = 1, nat
        arg = (g(1,ng) * tau(1,na) + g(2,ng) * tau(2,na) + &
               g(3,ng) * tau(3,na) ) * tpi
        rhostar = rhostar + zv(ityp(na)) * CMPLX(COS(arg), SIN(arg), KIND=DP)
      ENDDO
      rhostar = rhostar / omega
      sewald = fact * tpi * e2 * EXP(-g2a) / g2 * ABS(rhostar)**2
      sdewald = sdewald - sewald
      DO l = 1, 3
         DO m = 1, l
            sigmaewa(l,m) = sigmaewa(l,m) + sewald * tpiba2 * 2.d0 * &
                            g(l,ng) * g(m,ng) / g2 * (g2a + 1)
         ENDDO
      ENDDO
      !
    ENDDO
  ENDIF
  write(*,*) 'sdewald after subtracted (in Ha/bohr**2) = ', 0.5d0*sdewald

  !
  DO l = 1, 3
    sigmaewa(l,l) = sigmaewa(l,l) + sdewald
  ENDDO
  !
  ! R-space sum here (see ewald.f90 for details on parallelization)
  !
  rmax = 4.0d0 / SQRT(alpha) / alat
  write(*,*) 'rmax*alat = ', rmax*alat
  !
  ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
  !
  DO na = 1, nat
    DO nb = 1, nat
      dtau(:) = tau(:,na) - tau(:,nb)
      !
      !     generates nearest-neighbors shells r(i)=R(i)-dtau(i)
      !
      CALL rgen( dtau, rmax, mxr, at, bg, r, r2, nrm )
      !
      DO nr = 1, nrm
        rr = SQRT(r2 (nr) ) * alat
        fac = - e2 / 2.0_DP/ omega * alat**2 * zv(ityp(na)) * &
              zv(ityp(nb)) / rr**3 * (qe_erfc(SQRT(alpha) * rr) + &
              rr * SQRT(8.0_DP * alpha / tpi) * EXP( - alpha * rr**2) )
        DO l = 1, 3
          DO m = 1, l
            sigmaewa(l,m) = sigmaewa(l,m) + fac * r(l,nr) * r(m,nr)
          ENDDO
        ENDDO
      ENDDO
      !
    ENDDO
  ENDDO
  !
  DO l = 1, 3
    DO m = 1, l - 1
      sigmaewa(m,l) = sigmaewa(l,m)
    ENDDO
  ENDDO
  !
  DO l = 1, 3
    DO m = 1, 3
      sigmaewa(l,m) = -sigmaewa(l,m)
    ENDDO
  ENDDO

  write(*,*)
  write(*,*) 'Ewald stress, not symmetrized (Ry/bohr**3):'
  write(*,*)
  do l = 1,3
    write(*,'(1x,3F18.10)') sigmaewa(l,1), sigmaewa(l,2), sigmaewa(l,3)
  enddo

  !
  write(*,*)
  write(*,*) '*** EXIT my_stres_ewa'
  write(*,*)
  !
  RETURN
  !
END SUBROUTINE my_stress_ewa

