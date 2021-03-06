SUBROUTINE my_force_ew( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, &
                     g, gg, ngm, gstart, gamma_only, gcutm, strf, forceion )
  !-----------------------------------------------------------------------
  !! This routine computes the Ewald contribution to the forces,
  !! both the real- and reciprocal-space terms are present.
  !
  USE kinds
  USE constants,    ONLY : tpi, e2
  USE mp_bands,     ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  USE mp,           ONLY : mp_sum
  USE Coul_cut_2D,  ONLY : do_cutoff_2D, cutoff_force_ew 
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nat
  !! the number of atoms
  INTEGER, INTENT(IN) :: ntyp
  !! the number of types of atom
  INTEGER, INTENT(IN) :: ngm
  !! the number of G vectors
  INTEGER, INTENT(IN) :: ityp(nat)
  !! the type of each atom
  INTEGER, INTENT(IN) :: gstart
  !! first non-zero G vector
  LOGICAL, INTENT(IN) :: gamma_only
  !! gamma only switch
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! the coordinates of the atoms
  REAL(DP), INTENT(IN) :: g(3,ngm)
  !! the G vectors
  REAL(DP), INTENT(IN) :: gg(ngm)
  !! the moduli of G vectors
  REAL(DP), INTENT(IN) :: zv(ntyp)
  !! the charge of the atoms
  REAL(DP), INTENT(IN) :: at(3,3)
  !! the direct lattice vectors
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! the reciprocal lattice vectors
  REAL(DP), INTENT(IN) :: omega
  !! the volume of the unit cell
  REAL(DP), INTENT(IN) :: gcutm
  !! cut-off of g vectors
  REAL(DP), INTENT(IN) :: alat
  !! the edge of the cell
  COMPLEX(DP), INTENT(IN) :: strf(ngm,ntyp)
  !! the structure factor on the potential
  REAL(DP), INTENT(OUT) :: forceion(3,nat)
  !! the ewald part of the forces
  !
  ! ... local variables
  !
  INTEGER, PARAMETER :: mxr=50
  ! the maximum number of R vectors
  !
  INTEGER :: ig, n, na, nb, nt, nrm, ipol
  ! counter on G vectos
  ! counter on r vectors
  ! counter on atoms
  ! counter on atoms
  ! counter on atomic types
  ! the number of R vectors for real space su
  ! counter on polarization
  INTEGER :: na_s, na_e, mykey
  !
  REAL(DP) :: sumnb, arg, tpiba2, alpha, dtau(3), r(3,mxr), &
              r2(mxr), rmax, rr, charge, upperbound, fact
  ! auxiliary variable for speed
  ! the argument of the exponential
  ! 2 pi /alat
  ! the alpha parameter
  ! the difference of two tau
  ! the position of the atoms in the shell
  ! the square of r
  ! the maximum r
  ! the modulus of the r vectors
  ! the total charge
  ! used to determine alpha
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  real(8), allocatable :: forceion_R(:,:)
  ! auxiliary space
  REAL(DP), EXTERNAL :: qe_erfc

  !
  forceion(:,:) = 0.0_DP
  tpiba2 = (tpi / alat)**2
  charge = 0.0_DP
  DO na = 1, nat
     charge = charge + zv(ityp(na))
  ENDDO

  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error ON THE ENERGY
  alpha = 1.1_DP
10 alpha = alpha - 0.1_DP
  IF (alpha==0.d0) CALL errore( 'force_ew', 'optimal alpha not found', 1 )
  upperbound = e2 * charge**2 * SQRT(2._DP * alpha / tpi) * &
               qe_erfc( SQRT(tpiba2 * gcutm / 4.d0 / alpha) )
  IF ( upperbound > 1.0d-6 ) GOTO 10

  write(*,*) 'upperbound = ', upperbound
  write(*,*) 'alpha = ', alpha
  write(*,*) 'tpiba2*gcutm = ', tpiba2*gcutm


  !
  ! G-space sum here
  !
  ALLOCATE( aux(ngm) )
  aux(:) = (0.0_DP, 0.0_DP)
  !
  DO nt = 1, ntyp
     DO ig = gstart, ngm
        aux(ig) = aux(ig) + zv(nt) * CONJG( strf(ig,nt) )
     ENDDO
  ENDDO

  !write(*,*) 'ngm = ', ngm
  !write(*,*) 'sum(strf) = ', sum(strf)
  !write(*,*) 'sum(aux)  = ', sum(aux)

  !write(*,*) 'some gg and aux'
  !do ig = 2,10
  !  write(*,'(1x,I8,3F18.10)') ig, gg(ig)*tpiba2, aux(ig)
  !enddo

  !
  IF (do_cutoff_2D) THEN 
     CALL cutoff_force_ew( aux, alpha )
  ELSE
     DO ig = gstart, ngm
        aux(ig) = aux(ig) * EXP(-gg(ig) * tpiba2 / alpha / 4.d0) / &
                  (gg(ig) * tpiba2)
     ENDDO
  ENDIF

  !write(*,*) 'some gg and aux after modified'
  !do ig = 2,10
  !  write(*,'(1x,I8,3F18.10)') ig, gg(ig)*tpiba2, aux(ig)
  !enddo

  IF (gamma_only) THEN
     fact = 4.0_DP
  ELSE
     fact = 2.0_DP
  ENDIF
  !
  DO na = 1, nat
     DO ig = gstart, ngm
        arg = tpi * (g(1,ig) * tau(1,na) + g(2,ig) * tau(2,na) + g(3,ig) * tau(3,na))
        sumnb = COS(arg)*AIMAG(aux(ig)) - SIN(arg)*DBLE(aux(ig))
        forceion(1,na) = forceion(1,na) + g(1,ig) * sumnb
        forceion(2,na) = forceion(2,na) + g(2,ig) * sumnb
        forceion(3,na) = forceion(3,na) + g(3,ig) * sumnb
     ENDDO
     DO ipol = 1, 3
        forceion(ipol,na) = - zv (ityp(na) ) * fact * e2 * tpi**2 / &
                            omega / alat * forceion(ipol,na)
     ENDDO
  ENDDO

  write(*,*) 'G-space contribution: '
  DO na = 1,nat
    write(*,'(1x,A,3F18.10)') 'forceion: ', forceion(:,na)
  ENDDO

  write(*,*)
  write(*,*) 'gstart = ', gstart
  write(*,*) 'e2 = ', e2
  write(*,*) 'gcutm = ', gcutm
  write(*,*) 'sqrt(tpiba2) = ', sqrt(tpiba2)
  write(*,*) 'tpiba2 = ', tpiba2
  write(*,*) 'alat = ', alat
  write(*,*) 'omega = ', omega
  write(*,*) 'fact = ', fact
  write(*,*) 'G2 max = ', gcutm*tpiba2
  write(*,*) 'alpha = ', alpha
  write(*,'(1x,A,ES18.10)') 'factor = ', -1.d0*fact * e2 * tpi**2 / omega / alat
  write(*,*)

  DEALLOCATE( aux )


  allocate(forceion_R(3,nat))
  forceion_R(:,:) = 0.d0

  !
  ! R-space sum here (see ewald.f90 for details on parallelization)
  !
  CALL block_distribute( nat, me_bgrp, nproc_bgrp, na_s, na_e, mykey )
  IF ( mykey > 0 ) GOTO 100
  rmax = 5.d0 / (SQRT(alpha) * alat)
  !
  ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
  !
  DO na = na_s, na_e
     DO nb = 1, nat
        IF (nb == na) GOTO 50
        dtau(:) = tau(:,na) - tau(:,nb)
        !
        ! generates nearest-neighbors shells r(i)=R(i)-dtau(i)
        !
        CALL rgen( dtau, rmax, mxr, at, bg, r, r2, nrm )
        write(*,*) 'dtau = ', dtau
        write(*,*) 'nrm = ', nrm
        DO n = 1, nrm
           rr = SQRT(r2(n)) * alat
           fact = zv(ityp(na)) * zv(ityp(nb)) * e2 / rr**2 * &
                  (qe_erfc(SQRT(alpha) * rr) / rr +          &
                  SQRT(8.0d0 * alpha / tpi) * EXP(- alpha * rr**2) ) * alat
           DO ipol = 1, 3
              forceion_R(ipol,na) = forceion_R(ipol,na) - fact * r(ipol,n)
           ENDDO
        ENDDO
50      CONTINUE
     ENDDO
  ENDDO
100 CONTINUE


  write(*,*) 'R-space contribution: '
  DO na = 1,nat
    write(*,'(1x,A,3F18.10)') 'forceion_R: ', forceion_R(:,na)
  ENDDO
  

  forceion(:,:) = forceion(:,:) + forceion_R(:,:)

  write(*,*) 'Plus R-space contribution: '
  DO na = 1,nat
    write(*,'(1x,A,3F18.10)') 'forceion: ', forceion(:,na)
  ENDDO
  !
  CALL mp_sum( forceion, intra_bgrp_comm )

  deallocate(forceion_R)


  RETURN
  !
END SUBROUTINE


!-----------------------------------------------------
PROGRAM main
!-----------------------------------------------------

  USE kinds,             ONLY : DP
  USE cell_base,         ONLY : at, bg, alat, omega  
  USE ions_base,         ONLY : nat, ntyp => nsp, ityp, tau, zv
  USE vlocal,            ONLY : strf
  USE gvect,             ONLY : ngm, gstart, g, gg, gcutm
  USE control_flags,     ONLY : gamma_only

  IMPLICIT NONE
  INTEGER :: ia
  REAL(DP), ALLOCATABLE :: forceion(:,:)

  CALL prepare_all()

  ALLOCATE(forceion(3,nat))

  CALL my_force_ew( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, &
                     g, gg, ngm, gstart, gamma_only, gcutm, strf, forceion )

  write(*,*) 'Result: '
  DO ia = 1,nat
    write(*,'(1x,A,3F18.10)') 'forceion: ', forceion(:,ia)
  ENDDO

  DEALLOCATE( forceion )

END PROGRAM