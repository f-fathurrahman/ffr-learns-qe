!----------------------------------------------------------------------------
SUBROUTINE my_stress_us( ik, gk, sigmanlc )
!----------------------------------------------------------------------------
  !! nonlocal (separable pseudopotential) contribution to the stress
  !! NOTICE: sum of partial results over procs is performed in calling routine
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE constants,            ONLY : eps8
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wvfct,                ONLY : npwx, nbnd, wg, et
  USE control_flags,        ONLY : gamma_only
  USE uspp_param,           ONLY : upf, lmaxkb, nh, nhm
  USE uspp,                 ONLY : nkb, vkb, deeq, deeq_nc
  USE wavefunctions,        ONLY : evc
  USE spin_orb,             ONLY : lspinorb
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp_pools,             ONLY : me_pool, root_pool
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   bec_type, becp, calbec
  USE mp,                   ONLY : mp_sum, mp_get_comm_null, mp_circular_shift_left 
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: ik
  !! k-point index
  REAL(DP), INTENT(IN)    :: gk(npwx,3)
  !! wave function components for fixed k-point
  REAL(DP), INTENT(INOUT) :: sigmanlc(3,3)
  !! stress tensor, non-local contribution
  !
  REAL(DP), ALLOCATABLE   :: qm1(:)
  REAL(DP)                :: q
  INTEGER                 :: npw, i
  !
  !
  IF ( nkb == 0 ) RETURN
  !
  IF ( lsda ) current_spin = isk(ik)
  npw = ngk(ik)
  IF ( nks > 1 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
  !
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm ) 
  CALL calbec( npw, vkb, evc, becp )
  !
  ALLOCATE( qm1( npwx ) )
  DO i = 1, npw
    q = SQRT( gk(i, 1)**2 + gk(i, 2)**2 + gk(i, 3)**2 )
    IF ( q > eps8 ) THEN
      qm1(i) = 1.D0 / q
    ELSE
      qm1(i) = 0.D0
    ENDIF
  ENDDO
  !
  IF( gamma_only ) THEN
    !
    stop 'my_stress_us_gamma is disabled'
    !
  ELSE
    !
    CALL my_stress_us_k()
    !
  ENDIF
  !
  DEALLOCATE( qm1 )
  CALL deallocate_bec_type( becp ) 
  !
  RETURN
  !
  CONTAINS

! ffr: noncolin statements are removed
!----------------------------------------------------------------------
SUBROUTINE my_stress_us_k()
!----------------------------------------------------------------------  
  !! nonlocal contribution to the stress - k-points version
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  INTEGER :: na, np, ibnd, ipol, jpol, l, i, ipw, &
             ikb, jkb, ih, jh, ijkb0, is, js, ijs
  REAL(DP) :: fac, xyz (3, 3), evps, ddot
  COMPLEX(DP), ALLOCATABLE :: work1(:), work2(:), dvkb(:,:)
  REAL(DP), ALLOCATABLE :: deff(:,:,:)
  ! dvkb contains the derivatives of the kb potential
  COMPLEX(DP) :: ps
  ! xyz are the three unit vectors in the x,y,z directions
  DATA xyz / 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /
  !
  !
  ALLOCATE( deff(nhm,nhm,nat) )
  Deff(:,:,:) = 0.d0
  !
  ALLOCATE( work1(npwx), work2(npwx) )
  !
  evps = 0.D0
  ! diagonal contribution
  !
  IF ( me_bgrp /= root_bgrp ) GO TO 100
  !
  ! the contribution is calculated only on one processor because
  ! partial results are later summed over all processors
  !
  DO ibnd = 1, nbnd
    fac = wg(ibnd,ik)
    IF( ABS(fac) < 1.d-9 ) CYCLE
    !
    CALL compute_deff(deff, et(ibnd,ik))
    !write(*,*)
    !write(*,*) 'ibnd, ebands (in Ha), fac = ', ibnd, 0.5d0*et(ibnd,ik), fac
    ijkb0 = 0
    DO np = 1, ntyp
      DO na = 1, nat
        IF( ityp(na) == np ) THEN
          DO ih = 1, nh(np)
            ikb = ijkb0 + ih
            evps = evps + fac * deff(ih,ih,na) * ABS(becp%k(ikb,ibnd))**2
            !write(*,'(1x,A,I4,I4,F18.10)') 'na, ih, Deff (in Ha) = ', na, ih, deff(ih,ih,na)*0.5d0
            !write(*,*) 'Current value of evps (in Ha) = ', evps*0.5d0
            !
            IF ( upf(np)%tvanp .OR. upf(np)%is_multiproj ) THEN
              !
              ! only in the US case there is a contribution 
              ! for jh<>ih
              ! we use here the symmetry in the interchange of 
              ! ih and jh
              !
              DO jh = (ih + 1), nh(np)
                jkb = ijkb0 + jh
                evps = evps + deff(ih,jh,na) * fac * 2.D0 * DBLE(CONJG(becp%k(ikb,ibnd)) * becp%k(jkb,ibnd) )
                !write(*,*) 'Nondiagonal: evps = ', evps*0.5d0, "Deff = ", Deff(ih,jh,na)*0.5d0
              ENDDO
            ENDIF
          ENDDO
          ijkb0 = ijkb0 + nh(np)
        ENDIF
      ENDDO ! na
    ENDDO ! np
  ENDDO ! ibnd
  !write(*,*) 'evps (in Ha) = ', evps*0.5d0
  !
  DO l = 1, 3
    sigmanlc(l,l) = sigmanlc(l,l) - evps
  ENDDO
  !
  !
100 CONTINUE
  !
  ! non diagonal contribution - derivative of the bessel function
  !
  ALLOCATE( dvkb( npwx, nkb ) )
  !
  CALL my_gen_us_dj( ik, dvkb )
  !write(*,'(1x,A,I4,2F18.10)') 'ik, sum(dvkb Bessel) = ', ik, sum(dvkb)
  !
  DO ibnd = 1, nbnd
    work2 = (0.D0,0.D0)
    CALL compute_deff(deff,et(ibnd,ik))
    !
    ijkb0 = 0
    DO np = 1, ntyp
      DO na = 1, nat
        IF ( ityp(na) == np ) THEN
          DO ih = 1, nh(np)
            ikb = ijkb0 + ih
            IF( .NOT. ( upf(np)%tvanp .OR. upf(np)%is_multiproj ) ) THEN
              ps = becp%k(ikb, ibnd) * deeq(ih,ih,na,current_spin)
            ELSE
              !
              ! in the US case there is a contribution 
              ! also for jh<>ih
              !
              ps = (0.D0,0.D0)
              DO jh = 1, nh(np)
                jkb = ijkb0 + jh
                ps = ps + becp%k(jkb,ibnd) * deff(ih,jh,na)
              ENDDO
            ENDIF
            !
            DO ipw = 1, npw
              work2(ipw) = ps * dvkb(ipw, ikb) + work2(ipw)
            ENDDO
            !
          ENDDO
          ijkb0 = ijkb0 + nh(np)
        ENDIF
      ENDDO
    ENDDO
    !
    DO ipol = 1, 3
      DO jpol = 1, ipol
        DO i = 1, npw
           work1(i) = evc(i,ibnd)*gk(i, ipol)*gk(i, jpol)*qm1(i)
        ENDDO
        sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                              2.D0 * wg(ibnd,ik) * &
                              ddot( 2*npw, work1, 1, work2, 1 )
      ENDDO
    ENDDO
  ENDDO
  !write(*,*)
  !write(*,*) 'sigmanlc after adding deriv Bessel contrib, not symmetrized (Ry/bohr**3):'
  !write(*,*)
  !do l = 1,3
  !  write(*,'(1x,3F18.10)') sigmanlc(l,1), sigmanlc(l,2), sigmanlc(l,3)
  !enddo
  !write(*,*)


  !
  ! non diagonal contribution - derivative of the spherical harmonics
  ! (no contribution from l=0)
  !
  IF ( lmaxkb == 0 ) GO TO 10
  !
  DO ipol = 1, 3
    !
    CALL my_gen_us_dy( ik, xyz(1,ipol), dvkb )
    !
    DO ibnd = 1, nbnd
      work2 = (0.D0,0.D0)
      CALL compute_deff( deff, et(ibnd,ik) )
      ijkb0 = 0
      DO np = 1, ntyp
        DO na = 1, nat
          IF ( ityp(na) == np ) THEN
            DO ih = 1, nh(np)
              ikb = ijkb0 + ih
              IF ( .NOT. ( upf(np)%tvanp .OR. upf(np)%is_multiproj ) ) THEN
                ps = becp%k(ikb,ibnd) * deeq(ih,ih,na,current_spin)
              ELSE
                !
                ! in the US case there is a contribution 
                ! also for jh<>ih
                !
                ps = (0.D0,0.D0)
                DO jh = 1, nh(np)
                  jkb = ijkb0 + jh
                  ps = ps + becp%k(jkb,ibnd) * deff(ih,jh,na)
                ENDDO
              ENDIF
              !
              DO ipw = 1, npw
                work2(ipw) = ps * dvkb(ipw, ikb) + work2(ipw)
              ENDDO
            ENDDO
            ijkb0 = ijkb0 + nh(np)
          ENDIF
        ENDDO
      ENDDO
      !
      DO jpol = 1, ipol
        DO i = 1, npw
           work1(i) = evc(i,ibnd) * gk(i, jpol)
        ENDDO
        sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - 2.D0 * wg(ibnd,ik) * & 
                           ddot( 2 * npw, work1, 1, work2, 1 )
      ENDDO
    ENDDO
  ENDDO
  !
10 CONTINUE
  !
  DEALLOCATE( work2 )
  DEALLOCATE( deff )
  DEALLOCATE( dvkb )
  DEALLOCATE( work1 )
  !
  RETURN
  !
END SUBROUTINE my_stress_us_k


END SUBROUTINE my_stress_us
