!-----------------------------------------------------------------------
SUBROUTINE my_stress_knl( sigmanlc, sigmakin )
!-----------------------------------------------------------------------
  !! Computes the kinetic + nonlocal contribuition to the stress
  !
  USE kinds,                ONLY: DP
  USE constants,            ONLY: pi, e2
  USE cell_base,            ONLY: omega, alat, at, bg, tpiba
  USE gvect,                ONLY: g
  USE gvecw,                ONLY: qcutz, ecfixed, q2sigma
  USE klist,                ONLY: nks, xk, ngk, igk_k
  USE io_files,             ONLY: iunwfc, nwordwfc
  USE buffers,              ONLY: get_buffer
  USE symme,                ONLY: symmatrix
  USE wvfct,                ONLY: npwx, nbnd, wg
  USE control_flags,        ONLY: gamma_only
  USE noncollin_module,     ONLY: noncolin, npol
  USE wavefunctions,        ONLY: evc
  USE mp_pools,             ONLY: inter_pool_comm
  USE mp_bands,             ONLY: intra_bgrp_comm
  USE mp,                   ONLY: mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP) :: sigmanlc (3, 3)
  !! non-local contribution to stress
  REAL(DP) :: sigmakin (3, 3)
  !! kinetic contribution to stress
  !
  ! local variables
  !
  REAL(DP), ALLOCATABLE :: gk (:,:), kfac (:)
  REAL(DP) :: twobysqrtpi, gk2, arg
  INTEGER  :: npw, ik, l, m, i, ibnd, is
  
  write(*,*)
  write(*,*) '*** ENTER my_stress_knl'
  write(*,*)

  !
  ALLOCATE( gk(npwx,3) )
  ALLOCATE( kfac(npwx) )
  !
  sigmanlc(:,:) = 0.d0
  sigmakin(:,:) = 0.d0
  twobysqrtpi   = 2.d0/SQRT(pi)
  !
  kfac(:) = 1.d0
  !
  DO ik = 1, nks
    write(*,*)
    !write(*,*) 'ik = ', ik
    IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )
    npw = ngk(ik)
    DO i = 1, npw
      gk(i,1) = ( xk(1,ik) + g(1,igk_k(i,ik)) ) * tpiba
      gk(i,2) = ( xk(2,ik) + g(2,igk_k(i,ik)) ) * tpiba
      gk(i,3) = ( xk(3,ik) + g(3,igk_k(i,ik)) ) * tpiba
      IF( qcutz > 0.d0 ) THEN
        gk2 = gk(i,1)**2 + gk(i,2)**2 + gk(i,3)**2
        arg = ( (gk2-ecfixed)/q2sigma )**2
        kfac(i) = 1.d0 + qcutz / q2sigma * twobysqrtpi * EXP(-arg)
      ENDIF
    ENDDO
    !write(*,'(1x,A,I5,F18.10)') 'ik, sum(gk) = ', ik, sum(gk(1:npw,1:3))
    !write
    !
    ! kinetic contribution
    !
    DO l = 1, 3
      DO m = 1, l
        DO ibnd = 1, nbnd
          DO i = 1, npw
            IF (noncolin) THEN
              sigmakin(l,m) = sigmakin(l,m) + wg(ibnd,ik) * &
               gk(i,l) * gk(i, m) * kfac(i) * &
               ( DBLE (CONJG(evc(   i  ,ibnd))*evc(   i  ,ibnd)) + &
                 DBLE (CONJG(evc(i+npwx,ibnd))*evc(i+npwx,ibnd)))
            ELSE
              sigmakin(l,m) = sigmakin(l,m) + wg(ibnd,ik) * &
                   gk(i,l) * gk(i, m) * kfac(i) * &
                  DBLE( CONJG( evc(i,ibnd) ) * evc(i, ibnd) )
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      !
    ENDDO
    !write(*,*) 'wg ik = ', wg(:,ik)
    !
    ! contribution from the  nonlocal part
    !
    CALL my_stress_us( ik, gk, sigmanlc )
    !
  ENDDO
  !write(*,*) 'kfac = ', kfac
  !
  DEALLOCATE( kfac )
  DEALLOCATE(  gk  )
  !
  ! the kinetic term must be summed over PW's and over k-points
  !
  CALL mp_sum( sigmakin, intra_bgrp_comm )
  CALL mp_sum( sigmakin, inter_pool_comm )
  !
  ! the nonlocal term is summed here only over k-points, because we add
  ! to it the US term from augmentation charge derivatives
  !
  CALL mp_sum( sigmanlc, inter_pool_comm )
  !
  ! add US term from augmentation charge derivatives, sum result over PW's
  !
  CALL my_addusstress( sigmanlc )
  CALL mp_sum( sigmanlc, intra_bgrp_comm )
  !
  DO l = 1, 3
    DO m = 1, l-1
      sigmanlc(m,l) = sigmanlc(l,m)
      sigmakin(m,l) = sigmakin(l,m)
    ENDDO
  ENDDO
  !
  IF ( gamma_only ) THEN
    sigmakin(:,:) = 2.d0 * e2 / omega * sigmakin(:,:)
  ELSE
    sigmakin(:,:) = e2 / omega * sigmakin(:,:)
  ENDIF
  !
  write(*,*)
  write(*,*) 'Kinetic energy stress, not symmetrized (Ry/bohr**3):'
  write(*,*)
  do l = 1,3
    write(*,'(1x,3F18.10)') sigmakin(l,1), sigmakin(l,2), sigmakin(l,3)
  enddo
  write(*,*)
  !
  sigmanlc(:,:) = -1.d0 / omega * sigmanlc(:,:)
  !
  write(*,*)
  write(*,*) 'PsPonNL (NCPP, USPP, or/and PAW) energy stress, not symmetrized (Ry/bohr**3):'
  write(*,*)
  do l = 1,3
    write(*,'(1x,3F18.10)') sigmanlc(l,1), sigmanlc(l,2), sigmanlc(l,3)
  enddo
  write(*,*)

  !
  ! symmetrize stress
  !
  CALL symmatrix( sigmakin )
  CALL symmatrix( sigmanlc )
  !
  RETURN
  !
END SUBROUTINE my_stress_knl

