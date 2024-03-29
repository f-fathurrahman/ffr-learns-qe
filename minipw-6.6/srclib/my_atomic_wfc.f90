!-----------------------------------------------------------------------
SUBROUTINE my_atomic_wfc( ik, wfcatom )
!-----------------------------------------------------------------------
  !! This routine computes the superposition of atomic wavefunctions
  !! for k-point "ik" - output in "wfcatom".
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi, fpi, pi
  USE cell_base,        ONLY : omega, tpiba
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,            ONLY : natomwfc
  USE gvect,            ONLY : mill, eigts1, eigts2, eigts3, g
  USE klist,            ONLY : xk, igk_k, ngk
  USE wvfct,            ONLY : npwx
  USE us,               ONLY : tab_at, dq
  USE uspp_param,       ONLY : upf
  USE noncollin_module, ONLY : noncolin, npol, angle1, angle2
  USE spin_orb,         ONLY : lspinorb, rot_ylm, fcoef, lmaxx, domag, &
                               starting_spin_angle
  USE mp_bands,         ONLY : inter_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! k-point index
  COMPLEX(DP), INTENT(OUT) :: wfcatom( npwx, npol, natomwfc )
  !! Superposition of atomic wavefunctions
  !
  ! ... local variables
  !
  INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3, nwfcm, npw
  REAL(DP),    ALLOCATABLE :: qg(:), ylm (:,:), chiq(:,:,:), gk(:,:)
  COMPLEX(DP), ALLOCATABLE :: sk(:), aux(:)
  COMPLEX(DP) :: kphase, lphase
  REAL(DP)    :: arg, px, ux, vx, wx
  INTEGER     :: ig_start, ig_end

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  DO nt = 1, ntyp
    lmax_wfc = MAX( lmax_wfc, MAXVAL( upf(nt)%lchi(1:upf(nt)%nwfc) ) )
  ENDDO
  !
  nwfcm = MAXVAL( upf(1:ntyp)%nwfc )
  npw = ngk(ik)
  !
  ALLOCATE( ylm(npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,ntyp), &
             sk(npw), gk(3,npw), qg(npw) )
  !
  DO ig = 1, npw
    iig = igk_k (ig,ik)
    gk(1,ig) = xk(1,ik) + g(1,iig)
    gk(2,ig) = xk(2,ik) + g(2,iig)
    gk(3,ig) = xk(3,ik) + g(3,iig)
    qg(ig) = gk(1,ig)**2 +  gk(2,ig)**2 + gk(3,ig)**2
  ENDDO
  !
  !  ylm = spherical harmonics
  !
  CALL ylmr2( (lmax_wfc+1)**2, npw, gk, qg, ylm )

  ! from now to the end of the routine the ig loops are distributed across bgrp
  CALL divide( inter_bgrp_comm,npw,ig_start,ig_end )
  !
  ! set now q=|k+G| in atomic units
  !
  DO ig = ig_start, ig_end
     qg(ig) = SQRT( qg(ig) )*tpiba
  END DO
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  DO nt = 1, ntyp
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc (nb) >= 0.d0 ) THEN
           DO ig = ig_start, ig_end
              px = qg(ig) / dq - INT(qg(ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT( qg (ig) / dq ) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq (ig, nb, nt) = &
                     tab_at(i0, nb, nt) * ux * vx * wx / 6.d0 + &
                     tab_at(i1, nb, nt) * px * vx * wx / 2.d0 - &
                     tab_at(i2, nb, nt) * px * ux * wx / 2.d0 + &
                     tab_at(i3, nb, nt) * px * ux * vx / 6.d0
           END DO
        END IF
     END DO
  END DO

  DEALLOCATE( qg, gk )
  ALLOCATE( aux(npw) )
  !
  wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
  !
  DO na = 1, nat
    arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
    kphase = CMPLX( COS(arg), - SIN(arg) ,KIND=DP)
    !
    ! sk is the structure factor
    !
    DO ig = ig_start, ig_end
      iig = igk_k (ig,ik)
      sk (ig) = kphase * eigts1(mill(1,iig), na) * &
                         eigts2(mill(2,iig), na) * &
                         eigts3(mill(3,iig), na)
    ENDDO
    !
    nt = ityp (na)
    DO nb = 1, upf(nt)%nwfc
      IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
        l = upf(nt)%lchi(nb)
        lphase = (0.d0,1.d0)**l
        !
        !  the factor i^l MUST BE PRESENT in order to produce
        !  wavefunctions for k=0 that are real in real space
        !
        IF( noncolin ) THEN
          stop 'noncolin is not supported in my_atomic_wfc'
        ELSE
          !
          CALL my_atomic_wfc___( )
          !
        ENDIF
      ENDIF
    ENDDO ! loop over nwfc
  ENDDO ! loop over natoms

  IF( n_starting_wfc /= natomwfc) call errore ('my_atomic_wfc internal error: some wfcs were lost ', 1 )

  DEALLOCATE( aux, sk, chiq, ylm )

  ! collect results across bgrp
  CALL mp_sum( wfcatom, inter_bgrp_comm )

  RETURN

CONTAINS

! Inner subroutine
! LSDA or nonmagnetic case
!--------------------------
SUBROUTINE my_atomic_wfc___( )
!--------------------------
  DO m = 1, 2 * l + 1
    lm = l**2 + m
    n_starting_wfc = n_starting_wfc + 1
    IF ( n_starting_wfc > natomwfc) CALL errore('atomic_wfc___', 'internal error: too many wfcs', 1)
    DO ig = ig_start, ig_end
      wfcatom(ig, 1, n_starting_wfc) = lphase * sk(ig) * ylm(ig,lm) * chiq(ig,nb,nt)
    ENDDO
  ENDDO
END SUBROUTINE my_atomic_wfc___


END SUBROUTINE my_atomic_wfc
