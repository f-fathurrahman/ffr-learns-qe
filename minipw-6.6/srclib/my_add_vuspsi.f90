!----------------------------------------------------------------------------
SUBROUTINE my_add_vuspsi( lda, n, m, hpsi )
!----------------------------------------------------------------------------
  !! This routine applies the Ultra-Soft Hamiltonian to a
  !! vector psi and puts the result in hpsi. 
  !! It requires the products of psi with all beta functions
  !! in array becp(nkb,m) (calculated by calbec).
  !
  USE kinds,           ONLY: DP
  USE ions_base,       ONLY: nat, ntyp => nsp, ityp
  USE lsda_mod,        ONLY: current_spin
  USE control_flags,   ONLY: gamma_only
  USE noncollin_module
  USE uspp,            ONLY: vkb, nkb, deeq, deeq_nc, indv_ijkb0
  USE uspp_param,      ONLY: nh, nhm
  USE becmod,          ONLY: bec_type, becp
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda*npol,m)
  !! V_US|psi> is added to hpsi
  !
  ! ... here the local variables
  !
  INTEGER :: jkb, ikb, ih, jh, na, nt, ibnd ! counters
  !
  IF( gamma_only ) THEN
    stop 'not yet supported in my_add_vuspsi 40'
  ELSEIF( noncolin) THEN
    stop 'not yet supported in my_add_vuspsi 42'
    !
  ELSE
    CALL add_vuspsi_k()
  ENDIF

  RETURN


CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE add_vuspsi_k()
!-----------------------------------------------------------------------
  !! See add_vuspsi_gamma for comments
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), ALLOCATABLE :: ps(:,:), deeaux(:,:)
  INTEGER :: ierr
  !
  IF( nkb == 0 ) RETURN

  ALLOCATE( ps(nkb,m), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' add_vuspsi_k ', ' cannot allocate ps ', ABS( ierr ) )

  DO nt = 1, ntyp
    
    IF( nh(nt) == 0 ) CYCLE
    
    ALLOCATE( deeaux(nh(nt),nh(nt)) )
    DO na = 1, nat
      IF( ityp(na) == nt ) THEN
        ! deeq is real: copy it into a complex variable to perform
        ! a zgemm - simple but sub-optimal solution
        deeaux(:,:) = CMPLX(deeq(1:nh(nt),1:nh(nt),na,current_spin), 0.0_dp, KIND=dp )
        !
        CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
                   deeaux, nh(nt), becp%k(indv_ijkb0(na)+1,1), nkb, &
                  (0.0_dp, 0.0_dp), ps(indv_ijkb0(na)+1,1), nkb )
        !
      ENDIF
      !
    ENDDO
    DEALLOCATE( deeaux )
    !
  ENDDO
  !
  CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ) , vkb, &
              lda, ps, nkb, ( 1.D0, 0.D0 ) , hpsi, lda )
  !
  DEALLOCATE( ps )
  !
  RETURN
  !
END SUBROUTINE




END SUBROUTINE
