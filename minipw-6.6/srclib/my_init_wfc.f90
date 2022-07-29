!
!----------------------------------------------------------------------------
SUBROUTINE my_init_wfc( ik )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes starting wavefunctions for k-point ik
  !
  USE kinds,                ONLY : DP
  USE bp,                   ONLY : lelfield
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   bec_type, becp
  USE constants,            ONLY : tpi
  USE basis,                ONLY : natomwfc, starting_wfc
  USE gvect,                ONLY : g, gstart
  USE klist,                ONLY : xk, ngk, igk_k
  USE wvfct,                ONLY : nbnd, npwx, et
  USE uspp,                 ONLY : nkb, okvan
  USE noncollin_module,     ONLY : npol
  USE wavefunctions, ONLY : evc
  USE random_numbers,       ONLY : randy
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm, &
                                   nbgrp, root_bgrp_id
  USE mp,                   ONLY : mp_bcast
  USE funct,                ONLY : dft_is_hybrid, stop_exx
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ik
  !
  INTEGER :: ibnd, ig, ipol, n_starting_wfc, n_starting_atomic_wfc
  LOGICAL :: lelfield_save
  !
  REAL(DP) :: rr, arg
  REAL(DP), ALLOCATABLE :: etatom(:) ! atomic eigenvalues
  !
  COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:,:) ! atomic wfcs for initialization


  write(*,*) 'my_init_wfc is called for ik = ', ik

  IF( starting_wfc(1:6) == 'atomic' ) THEN
    !
    n_starting_wfc = MAX( natomwfc, nbnd )
    n_starting_atomic_wfc = natomwfc
    !
  ELSE IF( starting_wfc == 'random' ) THEN
    !
    n_starting_wfc = nbnd
    n_starting_atomic_wfc = 0
    !
  ELSE
    !
    ! ...case 'file' should not be done here
    !
    CALL errore ( 'init_wfc', &
         'invalid value for startingwfc: ' // TRIM ( starting_wfc ) , 1 )
    !
  ENDIF
  !
  ALLOCATE( wfcatom( npwx, npol, n_starting_wfc ) )
  !
  IF ( starting_wfc(1:6) == 'atomic' ) THEN

    CALL my_atomic_wfc( ik, wfcatom )
    !
    IF( starting_wfc == 'atomic+random' .AND. &
      n_starting_wfc == n_starting_atomic_wfc ) THEN

      ! in this case, introduce a small randomization of wavefunctions
      ! to prevent possible "loss of states"
      DO ibnd = 1, n_starting_atomic_wfc
        DO ipol = 1, npol
          DO ig = 1, ngk(ik)
            rr  = randy()
            arg = tpi * randy()
            wfcatom(ig,ipol,ibnd) = wfcatom(ig,ipol,ibnd) * &
               ( 1.0_DP + 0.05_DP * CMPLX( rr*COS(arg), rr*SIN(arg), kind=DP) )
          ENDDO ! ig
        ENDDO ! ipol
      ENDDO ! ibnd
    ENDIF ! starting_wfc = atomic+random
  ENDIF ! atomic

  ! if not enough atomic wfc are available,
  ! fill missing wfcs with random numbers
  DO ibnd = n_starting_atomic_wfc + 1, n_starting_wfc
    DO ipol = 1, npol
      wfcatom(:,ipol,ibnd) = (0.0_dp, 0.0_dp)
      DO ig = 1, ngk(ik)
        rr  = randy()
        arg = tpi * randy()
        wfcatom(ig,ipol,ibnd) = &
             CMPLX( rr*COS( arg ), rr*SIN( arg ) ,kind=DP) / &
                    ( ( xk(1,ik) + g(1,igk_k(ig,ik)) )**2 + &
                      ( xk(2,ik) + g(2,igk_k(ig,ik)) )**2 + &
                      ( xk(3,ik) + g(3,igk_k(ig,ik)) )**2 + 1.0_DP )
      ENDDO
    ENDDO
  ENDDO
  
  ! when band parallelization is active, the first band group distributes
  ! the wfcs to the others making sure all bgrp have the same starting wfc
  ! FIXME: maybe this should be done once evc are computed, not here?
  !
  IF( nbgrp > 1 ) CALL mp_bcast( wfcatom, root_bgrp_id, inter_bgrp_comm )

  ! Diagonalize the Hamiltonian on the basis of atomic wfcs
  ALLOCATE( etatom( n_starting_wfc ) )

  ! Allocate space for <beta|psi>
  CALL allocate_bec_type( nkb, n_starting_wfc, becp, intra_bgrp_comm )

  ! the following trick is for electric fields with Berry's phase:
  ! by setting lelfield = .false. one prevents the calculation of
  ! electric enthalpy in the Hamiltonian (cannot be calculated
  ! at this stage: wavefunctions at previous step are missing)
  lelfield_save = lelfield
  lelfield = .FALSE.

  ! subspace diagonalization (calls Hpsi)
  IF ( dft_is_hybrid()  ) CALL stop_exx() 

  CALL rotate_wfc( npwx, ngk(ik), n_starting_wfc, gstart, nbnd, wfcatom, npol, okvan, evc, etatom )

  lelfield = lelfield_save   ! ffr: set lelfield to its original value

  ! copy the first nbnd eigenvalues
  ! eigenvectors are already copied inside routine rotate_wfc
  et(1:nbnd,ik) = etatom(1:nbnd)

  CALL deallocate_bec_type ( becp )
  DEALLOCATE( etatom )
  DEALLOCATE( wfcatom )

  RETURN

END SUBROUTINE my_init_wfc
