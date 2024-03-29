
!-----------------------------------------------------------------------
SUBROUTINE my_compute_deff( deff, et )
!-----------------------------------------------------------------------
  !!  This routine computes the effective value of the D-eS coefficients
  !!  which appear often in many expressions in the US or PAW case. 
  !!  This routine is for the collinear case.
  !
  USE kinds,       ONLY: DP
  USE ions_base,   ONLY: nsp, nat, ityp
  USE uspp,        ONLY: deeq, qq_at, okvan
  USE uspp_param,  ONLY: nhm
  USE lsda_mod,    ONLY: current_spin
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: et
  !! The eigenvalues of the hamiltonian
  REAL(DP), INTENT(OUT) :: deff(nhm,nhm,nat)
  !! Effective values of the D-eS coefficients
  !
  ! ... local variables
  !
  INTEGER :: nt, na, is
  !
  deff(:,:,:) = deeq(:,:,:,current_spin)
  !
  IF(okvan) THEN
    !
    DO nt = 1, nsp
      DO na = 1, nat
        !
        IF( ityp(na) == nt ) THEN
          deff(:,:,na) = deff(:,:,na) - et*qq_at(:,:,na)
        ENDIF
        !
      ENDDO
    ENDDO
    !
    !write(*,*) 'sum qq_at = ', sum(qq_at)
    !write(*,*) 'sum Deeq = ', sum(Deeq)
    !write(*,*) 'sum Deff (in Ha) = ', 0.5d0*sum(Deff)
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE my_compute_deff
!
!
!---------------------------------------------------------------------------
SUBROUTINE my_compute_deff_nc( deff, et )
  !-------------------------------------------------------------------------
  !! This routine computes the effective value of the D-eS coefficients
  !! which appears often in many expressions. This routine is for the
  !! noncollinear case.
  !
  USE kinds,            ONLY: DP
  USE ions_base,        ONLY: nsp, nat, ityp
  USE spin_orb,         ONLY: lspinorb
  USE noncollin_module, ONLY: npol
  USE uspp,             ONLY: deeq_nc, qq_at, qq_so, okvan
  USE uspp_param,       ONLY: nhm
  USE lsda_mod,         ONLY: nspin
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: et
  !! The eigenvalues of the hamiltonian
  COMPLEX(DP), INTENT(OUT) :: deff(nhm,nhm,nat,nspin) 
  !! Effective values of the D-eS coefficients
  !
  ! ... local variables
  !
  INTEGER :: nt, na, is, js, ijs
  !
  deff=deeq_nc
  IF (okvan) THEN
     !
     DO nt = 1, nsp
        DO na = 1, nat
           !
           IF ( ityp(na) == nt ) THEN
              IF (lspinorb) THEN
                 deff(:,:,na,:) = deff(:,:,na,:) - et * qq_so(:,:,:,nt)
              ELSE
                 ijs=0
                 !
                 DO is=1,npol
                    DO js=1,npol
                       !
                       ijs=ijs+1
                       IF (is==js) deff(:,:,na,ijs)=deff(:,:,na,ijs)-et*qq_at(:,:,na)
                       !
                    ENDDO
                 ENDDO
                 !
              ENDIF
           ENDIF
           !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE my_compute_deff_nc

