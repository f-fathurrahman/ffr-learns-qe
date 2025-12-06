!------------------------------------------------------------------------------------
SUBROUTINE my_rotate_wfc( npwx, npw, nstart, gstart, nbnd, psi, npol, overlap, evc, e )
  !--------------------------------------------------------------------------------
  !! Driver routine (maybe it should be an interface) for Hamiltonian 
  !! diagonalization in the subspace spanned by nstart states 
  !! psi (atomic or random wavefunctions). 
  !! Produces on output nbnd eigenvectors ( nbnd <= nstart ) in evc.
  !! Calls \(\texttt{h_psi, s_psi}\) to calculate \(H|psi\rangle\) and
  !! \(S|psi\rangle\). 
  !! It only uses an auxiliary array of the same size as psi.
  !
  USE kinds,           ONLY: DP
  USE control_flags,   ONLY: use_para_diag, gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw
  !! dimension of the matrix to be diagonalized
  INTEGER, INTENT(IN) :: npwx
  !! leading dimension of matrix psi, as declared in the calling pgm unit
  INTEGER, INTENT(IN) :: nstart
  !! input number of states
  INTEGER, INTENT(IN) :: nbnd
  !! output number of states
  INTEGER, INTENT(IN) :: gstart
  !! first G with nonzero norm
  INTEGER, INTENT(IN) :: npol
  !! number of spin polarizations
  LOGICAL, INTENT(IN) :: overlap
  !! if .FALSE. \(S|psi\rangle\) not needed
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nstart)
  !! I/O eigenvectors (may overlap)
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx*npol,nbnd)
  !! I/O eigenvectors (may overlap)
  REAL(DP), INTENT(OUT) :: e(nbnd) !! eigenvalues
  EXTERNAL  h_psi, s_psi
  EXTERNAL  my_h_psi, my_s_psi
  ! h_psi(npwx,npw,nvec,psi,hpsi)
  !     calculates H|psi>
  ! s_psi(npwx,npw,nvec,spsi)
  !     calculates S|psi> (if needed)
  !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)
  !
  IF( use_para_diag ) THEN
    stop 'use_para_diag is disable in my_rotate_wfc.f90'
  ELSE
    !
    ! use serial subroutines
    !
    IF( gamma_only ) THEN
      !
      CALL rotate_wfc_gamma( h_psi, s_psi, overlap, &
                             npwx, npw, nstart, nbnd, psi, evc, e )
      !
    ELSE
      !
      CALL my_rotate_wfc_k( my_h_psi, my_s_psi, overlap, &
                            npwx, npw, nstart, nbnd, npol, psi, evc, e )
      !
     ENDIF
     !
  ENDIF
  !
END SUBROUTINE my_rotate_wfc
