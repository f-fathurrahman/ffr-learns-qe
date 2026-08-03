!-----------------------------------------------------------------------------
 SUBROUTINE my_aceinit( DoLoc, exex )
!----------------------------------------------------------------------------
  use kinds, only: dp
  use control_flags, only: gamma_only
  use noncollin_module, only: npol
  !! ACE Initialization
  !
  USE wvfct, ONLY : nbnd, npwx, current_k
  USE klist, ONLY : nks, xk, ngk, igk_k
  USE uspp, ONLY : nkb, vkb, okvan
  USE becmod, ONLY : allocate_bec_type, deallocate_bec_type, &
                                bec_type, calbec
  USE lsda_mod, ONLY : current_spin, lsda, isk
  USE io_files, ONLY : nwordwfc, iunwfc
  USE buffers, ONLY : get_buffer
  USE mp, ONLY : mp_sum
  USE wavefunctions, ONLY : evc
  !
  use exx, only: xi, domat, x_nbnd_occ, nbndproj, aceinit_gamma
  !
  IMPLICIT NONE
  !
  LOGICAL :: DoLoc
  !! if TRUE calculates exact exchange with SCDM orbitals
  REAL(DP) :: exex !ffr: This is originally
  !! ACE energy
  !
  ! local variables
  !
  REAL(DP) :: ee, eexx
  INTEGER :: ik, npw
  TYPE(bec_type) :: becpsi

  write(*,*)
  write(*,*) '<div> ENTER my_aceinit'
  write(*,*)

  !
  IF (nbndproj < x_nbnd_occ .OR. nbndproj > nbnd) THEN 
    WRITE(*, '(3(A,I4))') ' occ = ', x_nbnd_occ, ' proj = ', nbndproj, ' tot = ', nbnd
    CALL errore( 'aceinit', 'n_proj must be between occ and tot.', 1 )
  ENDIF
  !
  IF (.NOT. ALLOCATED(xi)) THEN 
    ALLOCATE( xi(npwx*npol,nbndproj,nks) )
  ENDIF
  !
  IF( okvan ) THEN 
    CALL allocate_bec_type( nkb, nbnd, becpsi )
  ENDIF
  !
  eexx = 0.0d0
  xi(:,:,:) = (0.0d0, 0.0d0)
  !
  DO ik = 1, nks
    npw = ngk(ik)
    current_k = ik
    IF ( lsda ) THEN 
      current_spin = isk(ik)
    ENDIF
    IF ( nks > 1 ) THEN 
      CALL get_buffer( evc, nwordwfc, iunwfc, ik )
    ENDIF
    IF ( okvan ) THEN
      CALL init_us_2( npw, igk_k(1,ik), xk(:,ik), vkb )
      CALL calbec( npw, vkb, evc, becpsi, nbnd )
    ENDIF
    
    IF (gamma_only) THEN
      CALL aceinit_gamma( DoLoc, npw, nbnd, evc, xi(1,1,ik), becpsi, ee )
    ELSE
      CALL my_aceinit_k( DoLoc, npw, nbnd, evc, xi(1,1,ik), becpsi, ee )
    ENDIF
    eexx = eexx + ee
  ENDDO
  !
  WRITE(*,'(/,5X,"ACE energy",f15.8)') eexx
  !
  !IF( PRESENT(exex) ) THEN 
    exex = eexx
  !ENDIF
  IF( okvan ) THEN 
    CALL deallocate_bec_type( becpsi )
  ENDIF
  !
  domat = .FALSE.
  
  write(*,*)
  write(*,*) '</div> EXIT my_aceinit'
  write(*,*)

  return
  !
END SUBROUTINE

!---------------------------------------------------------------------------------------------
SUBROUTINE my_aceinit_k( DoLoc, nnpw, nbnd, phi, xitmp, becpsi, exxe )
!-----------------------------------------------------------------------------------------
  !! Compute xi(npw,nbndproj) for the ACE method.
  !
  use kinds, only: dp
  !
  USE becmod,               ONLY : bec_type
  USE wvfct,                ONLY : current_k, npwx
  USE klist,                ONLY : wk
  USE noncollin_module,     ONLY : npol
  !
  use exx, only: nbndproj, vexx, domat, evc0, vexx_loc_k, vexxace_k
  !
  IMPLICIT NONE
  !
  LOGICAL :: DoLoc
  !! if TRUE calculates exact exchange with SCDM orbitals
  INTEGER :: nnpw
  !! number of PW
  INTEGER :: nbnd
  !! number of bands
  COMPLEX(DP) :: phi(npwx*npol,nbnd)
  !! wave function
  COMPLEX(DP) :: xitmp(npwx*npol,nbndproj)
  !! xi(nnpw,nbndproj)
  TYPE(bec_type), INTENT(IN) :: becpsi
  !! <beta|psi>
  REAL(DP) :: exxe
  !! exx energy
  !
  ! local variables
  !
  COMPLEX(DP), ALLOCATABLE :: mexx(:,:)
  REAL(DP) :: exxe0
  REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0
  INTEGER :: i
  LOGICAL :: domat0

  IF(nbndproj > nbnd) THEN 
    CALL errore( 'aceinit_k', 'nbndproj greater than nbnd.', 1 )
  ENDIF
  !
  IF(nbndproj<=0) THEN 
    CALL errore( 'aceinit_k', 'nbndproj le 0.', 1 )
  ENDIF
  !
  ALLOCATE( mexx(nbndproj,nbndproj) )
  xitmp = (Zero,Zero)
  mexx  = (Zero,Zero)
  IF ( DoLoc ) THEN
    CALL vexx_loc_k( nnpw, nbndproj, xitmp, mexx, exxe )
    CALL MatSymm_k( 'S', 'L', mexx, nbndproj )
  ELSE
    !write(*,'(1x,A)',advance='no') 'Will call vexx ...'
    ! |xi> = Vx[phi]|phi>
    !CALL vexx( npwx, nnpw, nbndproj, phi, xitmp, becpsi )
    write(*,*) ' ... done calling vexx'
    ! mexx = <phi|Vx[phi]|phi>
    CALL matcalc_k( 'exact', .TRUE., 0, current_k, npwx*npol, nbndproj, nbndproj, &
                    phi, xitmp, mexx, exxe )
  ENDIF
  !
  WRITE(*,'(3(A,I3),A,I9,A,f12.6)') 'aceinit_k: nbnd=', nbnd, ' nbndproj=', nbndproj, ' k=',current_k,' npw=',nnpw,' Ex(k)=',exxe
  ! Skip k-points that have exactly zero weight
  IF( wk(current_k) /= 0._dp )THEN
    ! |xi> = -One * Vx[phi]|phi> * rmexx^T
    CALL my_aceupdate_k( nbndproj, nnpw, xitmp, mexx )
  ENDIF
  !
  DEALLOCATE( mexx )
  !
  IF ( DoLoc ) THEN
    domat0 = domat
    domat = .TRUE.
    CALL vexxace_k( nnpw, nbnd, evc0(1,1,current_k), exxe )
    !CALL my_vexxace_k( nnpw, nbnd, evc0(1,1,current_k), exxe )
    evc0(:,:,current_k) = phi(:,:)
    domat = domat0
  ENDIF 
  !
END SUBROUTINE


!------------------------------------------------------------------------------
SUBROUTINE my_aceupdate_k( nbndproj, nnpw, xitmp, mexx )
!----------------------------------------------------------------------------
  !! Updates xi(npw,nbndproj) for the ACE method.
  !
  use kinds, only: dp
  USE wvfct,                ONLY : npwx
  USE noncollin_module,     ONLY : npol
  !
  IMPLICIT NONE
  !
  INTEGER :: nbndproj !! number of bands
  INTEGER :: nnpw !! number of PW
  COMPLEX(DP) :: mexx(nbndproj,nbndproj) !! mexx = -(Cholesky(mexx))^-1
  COMPLEX(DP) :: xitmp(npwx*npol,nbndproj) !! |xi> = -One * Vx[phi]|phi> * mexx^T

  !
  ! mexx = -(Cholesky(mexx))^-1
  mexx = -mexx
  CALL invchol_k( nbndproj, mexx )
  !
  ! |xi> = -One * Vx[phi]|phi> * mexx^T
  CALL ZTRMM( 'R', 'L', 'C', 'N', npwx*npol, nbndproj, (1.0_dp,0.0_dp), mexx,nbndproj, &
              xitmp, npwx*npol )
END SUBROUTINE
