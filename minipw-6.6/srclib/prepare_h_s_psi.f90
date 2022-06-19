!-----------------------------
subroutine prepare_h_s_psi(ik)
!-----------------------------
  USE wvfct, ONLY : nbnd, current_k
  USE lsda_mod, ONLY : current_spin, lsda, isk
  USE wavefunctions, ONLY : evc
  USE io_files, ONLY : iunwfc, nwordwfc
  USE becmod, ONLY : bec_type, becp, calbec, allocate_bec_type
  USE klist, ONLY : nks, ngk, xk, ngk, igk_k
  USE mp_bands, ONLY : intra_bgrp_comm
  USE bp, ONLY : lelfield
  USE uspp, ONLY : vkb, nkb
  USE buffers, ONLY: get_buffer
  !
  implicit none
  !
  integer :: ik

  ! Set k-point, spin, kinetic energy, needed by Hpsi
  current_k = ik
  
  !
  !IF (lda_plus_u .AND. lda_plus_u_kind .EQ. 2) CALL phase_factor(ik)

  IF ( lsda ) current_spin = isk(ik)

  CALL g2_kin( ik )

  ! More stuff needed by the hamiltonian: nonlocal projectors
  IF ( nkb > 0 ) CALL my_init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
  

  ! read in wavefunctions from the previous iteration
  IF ( nks > 1 .OR. lelfield ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )

  ! ... Needed for LDA+U
  !IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
  !     CALL get_buffer ( wfcU, nwordwfcU, iunhub, ik )

  ! From my_diag
  CALL allocate_bec_type( nkb, nbnd, becp, intra_bgrp_comm )

  return

end subroutine
