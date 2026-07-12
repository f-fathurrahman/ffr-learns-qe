FUNCTION my_exxenergy2_k_qe52()
  !
  use kinds, only : dp
  USE constants, ONLY : fpi, e2, pi
  USE io_files, ONLY : iunwfc_exx, nwordwfc
  USE buffers, ONLY : get_buffer
  USE cell_base, ONLY : alat, omega, bg, at, tpiba
  USE symm_base, ONLY : nsym, s
  USE gvect, ONLY : ngm, gstart, g
  USE wvfct, ONLY : nbnd, npwx, wg, current_k
  USE control_flags, ONLY : gamma_only
  USE wavefunctions, ONLY : evc
  USE klist, ONLY : xk, ngk, nks, nkstot
  USE lsda_mod, ONLY : lsda, current_spin, isk
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE mp, ONLY : mp_sum
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvect, ONLY : ecutrho
  USE klist, ONLY : wk
  use exx, only: dfftt, exxbuff, exxalfa, eps_occ, x_occupation, gt
  use exx_base, only: nqs, index_sym, index_xk, xkq_collect, index_xkq
  USE exx_band, ONLY : nwordwfc_exx, igk_exx, evc_exx
  !
  IMPLICIT NONE
  !
  REAL(DP) :: my_exxenergy2_k_qe52
  !
  ! local variables
  REAL(DP) :: energy 
  COMPLEX(DP), ALLOCATABLE :: temppsic(:)
  COMPLEX(DP), ALLOCATABLE :: temppsic_nc(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhoc(:)
  REAL(DP),    ALLOCATABLE :: fac(:)
  INTEGER  :: jbnd, ibnd, ik, ikk, ig, ikq, iq, ir
  integer :: npw
  INTEGER  :: h_ibnd, nrxxs, current_ik, ibnd_loop_start
  REAL(DP) :: x1, x2
  REAL(DP) :: xkq(3), xkp(3), vc
  ! temp array for vcut_spheric
  INTEGER, EXTERNAL :: global_kpoint_index
  !
  COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
  INTEGER, ALLOCATABLE :: igkt(:)
  !

  nrxxs = dfftt%nnr
  ALLOCATE( fac(dfftt%ngm) )

  !
  ALLOCATE(temppsic(nrxxs)) 
  ALLOCATE( rhoc(nrxxs) )
  !
  energy = 0.0_DP
  !
  !IF ( nks > 1 ) REWIND( iunigk )
  !
  DO ikk = 1,nks
    current_ik = global_kpoint_index( nkstot, ikk )
    xkp = xk(:,ikk)
    !
    IF( lsda ) THEN
      current_spin = isk(ikk)
    ENDIF
    npw = ngk(ikk)
    IF ( nks > 1 ) THEN
      CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ikk )
    END IF
    !
    !
    DO jbnd = 1, nbnd     !for each band of psi (the k cycle is outside band)
      !
      IF( ABS(wg(jbnd,ikk)) < eps_occ ) THEN
        CYCLE
      ENDIF
      temppsic = 0.0_DP
      DO ig = 1, npw
        temppsic(dfftt%nl(igk_exx(ig,current_ik))) = evc_exx(ig,jbnd)
      ENDDO
      CALL invfft('Wave', temppsic, dfftt)
      !
      DO iq = 1,nqs
        !
        ikq = index_xkq(current_ik,iq)
        ik = index_xk(ikq)
        xkq = xkq_collect(:,ikq)
        CALL my_g2_convolution_all( dfftt%ngm, gt, xkp, xkq, iq, current_k )
        !
        DO ibnd = 1, nbnd
          !
          IF(ABS(x_occupation(ibnd,ik)) < eps_occ) THEN
            CYCLE
          ENDIF
          !
          ! load the phi at this k+q and band
          !calculate rho in real space
          DO ir = 1, nrxxs
            rhoc(ir) = CONJG(exxbuff(ir,ibnd,ikq))*temppsic(ir) / omega
          ENDDO
          !
          ! bring rhoc to G-space
          CALL fwfft( 'Rho', rhoc, dfftt )
          !
          vc = 0.0_DP
          DO ig = 1,dfftt%ngm
            vc = vc + fac(ig) * DBLE(rhoc(dfftt%nl(ig)) * CONJG(rhoc(dfftt%nl(ig))))
          ENDDO
          vc = vc * omega * x_occupation(ibnd,ik) / nqs
          ! 
          energy = energy - exxalfa * vc * wg(jbnd,ikk)
          !
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  DEALLOCATE(temppsic) 
  !
  DEALLOCATE(rhoc, fac)
  !
  my_exxenergy2_k_qe52 = energy
  !
END FUNCTION
