!-----------------------------------------------------------------------
SUBROUTINE my_vexx_k_no_becpsi_qe52(lda, n, m, psi, hpsi)
!-----------------------------------------------------------------------
  use kinds, only: dp
  !
  ! generic, k-point version of vexx
  !
  USE constants,      ONLY : fpi, e2, pi
  USE cell_base,      ONLY : omega
  USE wvfct,          ONLY : npwx, current_k, nbnd
  USE klist,          ONLY : xk, nkstot
  USE fft_interfaces, ONLY : fwfft, invfft
  USE mp_bands,       ONLY : my_bgrp_id, nbgrp, inter_bgrp_comm
  USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast
  use noncollin_module, only: npol, noncolin
  USE exx, only: exxbuff, x_occupation, exxalfa, dfftt, eps_occ, gt
  USE exx_base, ONLY : nqs, xkq_collect, index_xkq, index_xk, coulomb_fac
  USE exx_band, ONLY : result_sum, igk_exx
  !
  !
  IMPLICIT NONE
  !
  INTEGER :: lda, n, m
  COMPLEX(DP) :: psi(lda*npol,m) 
  COMPLEX(DP) :: hpsi(lda*npol,m)
  !
  ! local variables
  COMPLEX(DP), ALLOCATABLE :: temppsic(:), result(:)
  COMPLEX(DP), ALLOCATABLE :: temppsic_nc(:,:),result_nc(:,:)
  !
  COMPLEX(DP), ALLOCATABLE :: rhoc(:), vc(:)
  REAL(DP), ALLOCATABLE :: fac(:)
  INTEGER :: ibnd, ik, im , ikq, iq
  INTEGER :: ir, ig
  INTEGER :: current_ik
  INTEGER :: nrxxs
  REAL(DP) :: xkp(3)
  REAL(DP) :: xkq(3)
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  !
  nrxxs = dfftt%nnr
  !
  IF (noncolin) THEN
    ALLOCATE( temppsic_nc(nrxxs,npol), result_nc(nrxxs,npol) )
  ELSE
    ALLOCATE( temppsic(nrxxs), result(nrxxs) )
  ENDIF
  !
  !ALLOCATE(fac(nrxxs))
  ALLOCATE(fac(dfftt%ngm))
  ALLOCATE(rhoc(nrxxs), vc(nrxxs))
  !
  current_ik = global_kpoint_index( nkstot, current_k )
  xkp = xk(:,current_ik)
  !
  ! This is to stop numerical inconsistencies creeping in through the band parallelization.
  !
  IF(my_bgrp_id > 0) THEN
    hpsi = 0.0_DP
    psi = 0.0_DP
  ENDIF
  IF(nbgrp > 1) THEN
    CALL mp_bcast(hpsi, 0, inter_bgrp_comm)
    CALL mp_bcast(psi, 0, inter_bgrp_comm)
  ENDIF
  !
  LOOP_ON_PSI_BANDS : &
  DO im = 1,m ! for each band of psi (the k cycle is outside band)
    !ffr: m is number of columns of input psi
    !
    IF (noncolin) THEN
      temppsic_nc = 0._DP
    ELSE
      temppsic = 0.0_DP
    ENDIF
    !
    !ffr: Bring psi(:,im) to real space using invfft
    IF (noncolin) THEN
      DO ig = 1, n
        temppsic_nc(dfftt%nl(igk_exx(ig,current_k)),1) = psi(ig,im)
        temppsic_nc(dfftt%nl(igk_exx(ig,current_k)),2) = psi(npwx+ig,im)
      ENDDO
      !ffr: im is band index
      !
      CALL invfft('Wave', temppsic_nc(:,1), dfftt)
      CALL invfft('Wave', temppsic_nc(:,2), dfftt)
    ELSE
      DO ig = 1, n
        temppsic(dfftt%nl(igk_exx(ig,current_k))) = psi(ig,im)
      ENDDO
      !
      CALL invfft( 'Wave', temppsic, dfftt )
    ENDIF
    !
    ! Zero out the results
    IF (noncolin) THEN
      result_nc = 0.0_DP
    ELSE
      result = 0.0_DP
    ENDIF
    !
    INTERNAL_LOOP_ON_Q : &
    DO iq = 1,nqs
      !
      ikq  = index_xkq(current_ik,iq)
      ik   = index_xk(ikq)
      xkq  = xkq_collect(:,ikq)
      ! calculate the 1/|r-r'| (actually, k+q+g) factor, result is in global variable coulomb_fac
      fac = 0.d0 ! need this?
      CALL my_g2_convolution( dfftt%ngm, gt, xkp, xkq, fac )
      !DO ig = 1, dfftt%ngm
      !  fac(dfftt%nl(ig)) = coulomb_fac(ig,iq,current_k)
      !ENDDO
      !
      IBND_LOOP_K : &
      DO ibnd = 1, nbnd !for each band of psi
        !write(*,*) 'ibnd = ', ibnd
        !
        IF( ABS(x_occupation(ibnd,ik)) < eps_occ) THEN 
          CYCLE IBND_LOOP_K
        ENDIF
        !
        !loads the phi from file
        !
        !   >>>> calculate rho in real space
        IF (noncolin) THEN
          DO ir = 1, nrxxs
            rhoc(ir) = ( CONJG(exxbuff(ir,ibnd,ikq))*temppsic_nc(ir,1) + &
                    &    CONJG(exxbuff(nrxxs+ir,ibnd,ikq))*temppsic_nc(ir,2) )/omega
          ENDDO
        ELSE
          DO ir = 1, nrxxs
            rhoc(ir) = CONJG(exxbuff(ir,ibnd,ikq))*temppsic(ir) / omega
          ENDDO
        ENDIF
        !
        !   >>>> brings it to G-space
        CALL fwfft( 'Rho', rhoc, dfftt )
        !   >>>> charge done
        !
        vc = 0._DP
        !
        DO ig = 1, dfftt%ngm
          vc(dfftt%nl(ig)) = fac(ig) * rhoc(dfftt%nl(ig)) * x_occupation(ibnd,ik) / nqs
        ENDDO
        !ffr: multiply by points?
        !DO ir = 1,nrxxs
        !  vc(ir) = fac(ir) * rhoc(ir)*x_occupation(ibnd,ik)/nqs
        !ENDDO
        !
        !brings back v in real space
        CALL invfft('Rho', vc, dfftt)
        !
        ! accumulates over bands and k points
        !
        IF (noncolin) THEN
          DO ir = 1, nrxxs
            result_nc(ir,1)= result_nc(ir,1) + vc(ir) * exxbuff(ir,ibnd,ikq)
          ENDDO
          DO ir = 1, nrxxs
            result_nc(ir,2)= result_nc(ir,2) + vc(ir) * exxbuff(ir+nrxxs,ibnd,ikq)
          ENDDO
        ELSE
          DO ir = 1, nrxxs
            result(ir) = result(ir) + vc(ir)*exxbuff(ir,ibnd,ikq)
          ENDDO
        ENDIF
          !
      ENDDO IBND_LOOP_K
      !
    ENDDO INTERNAL_LOOP_ON_Q
    !
    IF (noncolin) THEN
      CALL mp_sum( result_nc(1:nrxxs,1:npol), inter_bgrp_comm)
    ELSE
      CALL mp_sum( result(1:nrxxs), inter_bgrp_comm)
    ENDIF
    !
    !brings back result in G-space
    !
    IF (noncolin) THEN
      !brings back result in G-space
      CALL fwfft( 'Wave', result_nc(:,1), dfftt )
      CALL fwfft( 'Wave', result_nc(:,2), dfftt )
      !
      ! adds it to hpsi
      DO ig = 1, n
          hpsi(ig,im) = hpsi(ig,im) - exxalfa*result_nc(dfftt%nl(igk_exx(ig,current_k)),1)
      ENDDO
      DO ig = 1, n
          hpsi(lda+ig,im)= hpsi(lda+ig,im) - exxalfa*result_nc(dfftt%nl(igk_exx(ig,current_k)),2)
      ENDDO
        !
    ELSE
      !
      CALL fwfft('Wave', result, dfftt)
      !
      !adds it to hpsi
      DO ig = 1, n
        hpsi(ig,im) = hpsi(ig,im) - exxalfa*result(dfftt%nl(igk_exx(ig,current_k)))
      ENDDO
    ENDIF
    !
  ENDDO LOOP_ON_PSI_BANDS
  
  IF (noncolin) THEN
    DEALLOCATE(temppsic_nc, result_nc) 
  ELSE
    DEALLOCATE(temppsic, result) 
  ENDIF
  !
  DEALLOCATE(rhoc, vc, fac)
  !
  !flush(6) ! stdout
  !flush(0) ! stderr
  !stop 'Early stop 210 in my_vexx_k_no_becpsi_qe52'
  !
END SUBROUTINE
