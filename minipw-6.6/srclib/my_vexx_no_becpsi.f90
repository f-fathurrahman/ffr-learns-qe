! This is from exx module
! ffr: becpsi is removed
!--------------------------------------------
SUBROUTINE my_vexx_no_becpsi( lda, n, m, psi, hpsi )
!--------------------------------------------
  USE kinds, ONLY : DP
  USE noncollin_module, ONLY : npol
  !
  USE control_flags, ONLY : gamma_only
  USE fft_types, ONLY : fft_type_descriptor
  USE stick_base, ONLY : sticks_map, sticks_map_deallocate
  !
  !! Wrapper routine computing V_x\psi, V_x = exchange potential. 
  !! Calls generic version vexx_k or Gamma-specific one vexx_gamma.
  !
  USE becmod, ONLY : bec_type
  USE mp_exx, ONLY : negrp, init_index_over_band
  USE exx_band, ONLY : transform_psi_to_exx, transform_hpsi_to_local
  !
  IMPLICIT NONE
  !
  INTEGER :: lda
  !! input: leading dimension of arrays psi and hpsi
  INTEGER :: n
  !! input: true dimension of psi and hpsi
  INTEGER :: m
  !! input: number of states psi
  COMPLEX(DP) :: psi(lda*npol,m)
  !! input: m wavefunctions
  COMPLEX(DP) :: hpsi(lda*npol,m)
  !! output: V_x*psi
  
  ! !!! negrp stuffs are removed
  if(negrp > 1) then
    stop 'negrp > 1 is disabled here'
  endif
  !
  ! calculate the EXX contribution to hpsi
  IF ( gamma_only ) THEN
    stop 'This is not supported, label 48 in my_vexx_no_becpsi'
  ELSE
    !CALL my_vexx_k_no_becpsi( lda, n, m, psi, hpsi )
    CALL my_vexx_k_no_becpsi_qe52( lda, n, m, psi, hpsi )
  ENDIF
  !
END SUBROUTINE



!-------------------------------------------------------
SUBROUTINE my_vexx_k_no_becpsi( lda, n, m, psi, hpsi )
!-------------------------------------------------------
  USE kinds, ONLY : DP
  USE noncollin_module, ONLY : npol, noncolin
  !
  USE fft_types, ONLY : fft_type_descriptor
  USE stick_base, ONLY : sticks_map, sticks_map_deallocate
  !
  ! within my_exx module
  USE exx, only: exxbuff, x_occupation, exxalfa, dfftt, eps_occ, gt
  !! Generic, k-point version of vexx.
  !
  USE constants,      ONLY : fpi, e2, pi
  USE cell_base,      ONLY : omega
  USE wvfct,          ONLY : npwx, current_k
  USE klist,          ONLY : xk, nkstot
  USE fft_interfaces, ONLY : fwfft, invfft
  USE becmod,         ONLY : bec_type
  USE mp_exx,         ONLY : inter_egrp_comm, my_egrp_id, negrp, &
                              intra_egrp_comm, me_egrp, &
                              max_pairs, egrp_pairs, ibands, nibands, &
                              max_ibands, iexx_istart, iexx_iend, &
                              all_start, all_end, iexx_start, jblock
  USE mp,             ONLY : mp_sum, mp_barrier, mp_circular_shift_left
  USE uspp,           ONLY : nkb, okvan
  USE us_exx,         ONLY : bexg_merge, addusxx_g, addusxx_r, &
                              newdxx_g, newdxx_r, add_nlxx_pot, &
                              qvan_init, qvan_clean
  USE paw_exx,        ONLY : PAW_newdxx
  USE exx_base,       ONLY : nqs, xkq_collect, index_xkq, index_xk, &
                              coulomb_fac
  USE exx_band,       ONLY : result_sum, igk_exx
  !
  !
  IMPLICIT NONE
  !
  INTEGER :: lda
  !! input: leading dimension of arrays psi and hpsi
  INTEGER :: n
  !! input: true dimension of psi and hpsi
  INTEGER :: m
  !! input: number of states psi
  COMPLEX(DP) :: psi(lda*npol,max_ibands)
  !! input: m wavefunctions
  COMPLEX(DP) :: hpsi(lda*npol,max_ibands)
  !! output: V_x*psi
  !
  !TYPE(bec_type), OPTIONAL :: becpsi  ! or call a calbec(...psi) instead
  !! input: <beta|psi>, optional but needed for US and PAW case
  !
  ! local variables
  !
  COMPLEX(DP),ALLOCATABLE :: temppsic(:,:), result(:,:)
  COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:,:),result_nc(:,:,:)
  !INTEGER :: request_send, request_recv
  !
  COMPLEX(DP),ALLOCATABLE :: deexx(:,:)
  COMPLEX(DP),ALLOCATABLE,TARGET :: rhoc(:,:), vc(:,:)

  REAL(DP),   ALLOCATABLE :: fac(:), facb(:)
  INTEGER :: ibnd, ik, im , ikq, iq
  INTEGER :: ir, ig, ir_start, ir_end
  INTEGER :: irt, nrt, nblock
  INTEGER :: current_ik
  !INTEGER :: ibnd_loop_start
  INTEGER :: nrxxs
  REAL(DP) :: xkp(3), omega_inv, nqs_inv
  REAL(DP) :: xkq(3)
  INTEGER, EXTERNAL :: global_kpoint_index
  DOUBLE PRECISION :: max
  COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
  INTEGER :: ipair, jbnd
  INTEGER :: ii, jstart, jend, jcount
  INTEGER :: ialloc, ending_im
  INTEGER :: ijt, njt, jblock_start, jblock_end
  INTEGER :: iegrp, wegrp
  !
  ialloc = nibands(my_egrp_id+1) !ffr: what's this?
  !
  ALLOCATE( fac(dfftt%ngm) )
  nrxxs= dfftt%nnr
  ALLOCATE( facb(nrxxs) )
  !
  IF (noncolin) THEN
    ALLOCATE( temppsic_nc(nrxxs,npol,ialloc), result_nc(nrxxs,npol,ialloc) )
  ELSE
    ALLOCATE( temppsic(nrxxs,ialloc), result(nrxxs,ialloc) )
  ENDIF
  !
  current_ik = global_kpoint_index( nkstot, current_k )
  xkp = xk(:,current_k)
  !
  ALLOCATE( big_result(n*npol,m) )
  big_result = 0.0_DP
  !
  ! allocate arrays for rhoc and vc
  ALLOCATE( rhoc(nrxxs,jblock), vc(nrxxs,jblock) )

  write(*,*) 'ENTER my_vexx_k_no_becpsi' 

  write(*,*) 'my_egrp_id = ', my_egrp_id
  write(*,*) 'nibands(my_egrp_id+1) = ', nibands(my_egrp_id+1)
  !
  !stop 'early stop 153 in my_vexx_k_no_becpsi'
  !
  !
  DO ii = 1, nibands(my_egrp_id+1)
    !
    ibnd = ibands(ii,my_egrp_id+1)
    write(*,*) 'ii, ibnd = ', ii, ibnd
    !
    IF (ibnd==0 .OR. ibnd>m) CYCLE
    !
    IF (noncolin) THEN
      temppsic_nc(:,:,ii) = 0._DP
    ELSE
      DO ir = 1, nrxxs
        temppsic(ir,ii) = 0._DP
      ENDDO
    ENDIF
    !
    IF (noncolin) THEN
        !
        DO ig = 1, n
          temppsic_nc(dfftt%nl(igk_exx(ig,current_k)),1,ii) = psi(ig,ii)
          temppsic_nc(dfftt%nl(igk_exx(ig,current_k)),2,ii) = psi(npwx+ig,ii)
        ENDDO
        !
        CALL invfft( 'Wave', temppsic_nc(:,1,ii), dfftt )
        CALL invfft( 'Wave', temppsic_nc(:,2,ii), dfftt )
        !
      ELSE
        !
        DO ig = 1, n
          temppsic( dfftt%nl(igk_exx(ig,current_k)), ii ) = psi(ig,ii)
        ENDDO
        !
        CALL invfft( 'Wave', temppsic(:,ii), dfftt )
        !
      ENDIF
      !
      IF (noncolin) THEN
        DO ir = 1, nrxxs
          result_nc(ir,1,ii) = 0.0_DP
          result_nc(ir,2,ii) = 0.0_DP
        ENDDO
      ELSE
        DO ir = 1, nrxxs
          result(ir,ii) = 0.0_DP
        ENDDO
      ENDIF
      !
  ENDDO
  
  !stop 'early stop 204 my_vexx_k_no_becpsi'

  !
  !precompute these guys
  omega_inv = 1.0 / omega
  nqs_inv = 1.0 / nqs
  !
  !------------------------------------------------------------------------!
  ! Beginning of main loop
  !------------------------------------------------------------------------!
  DO iq = 1, nqs
    !
    ikq = index_xkq(current_ik,iq)
    ik  = index_xk(ikq)
    xkq = xkq_collect(:,ikq)
    !
    ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
    CALL my_g2_convolution_all( dfftt%ngm, gt, xkp, xkq, iq, current_k )
    !
    ! JRD - below not threaded
    facb = 0D0
    DO ig = 1, dfftt%ngm
      facb(dfftt%nl(ig)) = coulomb_fac(ig,iq,current_k)
    ENDDO
    !
    DO iegrp = 1, negrp
      write(*,*) 'iq, ikq, ik, iegrp, njt = ', iq, ikq, ik, iegrp
      !
      ! compute the id of group whose data is currently worked on
      wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
      njt = (all_end(wegrp)-all_start(wegrp)+jblock)/jblock
      write(*,*) 'wegrp, njt = ', wegrp, njt
      !
      DO ijt = 1, njt
        !
        jblock_start = (ijt - 1) * jblock + all_start(wegrp)
        jblock_end = MIN(jblock_start+jblock-1,all_end(wegrp))
        write(*,*) 'jblock_start, jblock_end = ', jblock_start, jblock_end
        !
        DO ii = 1, nibands(my_egrp_id+1)
          !
          ibnd = ibands(ii,my_egrp_id+1)
          !
          IF (ibnd==0 .OR. ibnd>m) CYCLE
          !
          !determine which j-bands to calculate
          jstart = 0
          jend = 0
          !
          DO ipair = 1, max_pairs
            IF (egrp_pairs(1,ipair,my_egrp_id+1) == ibnd) THEN
              IF (jstart == 0) THEN
                jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                jend = jstart
              ELSE
                jend = egrp_pairs(2,ipair,my_egrp_id+1)
              ENDIF
            ENDIF
          ENDDO
          !
          jstart = MAX( jstart, jblock_start )
          jend = MIN( jend, jblock_end )
          !
          !how many iters
          jcount = jend-jstart+1
          IF (jcount <= 0) CYCLE
          !
          !----------------------------------------------------------------------!
          !INNER LOOP START
          !----------------------------------------------------------------------!
          !
          nblock = 2048
          nrt = (nrxxs+nblock-1)/nblock
          !
          DO irt = 1, nrt
            DO jbnd = jstart, jend
              ir_start = (irt - 1) * nblock + 1
              ir_end = MIN(ir_start+nblock-1, nrxxs)
              IF (noncolin) THEN
                DO ir = ir_start, ir_end
                  rhoc(ir,jbnd-jstart+1) = ( CONJG(exxbuff(ir,jbnd-all_start(wegrp)+ &
                                              iexx_start,ikq))*temppsic_nc(ir,1,ii) + &
                                              CONJG(exxbuff(nrxxs+ir,jbnd-all_start(wegrp)+ &
                                              iexx_start,ikq))*temppsic_nc(ir,2,ii) )/omega
                ENDDO
              ELSE
                DO ir = ir_start, ir_end
                  rhoc(ir,jbnd-jstart+1) = CONJG(exxbuff(ir,jbnd-all_start(wegrp)+ &
                                            iexx_start,ikq))*temppsic(ir,ii)*omega_inv
                ENDDO
              ENDIF
            ENDDO
          ENDDO
          !
          !   brings it to G-space
          DO jbnd=jstart, jend
            CALL fwfft( 'Rho', rhoc(:,jbnd-jstart+1), dfftt )
          ENDDO
          !
          DO irt = 1, nrt
            DO jbnd = jstart, jend
              ir_start = (irt - 1) * nblock + 1
              ir_end = MIN(ir_start+nblock-1,nrxxs)
              DO ir = ir_start, ir_end
                vc(ir,jbnd-jstart+1) = facb(ir) * rhoc(ir,jbnd-jstart+1)*x_occupation(jbnd,ik) * nqs_inv
              ENDDO
            ENDDO
          ENDDO
          ! brings back v in real space
          DO jbnd = jstart, jend
            CALL invfft( 'Rho', vc(:,jbnd-jstart+1), dfftt )
          ENDDO
          !
          ! accumulates over bands and k points
          !
          DO irt = 1, nrt
            DO jbnd = jstart, jend
              ir_start = (irt - 1) * nblock + 1
              ir_end = MIN(ir_start+nblock-1, nrxxs)
              IF (noncolin) THEN
                DO ir = ir_start, ir_end
                  result_nc(ir,1,ii) = result_nc(ir,1,ii) + vc(ir,jbnd-jstart+1) * &
                                        exxbuff(ir,jbnd-all_start(wegrp)+iexx_start,ikq)
                  result_nc(ir,2,ii) = result_nc(ir,2,ii) + vc(ir,jbnd-jstart+1) * &
                                        exxbuff(ir+nrxxs,jbnd-all_start(wegrp)+iexx_start,ikq)
                ENDDO
              ELSE
                DO ir = ir_start, ir_end
                  result(ir,ii) = result(ir,ii) + vc(ir,jbnd-jstart+1)* &
                                  exxbuff(ir,jbnd-all_start(wegrp)+iexx_start,ikq)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
          !
          !----------------------------------------------------------------------!
          !INNER LOOP END
          !----------------------------------------------------------------------!
          !
        ENDDO !I-LOOP
      ENDDO !IJT
    ENDDO !iegrp
  ENDDO
  !
  stop 'early stop 356 my_vexx_k_no_becpsi'
  !
  DO ii = 1, nibands(my_egrp_id+1)
      !
      ibnd = ibands(ii,my_egrp_id+1)
      !
      IF (ibnd==0 .OR. ibnd>m) CYCLE
      !
      IF (okvan) THEN
        CALL mp_sum( deexx(:,ii), intra_egrp_comm )
      ENDIF
      !
      IF (noncolin) THEN
        !brings back result in G-space
        CALL fwfft( 'Wave', result_nc(:,1,ii), dfftt )
        CALL fwfft( 'Wave', result_nc(:,2,ii), dfftt )
        !
        DO ig = 1, n
            big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa* &
                                  result_nc(dfftt%nl(igk_exx(ig,current_k)),1,ii)
            big_result(n+ig,ibnd) = big_result(n+ig,ibnd) - exxalfa* &
                                    result_nc(dfftt%nl(igk_exx(ig,current_k)),2,ii)
        ENDDO
      ELSE
        !
        CALL fwfft( 'Wave', result(:,ii), dfftt )
        !
        DO ig = 1, n
            big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa* &
                                  result(dfftt%nl(igk_exx(ig,current_k)),ii)
        ENDDO
      ENDIF
      !
      ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
      IF (okvan) CALL add_nlxx_pot( lda, big_result(:,ibnd), xkp, n, igk_exx(:,current_k), &
                                  deexx(:,ii), eps_occ, exxalfa )
      !
  ENDDO
  !
  !deallocate temporary arrays
  DEALLOCATE( rhoc, vc )
  !
  !sum result
  CALL result_sum( n*npol, m, big_result )
  !
  IF (iexx_istart(my_egrp_id+1) > 0) THEN
      !
      IF (negrp == 1) THEN
        ending_im = m
      ELSE
        ending_im = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
      ENDIF
      !
      IF (noncolin) THEN
        DO im = 1, ending_im
            DO ig = 1, n
              hpsi(ig,im) = hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
            ENDDO
            DO ig = 1, n
              hpsi(lda+ig,im) = hpsi(lda+ig,im) + big_result(n+ig,im+iexx_istart(my_egrp_id+1)-1)
            ENDDO
        ENDDO
      ELSE
        DO im = 1, ending_im
            DO ig = 1, n
              hpsi(ig,im) = hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
            ENDDO
        ENDDO
      ENDIF
  ENDIF
  !
  IF (noncolin) THEN
    DEALLOCATE( temppsic_nc, result_nc )
  ELSE
    DEALLOCATE( temppsic, result )
  ENDIF
  !
  DEALLOCATE( big_result )
  DEALLOCATE( fac, facb )
  IF (okvan) DEALLOCATE( deexx )
  !
END SUBROUTINE
