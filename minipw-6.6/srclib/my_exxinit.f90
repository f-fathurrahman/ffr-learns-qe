!ffr: originally from exx module

!------------------------------------------------------------------------
SUBROUTINE my_exxinit( DoLoc )
!------------------------------------------------------------------------
  !! This subroutine is run before the first H_psi() of each iteration. 
  !! It saves the wavefunctions for the right density matrix, in real space.
  !
  !---------------------------------------------------------------
  USE kinds, ONLY : DP
  USE noncollin_module, ONLY : noncolin, npol
  !
  USE control_flags, ONLY : gamma_only
  USE fft_types, ONLY : fft_type_descriptor
  USE stick_base, ONLY : sticks_map, sticks_map_deallocate
  !
  ! within exx module
  USE exx, only: use_ace, exxbuff, npwt, x_occupation, x_nbnd_occ, locmat, &
                  & nbndproj, local_thr, locbuff, ibnd_start, ibnd_end, &
                  & ibnd_buff_start, ibnd_buff_end, exxmat, exxalfa, evc0, dfftt, &
                  & eps_occ
  USE exx, only: exx_fft_create
  !---------------------------------------------------------------
  !
  USE io_files, ONLY : iunwfc_exx
  USE buffers, ONLY : get_buffer
  USE wvfct, ONLY : nbnd, npwx, wg
  USE klist, ONLY : ngk, nks, nkstot, wk
  USE symm_base, ONLY : nsym, sr
  USE mp_exx, ONLY : init_index_over_band, inter_egrp_comm, &
                     iexx_start, iexx_end, all_start, all_end
  USE mp, ONLY : mp_sum, mp_bcast
  USE funct, ONLY : get_exx_fraction, start_exx,exx_is_active, &
                                    get_screening_parameter, get_gau_parameter
  USE scatter_mod, ONLY : gather_grid, scatter_grid
  USE fft_interfaces,       ONLY : invfft
  USE uspp,                 ONLY : okvan
  USE us_exx,               ONLY : rotate_becxx
  USE paw_variables,        ONLY : okpaw
  USE paw_exx,              ONLY : PAW_init_fock_kernel
  USE mp_orthopools,        ONLY : intra_orthopool_comm
  USE exx_base, ONLY : nkqs, xkq_collect, index_xk, index_sym,  &
                       exx_set_symm, rir, working_pool, exxdiv, &
                       erfc_scrlen, gau_scrlen, exx_divergence
  USE exx_band, ONLY : change_data_structure, nwordwfc_exx, &
                       transform_evc_to_exx, igk_exx, evc_exx
  use mp_exx, only: negrp
  !
  ! additional
  USE exx_base, only: use_coulomb_vcut_spheric, use_coulomb_vcut_ws, x_gamma_extrapolation, &
                    & exxdiv_treatment, use_regularization
  use exx_base, only: gau_scrlen, erf_scrlen, erfc_scrlen, exxdiv, yukawa, eps_qdiv


  !
  IMPLICIT NONE
  !
  LOGICAL :: DoLoc
  !! TRUE:  Real Array locbuff(ir, nbnd, nkqs);  
  !! FALSE: Complex Array exxbuff(ir, nbnd/2, nkqs).
  !
  ! local variables
  !
  INTEGER :: ik, ibnd, ir, isym, ikq, ig
  INTEGER :: ibnd_loop_start
  INTEGER :: ipol, jpol
  REAL(DP), ALLOCATABLE :: occ(:,:)
  COMPLEX(DP),ALLOCATABLE :: temppsic(:)
  COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:), psic_nc(:,:)
  COMPLEX(DP),ALLOCATABLE :: psic_exx(:)
  INTEGER :: nxxs, nrxxs
  COMPLEX(DP) :: d_spin(2,2,48)
  INTEGER :: npw, current_ik
  INTEGER, EXTERNAL :: global_kpoint_index
  INTEGER :: ibnd_start_new, ibnd_end_new, max_buff_bands_per_egrp
  INTEGER :: ibnd_exx, evc_offset

  write(*,*)
  write(*,*) '<div> ENTER my_exxinit()'
  write(*,*)

  if(negrp /= 1) then
    stop 'negrp parallezation is disabled here'
  endif

  IF ( DoLoc ) THEN
    WRITE(*,'(/,5X,"Using localization algorithm with threshold: ", D10.2)') local_thr
    ! IF (.NOT.gamma_only) CALL errore('exxinit','SCDM with K-points NYI',1)
    IF (okvan .OR. okpaw) CALL errore( 'exxinit','SCDM with USPP/PAW not implemented', 1 )
  ENDIF 
  IF( use_ace ) THEN
    WRITE(*,'(/,5X,"Using ACE for calculation of exact exchange")') 
  ENDIF
  !
  ! ffr: why need this? Only mandatory for parallelization?
  CALL transform_evc_to_exx( 2 )
  !
  ! prepare the symmetry matrices for the spin part
  !
  IF (noncolin) THEN
    DO isym = 1, nsym
      CALL find_u( sr(:,:,isym), d_spin(:,:,isym) )
    ENDDO
  ENDIF
  !
  CALL my_exx_fft_create()
  !
  ! Note that nxxs is not the same as nrxxs in parallel case
  nxxs = dfftt%nr1x * dfftt%nr2x * dfftt%nr3x
  nrxxs = dfftt%nnr
  IF(noncolin) THEN
    ALLOCATE( temppsic_nc(nrxxs, npol), psic_nc(nrxxs, npol) )
  ELSEIF( .NOT. gamma_only ) THEN
    ALLOCATE( temppsic(nrxxs) )
  ENDIF
  !
  ALLOCATE( psic_exx(nrxxs) )
  !
  !ffrL Start EXX if it is not active when this subroutine (exxinit) is called
  IF(.NOT. exx_is_active()) THEN
    ! ffr: what are these parameters?
    erfc_scrlen = get_screening_parameter()
    gau_scrlen = get_gau_parameter()
    exxdiv  = exx_divergence()
    exxalfa = get_exx_fraction()
    !
    CALL start_exx()
  ENDIF
  !ffr: what's this?
  IF( .NOT. gamma_only ) THEN
    CALL my_exx_set_symm( dfftt%nr1, dfftt%nr2,  dfftt%nr3, &
                          dfftt%nr1x, dfftt%nr2x, dfftt%nr3x )
  ENDIF
  !
  ! set occupations of wavefunctions used in the calculation of exchange term
  ! ffr: how is this different from the usual occupation numbers?
  IF(.NOT. ALLOCATED(x_occupation)) THEN
    ALLOCATE( x_occupation(nbnd,nkstot) )
  ENDIF
  !
  ALLOCATE( occ(nbnd,nks) )
  !ffr: do not include orbitals with small occupations (?)
  DO ik = 1, nks
    IF( ABS(wk(ik)) > eps_occ ) THEN
      occ(1:nbnd,ik) = wg(1:nbnd,ik) / wk(ik)
    ELSE
      occ(1:nbnd,ik) = 0._DP
    ENDIF
  ENDDO
  !ffr: wg is the weight of each band and kpt (combination of focc and wk)
  !ffr: x_occupation is now only focc
  !
  CALL poolcollect( nbnd, nks, occ, nkstot, x_occupation )
  !write(*,*) 'occ = ', occ
  !write(*,*) 'x_occupation = ', x_occupation
  !stop 'early stop in 154'
  !
  DEALLOCATE( occ )
  !
  ! find an upper bound to the number of bands with non zero occupation.
  ! Useful to distribute bands among band groups
  !
  x_nbnd_occ = 0
  DO ik = 1, nkstot
    DO ibnd = MAX(1,x_nbnd_occ), nbnd
      IF (ABS(x_occupation(ibnd,ik)) > eps_occ) x_nbnd_occ = ibnd
    ENDDO
  ENDDO
  write(*,*) 'nbnd = ', nbnd
  write(*,*) 'eps_occ = ', eps_occ
  write(*,*) 'x_nbnd_occ = ', x_nbnd_occ
  ! These should be small
  write(*,*) 'x_occupation(x_nbnd_occ:nbnd,1) = ', x_occupation(x_nbnd_occ:nbnd,1)
  !stop 'Early stop 170'
  !
  !ffr: by default we use all bands?
  IF (nbndproj == 0) nbndproj = nbnd
  write(*,*) 'nbndproj = ', nbndproj
  !stop 'Early stop 178'
  !
  CALL divide( inter_egrp_comm, x_nbnd_occ, ibnd_start, ibnd_end )
  CALL init_index_over_band( inter_egrp_comm, nbnd, nbnd )
  !
  ! this will cause exxbuff to be calculated for every band
  ibnd_start_new = iexx_start
  ibnd_end_new = iexx_end
  !
  IF ( gamma_only ) THEN
    ibnd_buff_start = (ibnd_start_new+1)/2
    ibnd_buff_end   = (ibnd_end_new+1)/2
    max_buff_bands_per_egrp = MAXVAL((all_end(:)+1)/2-(all_start(:)+1)/2)+1
  ELSE
    ibnd_buff_start = ibnd_start_new
    ibnd_buff_end   = ibnd_end_new
    max_buff_bands_per_egrp = MAXVAL(all_end(:)-all_start(:))+1
  ENDIF
  !
  !
  IF( DoLoc ) THEN
    !
    IF (gamma_only) THEN
      IF (.NOT. ALLOCATED(locbuff)) ALLOCATE( locbuff(nrxxs*npol,nbnd,nks) )
      IF (.NOT. ALLOCATED(locmat))  ALLOCATE( locmat(nbnd,nbnd,nks) )
      locbuff = 0.0d0
      locmat = 0.0d0
    ELSE 
      IF (.NOT. ALLOCATED(exxbuff)) ALLOCATE( exxbuff(nrxxs*npol,nbnd,nkqs) )
      IF (.NOT. ALLOCATED(exxmat) ) ALLOCATE( exxmat(nbnd,nkqs,nbnd,nks) )
      exxbuff = (0.0d0, 0.0d0)
      exxmat = 0.0d0
    ENDIF
    !
    IF (.NOT. ALLOCATED(evc0)) then 
      ALLOCATE( evc0(npwx*npol,nbndproj,nks) )
      evc0 = (0.0d0,0.0d0)
    ENDIF
    !
  ELSE
    !
    !ffr: This is the usual case
    !ffr: This is the case of DoLoc = .FALSE.
    !
    IF(.NOT. ALLOCATED(exxbuff)) THEN
      IF (gamma_only) THEN
        ALLOCATE( exxbuff(nrxxs*npol,ibnd_buff_start:ibnd_buff_start + &
                                      max_buff_bands_per_egrp-1,nkqs) ) ! THIS WORKS as for k
      ELSE
        ALLOCATE( exxbuff(nrxxs*npol,ibnd_buff_start:ibnd_buff_start + &
                                      max_buff_bands_per_egrp-1,nkqs) )
      ENDIF
    ENDIF
  ENDIF
  write(*,*) 'nrxxs = ', nrxxs
  write(*,*) 'npol = ', npol
  write(*,*) 'nkqs = ', nkqs
  write(*,*) 'shape exxbuff = ', shape(exxbuff)
  !stop 'Early stop 234'
  !
  !
  ! assign buffer
  IF(DoLoc) THEN
    IF(gamma_only) THEN
      DO ikq=1,SIZE(locbuff,3) 
        DO ibnd=1, x_nbnd_occ 
          DO ir=1,nrxxs*npol
            locbuff(ir,ibnd,ikq)=0.0_DP
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ELSE
    DO ikq = 1, SIZE(exxbuff,3) 
      DO ibnd = ibnd_buff_start, ibnd_buff_end
        DO ir = 1, nrxxs*npol
          exxbuff(ir,ibnd,ikq) = (0.0_DP,0.0_DP)
        ENDDO
      ENDDO
    ENDDO
    ! the above loops will replaced with the following line soon
    !CALL threaded_memset(exxbuff, 0.0_DP, nrxxs*npol*SIZE(exxbuff,2)*nkqs*2)
  ENDIF
  !
  ! This is parallelized over pools. Each pool computes only its k-points
  !
  KPOINTS_LOOP : &
  DO ik = 1, nks
    !
    IF( nks > 1 ) THEN
      CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ik )
    ENDIF
    !ffr: Now evc_exx is wave function for exx
    !ffr: why need to name it *_exx?
    !
    ! ik         = index of k-point in this pool
    ! current_ik = index of k-point over all pools
    !
    current_ik = global_kpoint_index( nkstot, ik )
    !
    IF( gamma_only ) THEN
      !
      IF (MOD(iexx_start,2) == 0) THEN
        ibnd_loop_start = iexx_start-1
      ELSE
        ibnd_loop_start = iexx_start
      ENDIF
      !
      evc_offset = 0
      DO ibnd = ibnd_loop_start, iexx_end, 2
        !
        psic_exx(:) = ( 0._DP, 0._DP )
        !
        IF ( ibnd < iexx_end ) THEN
          IF ( ibnd == ibnd_loop_start .AND. MOD(iexx_start,2) == 0 ) THEN
            DO ig = 1, npwt
              psic_exx(dfftt%nl(ig))  = ( 0._DP, 1._DP )*evc_exx(ig,1)
              psic_exx(dfftt%nlm(ig)) = ( 0._DP, 1._DP )*CONJG(evc_exx(ig,1))
            ENDDO
            evc_offset = -1
          ELSE
            DO ig = 1, npwt
              psic_exx(dfftt%nl(ig))  = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) &
                    + ( 0._DP, 1._DP ) * evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2)
              psic_exx(dfftt%nlm(ig)) = CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) ) &
                    + ( 0._DP, 1._DP ) * CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2) )
            ENDDO
          ENDIF
          !
        ELSE
          DO ig=1,npwt
            psic_exx(dfftt%nl (ig)) = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1)
            psic_exx(dfftt%nlm(ig)) = CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) )
          ENDDO
        ENDIF
        !
        CALL invfft( 'Wave', psic_exx, dfftt )
        !
        IF (DoLoc) THEN
          locbuff(1:nrxxs,ibnd-ibnd_loop_start+evc_offset+1,ik) = DBLE(  psic_exx(1:nrxxs) )
          IF (ibnd-ibnd_loop_start+evc_offset+2 <= nbnd) &
            locbuff(1:nrxxs,ibnd-ibnd_loop_start+evc_offset+2,ik) = AIMAG( psic_exx(1:nrxxs) )
        ELSE
          exxbuff(1:nrxxs,(ibnd+1)/2,current_ik)=psic_exx(1:nrxxs) 
        ENDIF
        !
      ENDDO
      !
    ELSE
      !
      !ffr: This is the usual case
      !ffr: This the case of gamma_only = .FALSE.
      !
      npw = ngk(ik)
      ! ffr: This is loop over all bands.
      !write(*,*) 'iexx_start = ', iexx_start
      !write(*,*) 'iexx_end = ', iexx_end
      !stop 'Stopped here 588'
      IBND_LOOP_K : &
      DO ibnd = iexx_start, iexx_end
        !
        ibnd_exx = ibnd
        IF (noncolin) THEN
          !ffr: noncollinear case
          DO ir = 1, nrxxs
            temppsic_nc(ir,1) = ( 0._DP, 0._DP )
            temppsic_nc(ir,2) = ( 0._DP, 0._DP )
          ENDDO
          DO ig = 1, npw
            temppsic_nc(dfftt%nl(igk_exx(ig,ik)),1) = evc_exx(ig,ibnd-iexx_start+1)
          ENDDO
          CALL invfft( 'Wave', temppsic_nc(:,1), dfftt )
          DO ig = 1, npw
            temppsic_nc(dfftt%nl(igk_exx(ig,ik)),2) = evc_exx(ig+npwx,ibnd-iexx_start+1)
          ENDDO
          CALL invfft( 'Wave', temppsic_nc(:,2), dfftt )
        ELSE
          DO ir = 1, nrxxs
            temppsic(ir) = ( 0._DP, 0._DP )
          ENDDO
          DO ig = 1, npw
            temppsic(dfftt%nl(igk_exx(ig,ik))) = evc_exx(ig,ibnd-iexx_start+1)
          ENDDO
          CALL invfft( 'Wave', temppsic, dfftt )
        ENDIF
        !ffr: dfftt is fft of the wave function, or wave function in real space
        !
        !ffr: At this point temppsic is ready, contains real space repr of evc_exx
        !
        !write(*,*) 'nkqs = ', nkqs
        !stop 'stopped here 623'
        !
        !ffr: This is loop over kpoints
        DO ikq = 1, nkqs
          !
          IF( index_xk(ikq) /= current_ik ) CYCLE
          !
          isym = ABS(index_sym(ikq)) ! ffr: what's this?
          !
          IF (noncolin) THEN ! noncolinear
            !
            DO ir = 1, nxxs
              DO ipol = 1, npol
                psic_nc(ir,ipol) = (0._DP,0._DP)
              ENDDO
            ENDDO
            !
            DO ir = 1, nxxs
              DO ipol = 1, npol
                DO jpol = 1, npol
                  psic_nc(ir,ipol) = psic_nc(ir,ipol) + CONJG(d_spin(jpol,ipol,isym))* &
                                     temppsic_nc(rir(ir,isym),jpol)
                ENDDO
              ENDDO
            ENDDO
            !
            IF (index_sym(ikq) > 0 ) THEN
              ! sym. op. without time reversal: normal case
              DO ir = 1, nrxxs
                exxbuff(ir,ibnd,ikq) = psic_nc(ir,1)
                exxbuff(ir+nrxxs,ibnd,ikq) = psic_nc(ir,2)
              ENDDO
            ELSE
              ! sym. op. with time reversal: spin 1->2*, 2->-1*
              DO ir = 1, nrxxs
                exxbuff(ir,ibnd,ikq) = CONJG(psic_nc(ir,2))
                exxbuff(ir+nrxxs,ibnd,ikq) = -CONJG(psic_nc(ir,1))
              ENDDO
            ENDIF ! index_sym
          !
          ELSE
            !
            !ffr This is the usual one
            !
            DO ir = 1, nrxxs
              psic_exx(ir) = temppsic(rir(ir,isym))
            ENDDO
            !
            DO ir = 1, nrxxs
              IF (index_sym(ikq) < 0 ) THEN
                psic_exx(ir) = CONJG(psic_exx(ir))
              ENDIF
              exxbuff(ir,ibnd,ikq) = psic_exx(ir)
            ENDDO
          ENDIF ! collinear vs noncolinear
          !
        ENDDO ! loop over ikq
      !
      ENDDO IBND_LOOP_K
      !
    ENDIF ! check if using gamma_only or not
    !
  ENDDO KPOINTS_LOOP
  !
  DEALLOCATE( psic_exx )
  IF (noncolin) THEN
    DEALLOCATE( temppsic_nc, psic_nc )
  ELSE IF ( .NOT. gamma_only ) THEN
    DEALLOCATE( temppsic )
  ENDIF
  !
  ! Each wavefunction in exxbuff is computed by a single pool, collect among 
  ! pools in a smart way (i.e. without doing all-to-all sum and bcast)
  ! See also the initialization of working_pool in exx_mp_init
  ! Note that in Gamma-only LSDA can be parallelized over two pools, and there
  ! is no need to communicate anything: each pools deals with its own spin
  !
  IF ( .NOT. gamma_only ) THEN
    DO ikq = 1, nkqs
      CALL mp_bcast( exxbuff(:,:,ikq), working_pool(ikq), intra_orthopool_comm ) 
    ENDDO
  ENDIF
  !
  ! For US/PAW only: compute <beta_I|psi_j,k+q> for the entire 
  ! de-symmetrized k+q grid by rotating the ones from the irreducible wedge
  ! ffr: why?
  IF (okvan) then
    CALL rotate_becxx( nkqs, index_xk, index_sym, xkq_collect )
  ENDIF
  !
  ! Initialize 4-wavefunctions one-center Fock integrals
  !    \int \psi_a(r)\phi_a(r)\phi_b(r')\psi_b(r')/|r-r'|
  !
  IF (okpaw) CALL PAW_init_fock_kernel()
  !
  CALL change_data_structure( .FALSE. )

  ! write out several parameter before exiting
  write(*,*) '----------------------------------------------------------------'
  write(*,*) 'DoLoc = ', DoLoc
  write(*,*) 'negrp = ', negrp
  write(*,*) 'eps_occ = ', eps_occ
  write(*,*) 'use_coulomb_vcut_ws = ', use_coulomb_vcut_ws
  write(*,*) 'use_coulomb_vcut_spheric = ', use_coulomb_vcut_spheric
  write(*,*) 'x_gamma_extrapolation = ', x_gamma_extrapolation
  write(*,*) 'exxdiv_treatment = ', trim(exxdiv_treatment)
  write(*,*) 'use_regularization = ', use_regularization
  write(*,*) 'exxdiv = ', exxdiv
  write(*,*) 'eps_qdiv = ', eps_qdiv
  write(*,*) 'gau_scrlen = ', gau_scrlen
  write(*,*) 'erf_scrlen = ', erf_scrlen
  write(*,*) 'erfc_scrlen = ', erfc_scrlen
  write(*,*) 'yukawa = ', yukawa
  write(*,*) '----------------------------------------------------------------'


  write(*,*)
  write(*,*) '</div> EXIT my_exxinit()'
  write(*,*)

  !stop 'early stop in my_exxinit'

END SUBROUTINE my_exxinit
