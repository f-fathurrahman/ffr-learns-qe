!----------------------------------------------------------------------------
SUBROUTINE my_sum_band()
!----------------------------------------------------------------------------
  !! Calculates the symmetrized charge density and related quantities.  
  !! Also computes the occupations and the sum of occupied eigenvalues.
  !
  USE kinds,                ONLY : DP
  USE ener,                 ONLY : eband
  USE control_flags,        ONLY : diago_full_acc, gamma_only, lxdm, tqr
  USE cell_base,            ONLY : at, bg, omega, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : ngm, g
  USE gvecs,                ONLY : doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k
  USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, is_hubbard_back
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho, rhoz_or_updw
  USE symme,                ONLY : sym_rho
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : nkb, vkb, becsum, ebecsum, nhtol, nhtoj, indv, okvan
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions,        ONLY : evc, psic, psic_nc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE spin_orb,             ONLY : lspinorb, domag, fcoef
  USE wvfct,                ONLY : nbnd, npwx, wg, et, btype
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
  USE mp,                   ONLY : mp_sum
  USE funct,                ONLY : dft_is_meta
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   becp
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ir,   &! counter on 3D r points
             is,   &! counter on spin polarizations
             ig,   &! counter on g vectors
             ibnd, &! counter on bands
             ik,   &! counter on k points
             nt,   &! counter on atomic types
             npol_,&! auxiliary dimension for noncolin case
             ibnd_start, ibnd_end, this_bgrp_nbnd ! first, last and number of band in this bgrp
  REAL (DP), ALLOCATABLE :: kplusg (:)
  !
  !

  write(*,*)
  write(*,*) 'Enter my_sum_band'
  write(*,*)

  !
  becsum(:,:,:) = 0.D0
  if (tqr) ebecsum(:,:,:) = 0.D0
  rho%of_r(:,:)      = 0.D0
  rho%of_g(:,:)      = (0.D0, 0.D0)
  if ( dft_is_meta() .OR. lxdm ) then
     rho%kin_r(:,:)      = 0.D0
     rho%kin_g(:,:)      = (0.D0, 0.D0)
  end if
  
  ! calculates weights of Kohn-Sham orbitals used in calculation of rho
  CALL weights()
  
  ! btype, used in diagonalization, is set here: a band is considered empty
  ! and computed with low accuracy only when its occupation is < 0.01, and
  ! only if option diago_full_acc is false; otherwise, use full accuracy
  btype(:,:) = 1
  IF( .NOT. diago_full_acc ) THEN
    !
    FORALL( ik = 1:nks, wk(ik) > 0.D0 )
       WHERE( wg(:,ik) / wk(ik) < 0.01D0 ) btype(:,ik) = 0
    END FORALL
    !
  ENDIF
  
  !
  ! ... Needed for DFT+U(+V): compute occupations of Hubbard states
  !
  IF (lda_plus_u) THEN
    stop 'Not supported'
  ENDIF
  
  !
  ! ... for band parallelization: set band computed by this processor
  !
  call divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  !
  ! ... Allocate (and later deallocate) arrays needed in specific cases
  !
  IF ( okvan ) CALL allocate_bec_type (nkb, this_bgrp_nbnd, becp, intra_bgrp_comm)
  IF (dft_is_meta() .OR. lxdm) ALLOCATE (kplusg(npwx))



  ! specialized routines are called to sum at Gamma or for each k point 
  ! the contribution of the wavefunctions to the charge
  ! The band energy contribution eband is computed together with the charge
  eband         = 0.D0
  IF(gamma_only) THEN
    stop 'gamma_only in my_sum_band is disabled'
  ELSE
    CALL my_sum_band_k()
  ENDIF

  CALL mp_sum( eband, inter_pool_comm )
  CALL mp_sum( eband, inter_bgrp_comm )
  !
  IF (dft_is_meta() .OR. lxdm) DEALLOCATE (kplusg)
  IF ( okvan ) CALL deallocate_bec_type ( becp )
  !
  ! ... sum charge density over pools (distributed k-points) and bands
  !
  CALL mp_sum( rho%of_r, inter_pool_comm )
  CALL mp_sum( rho%of_r, inter_bgrp_comm )
  IF ( noncolin .AND. .NOT. domag ) rho%of_r(:,2:4)=0.D0
  
  ! bring the unsymmetrized rho(r) to G-space (use psic as work array)
  DO is = 1, nspin
     psic(1:dffts%nnr) = rho%of_r(1:dffts%nnr,is)
     psic(dffts%nnr+1:) = 0.0_dp
     CALL fwfft('Rho', psic, dffts)
     rho%of_g(1:dffts%ngm,is) = psic(dffts%nl(1:dffts%ngm))
     rho%of_g(dffts%ngm+1:,is) = (0.0_dp,0.0_dp)
  END DO

  IF( okvan )  THEN
     !
     ! ... becsum is summed over bands (if bgrp_parallelization is done)
     ! ... and over k-points (but it is not symmetrized)
     !
     CALL mp_sum(becsum, inter_bgrp_comm )
     CALL mp_sum(becsum, inter_pool_comm )
     !
     ! ... same for ebecsum, a correction to becsum (?) in real space
     !
     IF (tqr) CALL mp_sum(ebecsum, inter_pool_comm )
     IF (tqr) CALL mp_sum(ebecsum, inter_bgrp_comm )
     !
     ! ... PAW: symmetrize becsum and store it
     ! ... FIXME: the same should be done for USPP as well (ffr: we need to check this!)
     !
     IF ( okpaw ) THEN
        rho%bec(:,:,:) = becsum(:,:,:)
        CALL PAW_symmetrize(rho%bec)
     END IF
     !
     ! ... Here we add the (unsymmetrized) Ultrasoft contribution to the charge
     !
     CALL my_addusdens( rho%of_g(:,:) )
  ENDIF

  ! symmetrize rho(G) 
  CALL sym_rho( nspin_mag, rho%of_g )

  ! synchronize rho%of_r to the calculated rho%of_g (use psic as work array)
  DO is = 1, nspin_mag
    psic(:) = ( 0.D0, 0.D0 )
    psic(dfftp%nl(:)) = rho%of_g(:,is)
    IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( rho%of_g(:,is) )
    CALL invfft('Rho', psic, dfftp)
    rho%of_r(:,is) = psic(:)
  END DO
  
  ! rho_kin(r): sum over bands, k-points, bring to G-space, symmetrize,
  ! synchronize with rho_kin(G)
  IF ( dft_is_meta() .OR. lxdm) THEN
     !
     CALL mp_sum( rho%kin_r, inter_pool_comm )
     CALL mp_sum( rho%kin_r, inter_bgrp_comm )
     DO is = 1, nspin
        psic(1:dffts%nnr) = rho%kin_r(1:dffts%nnr,is)
        psic(dffts%nnr+1:) = 0.0_dp
        CALL fwfft ('Rho', psic, dffts)
        rho%kin_g(1:dffts%ngm,is) = psic(dffts%nl(1:dffts%ngm))
     END DO
     !
     IF (.NOT. gamma_only) CALL sym_rho( nspin, rho%kin_g )
     !
     DO is = 1, nspin
        psic(:) = ( 0.D0, 0.D0 )
        psic(dfftp%nl(:)) = rho%kin_g(:,is)
        IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( rho%kin_g(:,is) )
        CALL invfft ('Rho', psic, dfftp)
        rho%kin_r(:,is) = psic(:)
     END DO
     !
  END IF
  !
  ! ... if LSDA rho%of_r and rho%of_g are converted from (up,dw) to
  ! ... (up+dw,up-dw) format.
  !
  IF ( nspin == 2 ) CALL rhoz_or_updw( rho, 'r_and_g', '->rhoz' )
  !
  RETURN


CONTAINS

! Inner subroutines


!-----------------------------------------------------------------------
SUBROUTINE my_sum_band_k()
!-----------------------------------------------------------------------
  !! \(\texttt{sum_band}\) - part for k-points version
  !
  USE mp_bands,     ONLY : me_bgrp
  USE mp,           ONLY : mp_sum, mp_get_comm_null
  USE fft_helper_subroutines, ONLY : fftx_ntgrp, fftx_tgpe, &
                     tg_reduce_rho, tg_get_nnr, tg_get_group_nr3
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  REAL(DP) :: w1
  ! weights
  INTEGER :: npw, ipol, na, np
  !
  INTEGER  :: idx, ioff, ioff_tg, nxyp, incr, v_siz, j, ir3
  COMPLEX(DP), ALLOCATABLE :: tg_psi(:), tg_psi_nc(:,:)
  REAL(DP),    ALLOCATABLE :: tg_rho(:), tg_rho_nc(:,:)
  LOGICAL  :: use_tg
  INTEGER :: right_nnr, right_nr3, right_inc, ntgrp
  !
  ! chunking parameters
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock

  ! here we sum for each k point the contribution
  ! of the wavefunctions to the charge
  use_tg = ( dffts%has_task_groups ) .AND. ( .NOT. (dft_is_meta() .OR. lxdm) )

  incr = 1

  IF( use_tg ) THEN
    stop 'use_tg = .true. is disabled'
  ENDIF


  k_loop: DO ik = 1, nks

    IF( lsda ) current_spin = isk(ik)
    npw = ngk(ik)

    IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )

    IF ( nkb > 0 ) CALL my_init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )

    ! here we compute the band energy: the sum of the eigenvalues
    DO ibnd = ibnd_start, ibnd_end, incr

      eband = eband + et( ibnd, ik ) * wg( ibnd, ik )

      ! the sum of eband and demet is the integral for e < ef of
      ! e n(e) which reduces for degauss=0 to the sum of the
      ! eigenvalues
      w1 = wg(ibnd,ik) / omega

      CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
      !
      DO j = 1, npw
         psic(dffts%nl(igk_k(j,ik))) = evc(j,ibnd)
      ENDDO

      CALL invfft('Wave', psic, dffts) ! to real space

      ! increment the charge density
      ! The operation is done on smooth grid 
      CALL get_rho( rho%of_r(:,current_spin), dffts%nnr, w1, psic )

      IF( dft_is_meta() .OR. lxdm) THEN
        DO j=1,3
          psic(:) = ( 0.D0, 0.D0 )
          !
          kplusg (1:npw) = (xk(j,ik)+g(j,igk_k(1:npw,ik))) * tpiba
          psic(dffts%nl(igk_k(1:npw,ik)))=CMPLX(0d0,kplusg(1:npw),kind=DP) * &
                                  evc(1:npw,ibnd)
          !
          CALL invfft ('Wave', psic, dffts)
          !
          ! ... increment the kinetic energy density ...
          !
          CALL get_rho(rho%kin_r(:,current_spin), dffts%nnr, w1, psic)
        ENDDO
      ENDIF ! dft_is_meta
    
    ENDDO

    ! If we have a US pseudopotential we compute here the becsum term
    IF( okvan ) CALL my_sum_bec( ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd ) 
          !
  END DO k_loop
  
  RETURN

END SUBROUTINE
     

! Inner subroutine
!-------------------------------------------------------
SUBROUTINE get_rho(rho_loc, nrxxs_loc, w1_loc, psic_loc)
!-------------------------------------------------------
  IMPLICIT NONE

  INTEGER :: nrxxs_loc
  REAL(DP) :: rho_loc(nrxxs_loc)
  REAL(DP) :: w1_loc
  COMPLEX(DP) :: psic_loc(nrxxs_loc)
  INTEGER :: ir

  DO ir = 1, nrxxs_loc
    rho_loc(ir) = rho_loc(ir) + w1_loc * (DBLE(psic_loc(ir))**2 + AIMAG( psic_loc(ir) )**2)
  ENDDO

END SUBROUTINE get_rho


END SUBROUTINE

