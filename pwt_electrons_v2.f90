! ffr: 18 Dec 2015
!
! A stripped version of pwscf
!
! Aimed to debug and 'refactor' self-consistent calculation
! Test the CheFSI implementation


! A driver for CheFSI
SUBROUTINE c_bands_cheby()
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, stdout
  USE io_files, ONLY : iunigk, iunwfc, nwordwfc
  USE buffers, ONLY : get_buffer, save_buffer
  USE klist, ONLY : nks, ngk, xk
  USE wvfct, ONLY : npw, current_k, igk, nbnd, et
  USE lsda_mod, ONLY : isk, lsda, current_spin
  USE wavefunctions_module, ONLY : evc
  USE uspp, ONLY : nkb, vkb
  USE bp, ONLY : lelfield  ! TODO: we don't really want to concentrate on this for the moment
  USE becmod, ONLY : becp, allocate_bec_type, deallocate_bec_type
  USE mp_bands, ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik_, ik
  INTEGER :: nlancz
  REAL(8) :: lb, ub
  ! ik : counter on k points
  ! ik_: k-point already done in a previous run

  CALL start_clock('c_bands_cheby')

  IF( ionode ) THEN
    WRITE(stdout,*)
    WRITE(stdout,*) '>>>>> Calling c_bands_cheby'
  ENDIF

  ik_ = 0
  ! restart case is skipped

  IF(nks > 1) REWIND( iunigk )
  
  ! ChebyFSI step for each k
  k_loop: DO ik = 1, nks
    !
    current_k = ik
    IF(lsda) current_spin = isk(ik)
    npw = ngk(ik)
    !
    ! Reads the lists of indices k+G <-> G of this k point
    IF( nks > 1) READ( iunigk ) igk

    ! Dirty restart trick is skipped
    
    ! Ultrasoft PP
    IF( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )

    ! Kinetic energy
    CALL g2_kin( ik )
    
    ! Read wavefun from previous iteration
    IF( nks > 1 .OR. lelfield ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )

    ! LDA+U is skipped 
    
    ! Required for h_psi
    ! FIXME: Is this required for each k?
    CALL allocate_bec_type( nkb, nbnd, becp, intra_bgrp_comm )

    !IF(ionode) WRITE(stdout,*) 'Before cheFSI'
    !CALL test_orthonorm_evc()

    ! call cheFSI here
    nlancz = min( nbnd, 10 )
    !CALL ef_lanczos( ik, npw, nlancz, lb, ub )
    CALL kstep_lanczos( ik, npw, nlancz, lb, ub )
    !
    CALL chefsi( ik, 9, lb, ub )
    !CALL chefsi_scaled( ik, 9, lb, ub, et(1,ik) )
    
    ! Save wavefun to be used as input for the next SCF iteration and for 
    ! electron density calculation
    IF( nks > 1 .OR. lelfield ) CALL save_buffer( evc, nwordwfc, iunwfc, ik )

    ! skipped check_stop_now() and save_in_cbands()

    ! FIXME: Is this required for each k?
    CALL deallocate_bec_type( becp )
  ENDDO k_loop

  ! avg_iter is skipped. It is probably won't be used anyway

  CALL stop_clock('c_bands_cheby')

END SUBROUTINE



!------------------------------------------------------------------------------
SUBROUTINE dbg_electrons_cheby
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, stdout
  USE io_files, ONLY : iunmix, output_drho
  USE cell_base, ONLY : alat, at, bg, omega
  USE ions_base, ONLY : nat, nsp, ityp, zv, tau
  USE gvect, ONLY : g, gg, gcutm, gstart, ngm
  USE gvecs, ONLY : doublegrid
  USE ener, ONLY : ewld, hwf_energy, eband, etxc, etxcc, ehart, demet, &
     deband, vtxc, etot
  USE klist, ONLY : nelec, nkstot, nks, lgauss
  USE wvfct, ONLY : et, nbnd
  USE fft_base, ONLY : dfftp
  USE dfunct, ONLY : newd
  !
  USE mp, ONLY : mp_sum, mp_bcast
  USE mp_bands, ONLY : intra_bgrp_comm
  USE mp_pools, ONLY : root_pool, inter_pool_comm
  !
  USE control_flags, ONLY : gamma_only, istep, ethr, niter, mixing_beta, nmix, &
    conv_elec, scf_must_converge, iprint, io_level
  USE ldaU, ONLY : lda_plus_u, eth
  USE noncollin_module, ONLY : noncolin
  USE lsda_mod, ONLY : nspin
  USE paw_variables, ONLY : okpaw, ddd_paw
  USE extfield, ONLY : etotefield
  USE io_rho_xml, ONLY : write_rho
  !
  USE vlocal, ONLY : strf
  !
  USE scf, ONLY : scf_type, create_scf_type, open_mix_file, scf_type_COPY, &
    bcast_scf_type, kedtau, close_mix_file, destroy_scf_type
  !
  USE scf, ONLY : rho, v, rho_core, rhog_core, vnew, vltot, vrs
  !
  IMPLICIT NONE
  !
  INTEGER :: printout
  LOGICAL :: exst, first
  INTEGER :: iter, idum, ik, ib
  REAL(DP) :: dr2, tr2_min
  REAL(DP) :: deband_hwf, descf
  REAL(DP) :: charge
  !
  TYPE(scf_type) :: rhoin  ! used to store rho_in of current/next iteration
  ! external functionst
  REAL(DP), EXTERNAL :: ewald, get_clock

  CALL start_clock('dbg_electrons_cheby')

  IF(ionode) THEN
    WRITE(stdout,*)
    WRITE(stdout,*) 'Calling dbg_electrons_cheby'
    WRITE(stdout,*) '***************************'
  ENDIF

  printout = 2

  iter = 0  ! the actual SCF iteration counter

  dr2 = 0.0_DP  ! the norm of the diffence between potential

  IF(ionode) WRITE(stdout,*) 'istep, ethr = ', istep, ethr
  ! Convergence threshold for iterative diagonalization
  ! for the first scf iteration of each ionic step (after the first),
  ! the threshold is fixed to a default value of 1.D-6
  IF ( istep > 0 ) ethr = 1.D-6

  ! Calculation of Ewald energy
  ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
  IF(ionode) WRITE(stdout,*) 'Ewald energy:', ewld


  CALL create_scf_type( rhoin )

  CALL open_mix_file( iunmix, 'mix', exst )
  ! it seems that the argument exst is not used anymore after this CALL

  niter = 60  ! FIXME overwritten for debugging purpose, use from input instead

  IF(ionode) THEN
    WRITE(stdout,*)
    WRITE(stdout,*) 'Starting eigenvalues (before SCF):'
    DO ik=1,nkstot
      WRITE(stdout,*) 'kpoint #', ik
      DO ib=1,nbnd
        WRITE(stdout,'(1x,I5,F18.10)') ib, et(ib,ik)
      ENDDO
    ENDDO
  ENDIF

  DO idum=1,niter

    ! idum is dummy SCF iteration counter, iter is the actual counter
    iter = iter + 1
    IF(ionode) THEN
      WRITE(stdout,*)
      WRITE(stdout,*) '======================='
      WRITE(stdout,'(1x,A,I5,A)') '=== ITERATION:', iter, ' ==='
      WRITE(stdout,*) '======================='
    ENDIF

    ! Probably this is not needed for ChebyFSI
    ! Convergence threshold for iterative diagonalization is
    ! automatically updated during self consistency
    IF ( iter > 1 ) THEN
      IF ( iter == 2 ) ethr = 1.D-2
      ethr = MIN( ethr, 0.1D0*dr2 / MAX( 1.D0, nelec ) )
      ! do not allow convergence threshold to become too small:
      ! iterative diagonalization may become unstable
      ethr = MAX( ethr, 1.D-13 )
    END IF

    first = ( iter == 1 )
    
    ! deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v> is calculated a
    ! first time here using the input density and potential ( to be
    ! used to calculate the Harris-Weinert-Foulkes energy )
    deband_hwf = delta_e()

    ! save input density to rhoin
    CALL scf_type_COPY( rho, rhoin )

    scf_step: DO

      ! tr2_min is set to an estimate of the error on the energy
      ! due to diagonalization - used only for the first scf iteration
      tr2_min = 0.0_DP
      IF(first) tr2_min = ethr*max(1.0_DP, nelec)

      ! iterative diagonalization of KS Hamiltonian
      IF(first) THEN
        CALL c_bands( iter )
      ELSE
        CALL c_bands_cheby()
      ENDIF

      ! xk, wk, isk, et, wg are distributed across pools;
      ! the first node has a complete copy of xk, wk, isk,
      ! while eigenvalues et and weights wg must be
      ! explicitely collected to the first node
      ! this is done here for et, in sum_band for wg
      CALL poolrecover( et, nbnd, nkstot, nks )
      !IF(ionode) THEN
      !  DO ik=1,nkstot
      !    WRITE(stdout,*) 'kpoint #', ik
      !    DO ib=1,nbnd
      !      WRITE(stdout,'(1x,I5,F18.10)') ib, et(ib,ik)
      !    ENDDO
      !  ENDDO
      !ENDIF

      ! the new density is computed here. For PAW:
      ! sum_band computes new becsum (stored in uspp modules)
      ! and a subtly different copy in rho%bec (scf module)
      CALL sum_band()

      ! the Harris-Weinert-Foulkes energy is computed here using only
      ! quantities obtained from the input density
      hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet
      ! okpaw and lda_plus_u cases are skipped
      !IF(ionode) THEN
      !  WRITE(stdout,*) 'hwf_energy =', hwf_energy
      !  WRITE(stdout,*) 'eband      =', eband
      !  WRITE(stdout,*) 'deband_hwf =', deband_hwf
      !  WRITE(stdout,*) 'etxc       =', etxc
      !  WRITE(stdout,*) 'etxcc      =', etxcc
      !  WRITE(stdout,*) 'ewld       =', ewld
      !  WRITE(stdout,*) 'ehart      =', ehart
      !  WRITE(stdout,*) 'demet      =', demet
      !ENDIF

      ! lsda or noncolin case is skipped

      ! eband  = \sum_v \epsilon_v    is calculated by sum_band
      ! deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v>
      ! eband + deband = \sum_v <\psi_v | T + Vion |\psi_v>
      deband = delta_e()

      ! mix_rho mixes several quantities: rho in g-space, tauk (for
      ! meta-gga), ns and ns_nc (for lda+u) and becsum (for paw)
      ! Results are broadcast from pool 0 to others to prevent trouble
      ! on machines unable to yield the same results from the same 
      ! calculation on same data, performed on different procs
      ! The mixing should be done on pool 0 only as well, but inside
      ! mix_rho there is a call to rho_ddot that in the PAW CASE 
      ! contains a hidden parallelization level on the entire image
      !
      ! IF ( my_pool_id == root_pool ) 
      CALL mix_rho ( rho, rhoin, mixing_beta, dr2, tr2_min, iter, nmix, &
                     iunmix, conv_elec )
      CALL bcast_scf_type ( rhoin, root_pool, inter_pool_comm )
      CALL mp_bcast ( dr2, root_pool, inter_pool_comm )
      CALL mp_bcast ( conv_elec, root_pool, inter_pool_comm )
      !
      IF (.not. scf_must_converge .and. idum == niter) conv_elec = .true.
      IF(ionode) WRITE(stdout,*) 'conv_elec =', conv_elec

      ! if convergence is achieved or if the self-consistency error
      ! (dr2) is smaller than the estimated error due to diagonalization
      ! (tr2_min), rhoin and rho are unchanged: rhoin contains the input
      !  density and rho contains the output density
      ! In the other cases rhoin contains the mixed charge density 
      ! (the new input density) while rho is unchanged
      IF ( first .and. nat > 0) THEN
         ! first scf iteration: check if the threshold on diagonalization
         ! (ethr) was small enough wrt the error in self-consistency (dr2)
         ! if not, perform a new diagonalization with reduced threshold
         first = .FALSE.
         !
         IF ( dr2 < tr2_min ) THEN
            !
            WRITE( stdout, '(/,5X,"Threshold (ethr) on eigenvalues was ", &
                             &    "too large:",/,5X,                      &
                             & "Diagonalizing with lowered threshold",/)' )
            !
            ethr = 0.1D0*dr2 / MAX( 1.D0, nelec )
            !
            CYCLE scf_step
            !
         END IF
         !
      END IF
      !
      IF ( .NOT. conv_elec ) THEN
         ! no convergence yet: calculate new potential from mixed
         ! charge density (i.e. the new estimate)
         CALL v_of_rho( rhoin, rho_core, rhog_core, &
                        ehart, etxc, vtxc, eth, etotefield, charge, v)
         !
         ! estimate correction needed to have variational energy:
         ! T + E_ion (eband + deband) are calculated in sum_band
         ! and delta_e using the output charge density rho;
         ! E_H (ehart) and E_xc (etxc) are calculated in v_of_rho
         ! above, using the mixed charge density rhoin%of_r.
         ! delta_escf corrects for this difference at first order
         descf = delta_escf()
         !
         ! now copy the mixed charge density in R- and G-space in rho
         !
         CALL scf_type_COPY( rhoin, rho )
         !
      ELSE 
         !
         ! convergence reached:
         ! 1) the output HXC-potential is saved in v
         ! 2) vnew contains V(out)-V(in) ( used to correct the forces ).
         !
         vnew%of_r(:,:) = v%of_r(:,:)
         CALL v_of_rho( rho,rho_core,rhog_core, &
                        ehart, etxc, vtxc, eth, etotefield, charge, v)
         vnew%of_r(:,:) = v%of_r(:,:) - vnew%of_r(:,:)
         !
         ! note that rho is here the output, not mixed, charge density
         ! so correction for variational energy is no longer needed
         !
         descf = 0._dp
         !
      END IF

      EXIT scf_step
    ENDDO scf_step

    !
    ! ... define the total local potential (external + scf)
    !
    CALL sum_vrs( dfftp%nnr, nspin, vltot, v%of_r, vrs )
    !
    ! ... interpolate the total local potential
    !
    CALL interpolate_vrs( dfftp%nnr, nspin, doublegrid, kedtau, v%kin_r, vrs )
    !
    ! ... in the US case we have to recompute the self-consistent
    ! ... term in the nonlocal potential
    ! ... PAW: newd contains PAW updates of NL coefficients
    !
    CALL newd()
    !
    IF ( conv_elec ) WRITE( stdout, 9101 )
    !
    IF ( conv_elec .OR. MOD( iter, iprint ) == 0 ) THEN
       !
       CALL print_ks_energies()
       !
    END IF
    !
    IF ( ABS( charge - nelec ) / charge > 1.D-7 ) THEN
       WRITE( stdout, 9050 ) charge, nelec
       IF ( ABS( charge - nelec ) / charge > 1.D-3 ) THEN
          IF (.not.lgauss) THEN
             CALL errore( 'electrons', 'charge is wrong: smearing is needed', 1 )
          ELSE
             CALL errore( 'electrons', 'charge is wrong', 1 )
          END IF
       END IF
    END IF
    !
    etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet + descf
    !
    CALL print_energies ( printout )
    !
    IF ( conv_elec ) THEN
       !
       WRITE( stdout, 9110 ) iter
       !
       ! ... jump to the end
       !
       GO TO 10
       !
    END IF
    !
    ! ... uncomment the following line if you wish to monitor the evolution
    ! ... of the force calculation during self-consistency
    !
    !CALL forces()
    !
    ! ... it can be very useful to track internal clocks during
    ! ... self-consistency for benchmarking purposes
  END DO
  !
  WRITE( stdout, 9101 )
  WRITE( stdout, 9120 ) iter
  !
10  FLUSH( stdout )
  !
  ! ... exiting: write (unless disables) the charge density to file
  ! ... (also write ldaU ns coefficients and PAW becsum)
  !
  IF ( io_level > -1 ) CALL write_rho( rho, nspin )
  !
  ! ... delete mixing info if converged, keep it if not
  !
  IF ( conv_elec ) THEN
     CALL close_mix_file( iunmix, 'delete' )
  ELSE
     CALL close_mix_file( iunmix, 'keep' )
  END IF
  !
  IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
  !
  CALL destroy_scf_type ( rhoin )
  !
  CALL stop_clock( 'dbg_electrons_cheby' )
  !
  RETURN
  !
  ! ... formats
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9001 FORMAT(/'     per-process dynamical memory: ',f7.1,' Mb' )
9002 FORMAT(/'     Self-consistent Calculation' )
9010 FORMAT(/'     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F4.2 )
9050 FORMAT(/'     WARNING: integrated charge=',F15.8,', expected=',F15.8 )
9101 FORMAT(/'     End of self-consistent calculation' )
9110 FORMAT(/'     convergence has been achieved in ',i3,' iterations' )
9120 FORMAT(/'     convergence NOT achieved after ',i3,' iterations: stopping' )



! inner functions and subroutines ---------------------------------------------

CONTAINS

!TODO : simplify these subroutines and functions


     !-----------------------------------------------------------------------
     FUNCTION delta_e()
     !-----------------------------------------------------------------------
       ! ... delta_e = - \int rho%of_r(r)  v%of_r(r)
       !               - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !               - \sum rho%ns       v%ns       [for LDA+U]
       !               - \sum becsum       D1_Hxc     [for PAW]
       USE funct,  ONLY : dft_is_meta
       IMPLICIT NONE
       REAL(DP) :: delta_e, delta_e_hub
       !
       delta_e = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
       !
       IF ( dft_is_meta() ) &
          delta_e = delta_e - SUM( rho%kin_r(:,:)*v%kin_r(:,:) )
       !
       delta_e = omega * delta_e / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
       CALL mp_sum( delta_e, intra_bgrp_comm )
       !
       if (lda_plus_u) then
         if (noncolin) then
           delta_e_hub = - SUM (rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:))
           delta_e = delta_e + delta_e_hub
         else
           delta_e_hub = - SUM (rho%ns(:,:,:,:)*v%ns(:,:,:,:))
           if (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
           delta_e = delta_e + delta_e_hub
         endif
       end if
       !
       IF (okpaw) delta_e = delta_e - SUM(ddd_paw(:,:,:)*rho%bec(:,:,:))
       !
       RETURN
       !
     END FUNCTION delta_e


     !-----------------------------------------------------------------------
     FUNCTION delta_escf()
       !-----------------------------------------------------------------------
       !
       ! ... delta_escf = - \int \delta rho%of_r(r)  v%of_r(r)
       !                  - \int \delta rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !                  - \sum \delta rho%ns       v%ns       [for LDA+U]
       !                  - \sum \delta becsum       D1         [for PAW]
       ! ... calculates the difference between the Hartree and XC energy
       ! ... at first order in the charge density difference \delta rho(r)
       !
       USE funct,  ONLY : dft_is_meta
       IMPLICIT NONE
       REAL(DP) :: delta_escf, delta_escf_hub
       !
       delta_escf = - SUM( ( rhoin%of_r(:,:)-rho%of_r(:,:) )*v%of_r(:,:) )
       !
       IF ( dft_is_meta() ) &
          delta_escf = delta_escf - &
                       SUM( (rhoin%kin_r(:,:)-rho%kin_r(:,:) )*v%kin_r(:,:))
       !
       delta_escf = omega * delta_escf / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
       CALL mp_sum( delta_escf, intra_bgrp_comm )
       !
       if (lda_plus_u) then
         if (noncolin) then
           delta_escf_hub = - SUM((rhoin%ns_nc(:,:,:,:)-rho%ns_nc(:,:,:,:))*v%ns_nc(:,:,:,:))
           delta_escf = delta_escf + delta_escf_hub
         else
           delta_escf_hub = - SUM((rhoin%ns(:,:,:,:)-rho%ns(:,:,:,:))*v%ns(:,:,:,:))
           if (nspin==1) delta_escf_hub = 2.d0 * delta_escf_hub
           delta_escf = delta_escf + delta_escf_hub
         endif
       end IF

       IF (okpaw) delta_escf = delta_escf - &
                               SUM(ddd_paw(:,:,:)*(rhoin%bec(:,:,:)-rho%bec(:,:,:)))

       RETURN
       !
     END FUNCTION delta_escf


     !-----------------------------------------------------------------------
     SUBROUTINE print_energies ( printout )
       !-----------------------------------------------------------------------
       !
       USE constants, ONLY : eps8
       INTEGER, INTENT (IN) :: printout
       !
       IF ( printout == 0 ) RETURN
       IF ( ( conv_elec .OR. MOD(iter,iprint) == 0 ) .AND. printout > 1 ) THEN
          !
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9081 ) etot, hwf_energy, dr2
          ELSE
             WRITE( stdout, 9083 ) etot, hwf_energy, dr2
          END IF
          !
          WRITE( stdout, 9060 ) &
               ( eband + deband ), ehart, ( etxc - etxcc ), ewld
          !
          IF ( ABS (descf) > eps8 ) WRITE( stdout, 9069 ) descf
          !
          ! ... With Fermi-Dirac population factor, etot is the electronic
          ! ... free energy F = E - TS , demet is the -TS contribution
          !
          IF ( lgauss ) WRITE( stdout, 9070 ) demet
          !
       ELSE IF ( conv_elec ) THEN
          !
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9081 ) etot, hwf_energy, dr2
          ELSE
             WRITE( stdout, 9083 ) etot, hwf_energy, dr2
          END IF
          !
       ELSE
          !
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9080 ) etot, hwf_energy, dr2
          ELSE
             WRITE( stdout, 9082 ) etot, hwf_energy, dr2
          END IF
       END IF
       !
       FLUSH( stdout )
       !
       RETURN
       !
9060 FORMAT(/'     The total energy is the sum of the following terms:',/,&
            /'     one-electron contribution =',F17.8,' Ry' &
            /'     hartree contribution      =',F17.8,' Ry' &
            /'     xc contribution           =',F17.8,' Ry' &
            /'     ewald contribution        =',F17.8,' Ry' )
9069 FORMAT( '     scf correction            =',F17.8,' Ry' )
9070 FORMAT( '     smearing contrib. (-TS)   =',F17.8,' Ry' )
9077 FORMAT( '     External forces energy    =',F17.8,' Ry' )
9080 FORMAT(/'     total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',0PF17.8,' Ry' )
9081 FORMAT(/'!    total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',0PF17.8,' Ry' )
9082 FORMAT(/'     total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',1PE17.1,' Ry' )
9083 FORMAT(/'!    total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',1PE17.1,' Ry' )
9085 FORMAT(/'     total all-electron energy =',0PF17.6,' Ry' )

  END SUBROUTINE print_energies



END SUBROUTINE



!------------------------------------------------------------------------------
PROGRAM t_pwscf
!------------------------------------------------------------------------------
  USE environment, ONLY: environment_start
  USE mp_global, ONLY: mp_startup
  USE read_input, ONLY: read_input_file
  USE command_line_options, ONLY: input_file_
  !
  USE check_stop, ONLY: check_stop_init
  IMPLICIT NONE

  CALL mp_startup()
  CALL environment_start( 'PWSCF' )
  CALL read_input_file( 'PW', input_file_ )
  
  ! Details of run_pwscf
  CALL iosys()
  CALL check_stop_init()
  CALL setup()
  CALL init_run()

  CALL dbg_electrons_cheby()

  CALL stop_run( 0 )
  CALL do_stop( 0 )
  STOP
END PROGRAM

