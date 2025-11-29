! remove clock related calls
! remove lda_plus_u stuffs


!----------------------------------------------------------------------------
SUBROUTINE my_electrons_scf( printout, exxen )
  !----------------------------------------------------------------------------
  !! This routine is a driver of the self-consistent cycle.  
  !! It uses the routine \(\texttt{c_bands}\) for computing the bands at fixed
  !! Hamiltonian, the routine \(\texttt{sum_band}\) to compute the charge density,
  !! the routine \(\texttt{v_of_rho}\) to compute the new potential and the routine
  !! \(\text{mix_rho}\) to mix input and output charge densities.
  !
  USE kinds,                ONLY : DP
  USE check_stop,           ONLY : check_stop_now, stopped_by_user
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : at, bg, alat, omega
  USE ions_base,            ONLY : zv, nat, nsp, ityp, tau, compute_eextfor, atm, &
                                   ntyp => nsp
  USE bp,                   ONLY : lelfield
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, gstart, g, gg, gcutm
  USE gvecs,                ONLY : doublegrid
  USE klist,                ONLY : nelec, nks, nkstot, lgauss, &
                                   tot_charge
  USE fixed_occ,            ONLY : one_atom_occupations
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag
  USE vlocal,               ONLY : strf
  USE wvfct,                ONLY : nbnd, et
  USE gvecw,                ONLY : ecutwfc
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon, edftd3, exdm, ef
  USE scf,                  ONLY : scf_type, scf_type_COPY, bcast_scf_type,&
                                   create_scf_type, destroy_scf_type, &
                                   open_mix_file, close_mix_file, &
                                   rho, rho_core, rhog_core, v, vltot, vrs, &
                                   kedtau, vnew
  USE control_flags,        ONLY : mixing_beta, ethr, niter, nmix, &
                                   iprint, conv_elec, &
                                   restart, io_level, do_makov_payne,  &
                                   gamma_only, iverbosity, textfor,     &
                                   llondon, ldftd3, scf_must_converge, lxdm, ts_vdw
  USE control_flags,        ONLY : n_scf_steps, scf_error

  USE io_files,             ONLY : iunmix, output_drho
  USE ldaU,                 ONLY : eth
  USE extfield,             ONLY : tefield, etotefield, gate, etotgatefield !TB
  USE noncollin_module,     ONLY : noncolin, magtot_nc, i_cons,  bfield, &
                                   lambda, report
  USE spin_orb,             ONLY : domag
  USE io_rho_xml,           ONLY : write_scf
  USE mp_pools,             ONLY : root_pool, &
                                   inter_pool_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  !
  USE london_module,        ONLY : energy_london
  USE dftd3_api,            ONLY : dftd3_pbc_dispersion, &
                                   dftd3_init, dftd3_set_functional, &
                                   get_atomic_number, dftd3_input, &
                                   dftd3_calc
  USE dftd3_qe,             ONLY : dftd3, energy_dftd3
  USE xdm_module,           ONLY : energy_xdm
  USE tsvdw_module,         ONLY : EtsvdW
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  !USE my_dfunct,            ONLY : my_newd
  USE dfunct,            ONLY : newd
  USE esm,                  ONLY : do_comp_esm, esm_printpot, esm_ewald
  USE fcp_variables,        ONLY : lfcpopt, lfcpdyn
  USE wrappers,             ONLY : memstat
  USE iso_c_binding,        ONLY : c_int
  !
  USE plugin_variables,     ONLY : plugin_etot
  !
  IMPLICIT NONE

  INTEGER, INTENT (IN) :: printout
  !! * If printout>0, prints on output the total energy;
  !! * if printout>1, also prints decomposition into energy contributions.

  REAL(DP),INTENT (IN) :: exxen !! current estimate of the exchange energy

  ! local variables
  REAL(DP) :: dr2 !! the norm of the diffence between potential
  REAL(DP) :: charge !! the total charge
  REAL(DP) :: deband_hwf !! deband for the Harris-Weinert-Foulkes functional
  INTEGER :: idum !! dummy counter on iterations
  INTEGER :: iter !! counter on iterations

  INTEGER :: kilobytes
  !
  REAL(DP) :: tr2_min !! estimated error on energy coming from diagonalization
  REAL(DP) :: descf !! correction for variational energy
  REAL(DP) :: en_el=0.0_DP !! electric field contribution to the total energy
  
  REAL(DP) :: eext=0.0_DP !! external forces contribution to the total energy

  !! auxiliary variables for calculating and storing temporary copies of
  !! the charge density and of the HXC-potential
  LOGICAL :: first, exst

  !! used to store rho_in of current/next iteration
  TYPE(scf_type) :: rhoin

  !
  ! external functions
  !
  REAL(DP), EXTERNAL :: ewald, my_delta_e, my_delta_escf, my_calc_pol
  REAL(DP) :: etot_cmp_paw(nat,2,2)

  !! auxiliary variables for grimme-d3  
  REAL(DP) :: latvecs(3,3)
  INTEGER:: atnum(1:nat), na

  !! if .TRUE. then background states are present (DFT+U)
  LOGICAL :: lhb 
  

  write(*,*)
  write(*,*) '<div> ENTER my_electrons_scf'
  write(*,*)

  lhb = .FALSE.
  iter = 0
  dr2  = 0.0_dp
  IF( restart ) CALL restart_in_electrons( iter, dr2, ethr, et )

  !
  CALL memstat( kilobytes )
  IF( kilobytes > 0 ) WRITE( stdout, 9001 ) kilobytes/1000.0
  !
  FLUSH( stdout )

  ! calculates the ewald contribution to total energy
  IF ( do_comp_esm ) THEN
    ewld = esm_ewald()
  ELSE
    ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                  omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
  ENDIF

  IF ( llondon ) THEN
    elondon = energy_london( alat , nat , ityp , at ,bg , tau )
  ELSE
    elondon = 0.d0
  ENDIF
  
  !
  ! Grimme-D3 correction to the energy
  !
  IF(ldftd3) THEN
    latvecs(:,:) = at(:,:)*alat
    tau(:,:) = tau(:,:)*alat
    DO na = 1, nat
       atnum(na) = get_atomic_number(TRIM(atm(ityp(na))))
    ENDDO
    call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,energy_dftd3)
    edftd3 = energy_dftd3*2.d0 ! to Ry
    tau(:,:) = tau(:,:)/alat
  ELSE
    edftd3 =  0.0
  ENDIF


  CALL create_scf_type( rhoin )
  
  WRITE( stdout, 9002 )
  9002 FORMAT(/'     Self-consistent Calculation' )
  FLUSH( stdout )

  CALL open_mix_file( iunmix, 'mix', exst )

  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%          iterate !          %%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  DO idum = 1, niter

    !write(*,*) 'idum = ', idum
    !
    IF ( check_stop_now() ) THEN
      conv_elec=.FALSE.
      CALL save_in_electrons (iter, dr2, ethr, et )
      GOTO 10
    ENDIF
    
    iter = iter + 1

    WRITE( stdout, 9010 ) iter, ecutwfc, mixing_beta
    9010 FORMAT(/'     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F5.2 )
    FLUSH( stdout )

    ! Convergence threshold for iterative diagonalization is
    ! automatically updated during self consistency
    IF( iter > 1 ) THEN
      !
      IF ( iter == 2 ) ethr = 1.D-2
      ethr = MIN( ethr, 0.1D0*dr2 / MAX( 1.D0, nelec ) )
      ! do not allow convergence threshold to become too small:
      ! iterative diagonalization may become unstable
      ethr = MAX( ethr, 1.D-13 )
      !
    ENDIF

    first = ( iter == 1 )

    ! deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v> is calculated a
    ! first time here using the input density and potential ( to be
    ! used to calculate the Harris-Weinert-Foulkes energy )
    deband_hwf = my_delta_e()
    ! ffr: what are the inputs  for delta_e() ?

    ! save input density to rhoin
    CALL scf_type_COPY( rho, rhoin )

    scf_step: DO

      ! tr2_min is set to an estimate of the error on the energy
      ! due to diagonalization - used only for the first scf iteration
      tr2_min = 0.D0
      IF ( first ) tr2_min = ethr*MAX( 1.D0, nelec ) 

      write(*,*) 'first = ', first
      write(*,*) 'tr2_min = ', tr2_min
      write(*,*) 'ethr = ', ethr

      !
      ! diagonalization of the KS hamiltonian
      !
      IF ( lelfield ) THEN
        CALL c_bands_efield( iter )
      ELSE
        CALL my_c_bands( iter )
      ENDIF


      IF( stopped_by_user ) THEN
        conv_elec=.FALSE.
        CALL save_in_electrons( iter-1, dr2, ethr, et )
        GO TO 10
      ENDIF

      IF(one_atom_occupations) CALL new_evc()

      ! xk, wk, isk, et, wg are distributed across pools;
      ! the first node has a complete copy of xk, wk, isk,
      ! while eigenvalues et and weights wg must be
      ! explicitly collected to the first node
      ! this is done here for et, in sum_band for wg
      CALL poolrecover( et, nbnd, nkstot, nks )

      ! the new density is computed here. For PAW:
      ! sum_band computes new becsum (stored in uspp modules)
      ! and a subtly different copy in rho%bec (scf module)
      !CALL sum_band()
      call my_sum_band()

      ! the Harris-Weinert-Foulkes energy is computed here using only
      ! quantities obtained from the input density
      write(*,*)
      write(*,*) 'BEGIN: Calculating hwf_energy terms: '
      hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet
      write(*,*) 'eband = ', eband
      write(*,*) 'deband_hwf = ', deband_hwf
      write(*,*) 'etxc = ', etxc
      write(*,*) 'etxcc = ', etxcc ! ffr
      write(*,*) 'ewld = ', ewld
      write(*,*) 'ehart = ', ehart
      write(*,*) 'demet = ', demet
      write(*,*) 'END: Calculating hwf_energy terms: '


      IF( okpaw ) hwf_energy = hwf_energy + epaw


      ! calculate total and absolute magnetization
      IF ( lsda .OR. noncolin ) CALL my_compute_magnetization()

      ! eband  = \sum_v \epsilon_v    is calculated by sum_band
      ! deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v>
      ! eband + deband = \sum_v <\psi_v | T + Vion |\psi_v>
      deband = my_delta_e()

      ! mix_rho mixes several quantities: rho in g-space, tauk (for
      ! meta-gga), ns and ns_nc (for lda+u) and becsum (for paw)
      ! The mixing could in principle be done on pool 0 only, but
      ! mix_rho contains a call to rho_ddot that in the PAW case
      ! is parallelized on the entire image
      !
      CALL my_mix_rho( rho, rhoin, mixing_beta, dr2, tr2_min, iter, nmix, &
                       iunmix, conv_elec )

      ! Results are broadcast from pool 0 to others to prevent trouble
      ! on machines unable to yield the same results for the same 
      ! calculations on the same data, performed on different procs
      CALL bcast_scf_type( rhoin, root_pool, inter_pool_comm )
      CALL mp_bcast( dr2, root_pool, inter_pool_comm )
      CALL mp_bcast( conv_elec, root_pool, inter_pool_comm )

      IF (.NOT. scf_must_converge .AND. idum == niter) conv_elec = .TRUE.


      ! if convergence is achieved or if the self-consistency error
      ! (dr2) is smaller than the estimated error due to diagonalization
      ! (tr2_min), rhoin and rho are unchanged: rhoin contains the input
      ! density and rho contains the output density.
      ! In all other cases, rhoin contains the mixed charge density 
      ! (the new input density) while rho is unchanged
      !
      IF( first .and. nat > 0) THEN
        !
        ! first scf iteration: check if the threshold on diagonalization
        ! (ethr) was small enough wrt the error in self-consistency (dr2)
        ! if not, perform a new diagonalization with reduced threshold
        !
        first = .FALSE.
        !
        IF( dr2 < tr2_min ) THEN
          !
          WRITE( stdout, '(/,5X,"Threshold (ethr) on eigenvalues was ", &
                           &    "too large:",/,5X,                      &
                           & "Diagonalizing with lowered threshold",/)' )
          !
          ethr = 0.1D0*dr2 / MAX( 1.D0, nelec )
          !
          CYCLE scf_step
          !
        ENDIF
        !
      ENDIF
      !
      IF( .NOT. conv_elec ) THEN
        !
        ! no convergence yet: calculate new potential from mixed
        ! charge density (i.e. the new estimate)
        !
        CALL my_v_of_rho( rhoin, rho_core, rhog_core, &
                       ehart, etxc, vtxc, eth, etotefield, charge, v )
        !
        IF (okpaw) THEN
           CALL PAW_potential( rhoin%bec, ddd_paw, epaw,etot_cmp_paw )
           CALL PAW_symmetrize_ddd( ddd_paw )
        ENDIF
        !
        ! estimate correction needed to have variational energy:
        ! T + E_ion (eband + deband) are calculated in sum_band
        ! and delta_e using the output charge density rho;
        ! E_H (ehart) and E_xc (etxc) are calculated in v_of_rho
        ! above, using the mixed charge density rhoin%of_r.
        ! delta_escf corrects for this difference at first order
        !
        descf = my_delta_escf( rhoin, rho )
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
        write(*,*)
        write(*,*) '-------------------------------------------------------'
        write(*,*) 'VNEW: Convergence achived: calculating vnew to correct forces'

        vnew%of_r(:,:) = v%of_r(:,:)
        write(*,*) 'VNEW: Before my_v_of_rho: sum vnew (in Ha) = ', sum(vnew%of_r)*0.5d0
        !
        CALL my_v_of_rho( rho,rho_core,rhog_core, &
                       ehart, etxc, vtxc, eth, etotefield, charge, v )
        !
        write(*,*) 'VNEW: After my_v_of_rho: sum v (in Ha) = ', sum(v%of_r)*0.5d0
        !
        vnew%of_r(:,:) = v%of_r(:,:) - vnew%of_r(:,:)
        !
        write(*,*) 'VNEW: After subtraction: sum vnew%of_r (in Ha) = ', sum(vnew%of_r)*0.5d0
        write(*,*) '-------------------------------------------------------'
        write(*,*)
        !
        IF(okpaw) THEN
          CALL PAW_potential( rho%bec, ddd_paw, epaw, etot_cmp_paw )
          CALL PAW_symmetrize_ddd( ddd_paw )
        ENDIF
        !
        ! note that rho is here the output, not mixed, charge density
        ! so correction for variational energy is no longer needed
        !
        descf = 0._dp
        !
      ENDIF 
      !
      ! if we didn't cycle before we can exit the do-loop
      !
      EXIT scf_step
      !
    ENDDO scf_step
    !
    plugin_etot = 0.0_dp
    !
    CALL plugin_scf_energy(plugin_etot,rhoin)
    !
    CALL plugin_scf_potential(rhoin,conv_elec,dr2,vltot)
    !
    ! define the total local potential (external + scf)
    !
    CALL sum_vrs( dfftp%nnr, nspin, vltot, v%of_r, vrs )
    !
    ! interpolate the total local potential
    !
    CALL interpolate_vrs( dfftp%nnr, nspin, doublegrid, kedtau, v%kin_r, vrs )
    !
    ! in the US case we have to recompute the self-consistent
    ! term in the nonlocal potential
    ! PAW: newd contains PAW updates of NL coefficients
    !
    !call my_newd()
    call newd()
    !
    IF( lelfield ) en_el =  my_calc_pol( )
    !
    IF( report > 0 ) THEN
      IF( conv_elec .OR.  MOD(iter,report) == 0 ) CALL report_mag()
    ELSEIF ( report < 0 ) THEN
      IF( conv_elec ) CALL report_mag()
    ENDIF

    IF( conv_elec ) WRITE( stdout, 9101 )
 
    IF( conv_elec ) THEN 
      scf_error = dr2
      n_scf_steps = iter
    ENDIF  

    !
    IF( conv_elec .OR. MOD( iter, iprint ) == 0 ) THEN
      ! iverbosity == 0 for the PW code
      ! iverbosity >  2 for the HP code
      CALL print_ks_energies()
    ENDIF
    !
    IF( ABS( charge - nelec ) / charge > 1.D-7 ) THEN
      WRITE( stdout, 9050 ) charge, nelec
      IF( ABS( charge - nelec ) / charge > 1.D-3 ) THEN
        IF(.not.lgauss) THEN
          CALL errore( 'electrons', 'charge is wrong: smearing is needed', 1 )
        ELSE
          CALL errore( 'electrons', 'charge is wrong', 1 )
        ENDIF
      ENDIF
    ENDIF
    !
    etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet + descf
    ! for hybrid calculations, add the current estimate of exchange energy
    ! (it will subtracted later if exx_is_active to be replaced with a better estimate)
    etot = etot - exxen
    hwf_energy = hwf_energy - exxen ! [LP]
    !
    IF (okpaw) etot = etot + epaw

    !
    IF ( lelfield ) etot = etot + en_el
    ! not sure about the HWF functional in the above case
    IF( textfor ) THEN
      eext = alat*compute_eextfor()
      etot = etot + eext
      hwf_energy = hwf_energy + eext
    ENDIF
    
    IF( llondon ) THEN
        etot = etot + elondon
        hwf_energy = hwf_energy + elondon
    ENDIF
    !
    ! grimme-d3 dispersion energy
    IF(ldftd3) THEN
      etot = etot + edftd3
      hwf_energy = hwf_energy + edftd3
    ENDIF
     !
     ! calculate the xdm energy contribution with converged density
     IF (lxdm .and. conv_elec) THEN
        exdm = energy_xdm()  
        etot = etot + exdm
        hwf_energy = hwf_energy + exdm
     ENDIF
     IF (ts_vdw) THEN
        ! factor 2 converts from Ha to Ry units
        etot = etot + 2.0d0*EtsvdW
        hwf_energy = hwf_energy + 2.0d0*EtsvdW
     ENDIF
     !
     IF ( tefield ) THEN
        etot = etot + etotefield
        hwf_energy = hwf_energy + etotefield
     ENDIF
     ! TB gate energy
     IF ( gate) THEN
        etot = etot + etotgatefield
        hwf_energy = hwf_energy + etotgatefield
     ENDIF
     !
     IF ( lfcpopt .or. lfcpdyn ) THEN
        etot = etot + ef * tot_charge
        hwf_energy = hwf_energy + ef * tot_charge
     ENDIF
     !
     ! adds possible external contribution from plugins to the energy
     !
     etot = etot + plugin_etot 
     !
     CALL print_energies ( printout )
     !
     IF ( conv_elec ) THEN
        !
        ! if system is charged add a Makov-Payne correction to the energy
        ! (not in case of hybrid functionals: it is added at the end)
        !
        IF ( do_makov_payne .AND. printout/= 0 ) CALL makov_payne( etot )
        !
        ! print out ESM potentials if desired
        !
        IF ( do_comp_esm ) CALL esm_printpot( rho%of_g )
        !
        WRITE( stdout, 9110 ) iter
        !
        ! jump to the end
        !
        GO TO 10
        !
     ENDIF
     !
     ! uncomment the following line if you wish to monitor the evolution
     ! of the force calculation during self-consistency
     !
     !CALL forces()
     
     !
  ENDDO
  n_scf_steps = iter
  scf_error = dr2
  !
  WRITE( stdout, 9101 )
  WRITE( stdout, 9120 ) iter
  !
10  FLUSH( stdout )
  !
  ! exiting: write (unless disabled) the charge density to file
  ! (also write ldaU ns coefficients and PAW becsum)
  !
  IF ( io_level > -1 ) CALL write_scf( rho, nspin )
  !
  ! delete mixing info if converged, keep it if not
  !
  IF ( conv_elec ) THEN
     CALL close_mix_file( iunmix, 'delete' )
  ELSE
     CALL close_mix_file( iunmix, 'keep' )
  ENDIF
  !
  IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
  call destroy_scf_type ( rhoin )
  
  write(*,*)
  write(*,*) '</div> EXIT my_electrons_scf'
  write(*,*)

  !
  RETURN
  !
  ! formats
  !

9001 FORMAT(/'     per-process dynamical memory: ',f7.1,' Mb' )
9050 FORMAT(/'     WARNING: integrated charge=',F15.8,', expected=',F15.8 )
9101 FORMAT(/'     End of self-consistent calculation' )
9110 FORMAT(/'     convergence has been achieved in ',i3,' iterations' )
9120 FORMAT(/'     convergence NOT achieved after ',i3,' iterations: stopping' )
  !
  CONTAINS





! Inner subroutine
 
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
     WRITE( stdout, 9081 ) etot
     IF ( only_paw ) WRITE( stdout, 9085 ) etot+total_core_energy
     IF ( iverbosity > 1 ) WRITE( stdout, 9082 ) hwf_energy
     IF ( dr2 > eps8 ) THEN
        WRITE( stdout, 9083 ) dr2
     ELSE
        WRITE( stdout, 9084 ) dr2
     END IF
     IF ( lgauss ) then
        WRITE( stdout, 9070 ) demet
        WRITE( stdout, 9170 ) etot-demet
        WRITE( stdout, 9061 )
     ELSE
        WRITE( stdout, 9060 )
     END IF
     WRITE( stdout, 9062 ) (eband + deband), ehart, ( etxc - etxcc ), ewld
     !
     IF ( llondon ) WRITE ( stdout , 9074 ) elondon
     IF ( ldftd3 )  WRITE ( stdout , 9078 ) edftd3
     IF ( lxdm )    WRITE ( stdout , 9075 ) exdm
     IF ( ts_vdw )  WRITE ( stdout , 9076 ) 2.0d0*EtsvdW
     IF ( textfor)  WRITE ( stdout , 9077 ) eext
     IF ( tefield )            WRITE( stdout, 9064 ) etotefield
     IF ( gate )               WRITE( stdout, 9065 ) etotgatefield
     IF ( ABS (descf) > eps8 ) WRITE( stdout, 9069 ) descf
     IF ( okpaw ) THEN
       WRITE( stdout, 9067 ) epaw
       ! Detailed printout of PAW energy components, if verbosity is high
       IF(iverbosity>0)THEN
       WRITE( stdout, 9068) SUM(etot_cmp_paw(:,1,1)), &
                            SUM(etot_cmp_paw(:,1,2)), &
                            SUM(etot_cmp_paw(:,2,1)), &
                            SUM(etot_cmp_paw(:,2,2)), &
       SUM(etot_cmp_paw(:,1,1))+SUM(etot_cmp_paw(:,1,2))+ehart, &
       SUM(etot_cmp_paw(:,2,1))+SUM(etot_cmp_paw(:,2,2))+etxc-etxcc
       ENDIF
     ENDIF
     !
     ! With Fermi-Dirac population factor, etot is the electronic
     ! free energy F = E - TS , demet is the -TS contribution
     !
     !
     ! With Fictitious charge particle (FCP), etot is the grand
     ! potential energy Omega = E - muN, -muN is the potentiostat
     ! contribution.
     !
     IF ( lfcpopt .OR. lfcpdyn ) WRITE( stdout, 9072 ) ef*tot_charge
     !
  ELSE IF ( conv_elec ) THEN
     !
     WRITE( stdout, 9081 ) etot
     IF ( iverbosity > 1 ) WRITE( stdout, 9082 ) hwf_energy
     IF ( dr2 > eps8 ) THEN
        WRITE( stdout, 9083 ) dr2
     ELSE
        WRITE( stdout, 9084 ) dr2
     END IF
     IF ( lgauss ) then
        WRITE( stdout, 9070 ) demet
        WRITE( stdout, 9170 ) etot-demet
     ENDIF
     !
  ELSE
     !
     WRITE( stdout, 9080 ) etot
     IF ( iverbosity > 1 ) WRITE( stdout, 9082 ) hwf_energy
     IF ( dr2 > eps8 ) THEN
        WRITE( stdout, 9083 ) dr2
     ELSE
        WRITE( stdout, 9084 ) dr2
     END IF
  ENDIF
  !
  CALL plugin_print_energies()
  !
  IF ( lsda ) WRITE( stdout, 9017 ) magtot, absmag
  !
  IF ( noncolin .AND. domag ) &
       WRITE( stdout, 9018 ) magtot_nc(1:3), absmag
  !
  IF ( i_cons == 3 .OR. i_cons == 4 )  &
       WRITE( stdout, 9071 ) bfield(1), bfield(2), bfield(3)
  IF ( i_cons /= 0 .AND. i_cons < 4 ) &
       WRITE( stdout, 9073 ) lambda
  !
  FLUSH( stdout )
  !
  RETURN
       !
9017 FORMAT(/'     total magnetization       =', F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9018 FORMAT(/'     total magnetization       =',3F9.2,' Bohr mag/cell' &
       &   ,/'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9060 FORMAT(/'     The total energy is the sum of the following terms:' )
9061 FORMAT(/'     The total energy is F=E-TS. E is the sum of the following terms:' )
9062 FORMAT( '     one-electron contribution =',F17.8,' Ry' &
            /'     hartree contribution      =',F17.8,' Ry' &
            /'     xc contribution           =',F17.8,' Ry' &
            /'     ewald contribution        =',F17.8,' Ry' )
9064 FORMAT( '     electric field correction =',F17.8,' Ry' )
9065 FORMAT( '     gate field correction     =',F17.8,' Ry' ) ! TB

9067 FORMAT( '     one-center paw contrib.   =',F17.8,' Ry' )
9068 FORMAT( '      -> PAW hartree energy AE =',F17.8,' Ry' &
            /'      -> PAW hartree energy PS =',F17.8,' Ry' &
            /'      -> PAW xc energy AE      =',F17.8,' Ry' &
            /'      -> PAW xc energy PS      =',F17.8,' Ry' &
            /'      -> total E_H with PAW    =',F17.8,' Ry' &
            /'      -> total E_XC with PAW   =',F17.8,' Ry' )
9069 FORMAT( '     scf correction            =',F17.8,' Ry' )
9070 FORMAT( '     smearing contrib. (-TS)   =',F17.8,' Ry' )
9071 FORMAT( '     Magnetic field            =',3F12.7,' Ry' )
9072 FORMAT( '     pot.stat. contrib. (-muN) =',F17.8,' Ry' )
9073 FORMAT( '     lambda                    =',F11.2,' Ry' )
9074 FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
9075 FORMAT( '     Dispersion XDM Correction =',F17.8,' Ry' )
9076 FORMAT( '     Dispersion T-S Correction =',F17.8,' Ry' )
9077 FORMAT( '     External forces energy    =',F17.8,' Ry' )
9078 FORMAT( '     DFT-D3 Dispersion         =',F17.8,' Ry' )
9080 FORMAT(/'     total energy              =',0PF17.8,' Ry' )
9081 FORMAT(/'!    total energy              =',0PF17.8,' Ry' )
9082 FORMAT( '     Harris-Foulkes estimate   =',0PF17.8,' Ry' )
9083 FORMAT( '     estimated scf accuracy    <',0PF17.8,' Ry' )
9084 FORMAT( '     estimated scf accuracy    <',1PE17.1,' Ry' )
9085 FORMAT(/'     total all-electron energy =',0PF17.6,' Ry' )
9170 FORMAT( '     internal energy E=F+TS    =',0PF17.8,' Ry' )
  END SUBROUTINE print_energies
  !
END SUBROUTINE