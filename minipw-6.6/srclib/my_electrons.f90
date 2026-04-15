SUBROUTINE my_electrons()
!----------------------------------------------------------------------------
!! General self-consistency loop, also for hybrid functionals.  
!! For non-hybrid functionals it just calls \(\texttt{electron_scf}\).

  USE kinds,                ONLY : DP
  USE check_stop,           ONLY : check_stop_now, stopped_by_user
  USE io_global,            ONLY : stdout, ionode
  USE fft_base,             ONLY : dfftp
  USE gvecs,                ONLY : doublegrid
  USE lsda_mod,             ONLY : nspin
  USE ener,                 ONLY : etot, hwf_energy, eband, ehart, &
                                   vtxc, etxc, etxcc, ewld, epaw, &
                                   elondon, edftd3, ef_up, ef_dw
  USE scf,                  ONLY : rho, rho_core, rhog_core, v, vltot, vrs, &
                                   kedtau, vnew
  USE control_flags,        ONLY : tr2, niter, conv_elec, restart, do_makov_payne
  USE io_files,             ONLY : iunres, seqopn
  USE ldaU,                 ONLY : eth
  USE extfield,             ONLY : etotefield
  USE wvfct,                ONLY : nbnd, wg, et
  USE klist,                ONLY : nks
  USE my_exx,               ONLY : aceinit,exxinit, exxenergy2, exxbuff, &
                                   fock0, fock1, fock2, fock3, dexx, use_ace, local_thr 
  USE funct,                ONLY : dft_is_hybrid, exx_is_active
  USE control_flags,        ONLY : adapt_thr, tr2_init, tr2_multi, gamma_only
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE ions_base,            ONLY : nat
  USE loc_scdm_k,           ONLY : localize_orbitals_k
  !
  IMPLICIT NONE
  REAL(DP) :: charge
  !! the total charge
  REAL(DP) :: exxen
  !! used to compute exchange energy
  REAL(DP), EXTERNAL :: exxenergyace
  INTEGER :: idum
  !! dummy counter on iterations
  INTEGER :: iter
  !! counter on iterations
  INTEGER :: printout, ik
  !
  REAL(DP) :: tr2_final
  !! final threshold for exx minimization 
  !! when using adaptive thresholds.
  LOGICAL :: first, exst
  REAL(DP) :: etot_cmp_paw(nat,2,2)
  LOGICAL :: DoLoc
  !
  !
  DoLoc = local_thr > 0.0d0
  exxen = 0.0d0
  iter = 0
  first = .true.
  tr2_final = tr2

  IF( restart ) THEN
   ! unsupported
   stop 'restart is not supported'
  ENDIF

  DO idum = 1,niter
    !
    iter = iter + 1
    !
    ! Self-consistency loop. For hybrid functionals the exchange potential
    ! is calculated with the orbitals at previous step (none at first step)
    !
    CALL my_electrons_scf (printout, exxen)
    ! Early return
    IF( .NOT. dft_is_hybrid() ) then
      RETURN
    endif
    !
    ! From now on: hybrid DFT only
    !
    IF( stopped_by_user .OR. .NOT. conv_elec ) THEN
      conv_elec = .FALSE.
      IF ( .NOT. first) THEN
        WRITE(stdout,'(5x,"Calculation (EXX) stopped during iteration #", &
                     & i6)') iter
        CALL seqopn (iunres, 'restart_e', 'formatted', exst)
        WRITE (iunres, *) iter-1, tr2, dexx
        WRITE (iunres, *) exxen, fock0, fock1, fock2
        WRITE (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
        WRITE (iunres, *) (et(1:nbnd,ik),ik=1,nks)
        CLOSE (unit=iunres, status='keep')
        CALL seqopn(iunres, 'restart_exx', 'unformatted', exst)
        WRITE (iunres) exxbuff
        CLOSE (unit=iunres, status='keep')
      ENDIF
      RETURN
    ENDIF
    !
    first =  first .AND. .NOT. exx_is_active ( )
    !
    ! "first" is true if the scf step was performed without exact exchange
    !
    IF ( first ) THEN
      !
      first = .false.
      !
      ! Activate exact exchange, set orbitals used in its calculation,
      ! then calculate exchange energy (will be useful at next step)
      !
      CALL exxinit(DoLoc)
      IF( DoLoc .and. gamma_only) THEN
        CALL localize_orbitals( )
      ELSEIF( DoLoc ) THEN
        CALL localize_orbitals_k( )
      ENDIF 
      IF( use_ace )THEN
        CALL aceinit ( DoLoc ) 
        fock2 = exxenergyace()
      ELSE
        fock2 = exxenergy2()
      ENDIF
      exxen = 0.50d0*fock2 
      etot = etot - etxc 
      !
      ! Recalculate potential because XC functional has changed,
      ! start self-consistency loop on exchange
      CALL v_of_rho( rho, rho_core, rhog_core, ehart, etxc, vtxc, eth, etotefield, charge, v)
      etot = etot + etxc + exxen
      !
      IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, epaw,etot_cmp_paw)
      CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
      !
    ELSE
      !
      ! fock1 is the exchange energy calculated for orbitals at step n,
      !       using orbitals at step n-1 in the expression of exchange
      !
      IF (use_ace) THEN
        fock1 = exxenergyace()
      ELSE
        fock1 = exxenergy2()
      ENDIF
      !
      ! Set new orbitals for the calculation of the exchange term
      !
      CALL exxinit( DoLoc )
      IF( DoLoc .and. gamma_only) THEN
        CALL localize_orbitals( )
      ELSE IF (DoLoc) THEN
        CALL localize_orbitals_k( )
      ENDIF 
      IF( use_ace ) then
        CALL aceinit( DoLoc, fock3 )
      ENDIF
      !
      ! fock2 is the exchange energy calculated for orbitals at step n,
      !       using orbitals at step n in the expression of exchange 
      ! fock0 is fock2 at previous step
      !
      fock0 = fock2
      IF( use_ace ) THEN
        fock2 = exxenergyace()
      ELSE
        fock2 = exxenergy2()
      ENDIF
      !
      ! check for convergence. dexx is positive definite: if it isn't,
      ! there is some numerical problem. One such cause could be that
      ! the treatment of the divergence in exact exchange has failed. 
      ! FIXME: to be properly implemented for all cases
      !
      IF ( DoLoc ) THEN
        dexx =  0.5D0 *( (fock1-fock0) + (fock3-fock2) ) 
      ELSE
        dexx = fock1 - 0.5D0*(fock0+fock2)
      ENDIF 
      !
      IF ( dexx < 0.0_dp ) THEN
        IF( Doloc ) THEN
          WRITE(stdout,'(5x,a,1e12.3)') "BEWARE: negative dexx:", dexx
          dexx = ABS ( dexx )
        ELSE
          CALL errore( 'electrons', 'dexx is negative! &
                  & Check that exxdiv_treatment is appropriate for the system,&
                  & or ecutfock may be too low', 1 )
        ENDIF
      ENDIF
      !
      ! remove the estimate exchange energy exxen used in the inner SCF
      !
      etot = etot + exxen + 0.5D0*fock2 - fock1
      hwf_energy = hwf_energy + exxen + 0.5D0*fock2 - fock1 ! [LP]
      exxen = 0.5D0*fock2 
      !
      IF ( dexx < tr2_final ) THEN
        WRITE( stdout, 9066 ) '!!', etot, hwf_energy
      ELSE
        WRITE( stdout, 9066 ) '  ', etot, hwf_energy
      ENDIF
      IF( dexx > 1.d-8 ) THEN
        WRITE( stdout, 9067 ) dexx
      ELSE
        WRITE( stdout, 9068 ) dexx
      ENDIF
      !
      WRITE( stdout, 9062 ) - fock1
      IF (use_ace) THEN
        WRITE( stdout, 9063 ) 0.5D0*fock2
      ELSE
        WRITE( stdout, 9064 ) 0.5D0*fock2
      ENDIF
      !
      IF ( dexx < tr2_final ) THEN
        IF ( do_makov_payne ) CALL makov_payne( etot )
        WRITE( stdout, 9101 )
        RETURN
      ENDIF
      !
      IF( adapt_thr ) THEN
        tr2 = MAX(tr2_multi * dexx, tr2_final)
        WRITE( stdout, 9121 ) tr2
      ENDIF
    ENDIF
    !
    WRITE( stdout,'(/5x,"EXX: now go back to refine exchange calculation")')
    !
    IF ( check_stop_now() ) THEN
      WRITE(stdout,'(5x,"Calculation (EXX) stopped after iteration #", i6)') iter
      conv_elec=.FALSE.
      CALL seqopn (iunres, 'restart_e', 'formatted', exst)
      WRITE (iunres, *) iter, tr2, dexx
      WRITE (iunres, *) exxen, fock0, fock1, fock2
      ! FIXME: et and wg are written to xml file
      WRITE (iunres, *) (wg(1:nbnd,ik),ik=1,nks)
      WRITE (iunres, *) (et(1:nbnd,ik),ik=1,nks)
      CLOSE (unit=iunres, status='keep')
      RETURN
    ENDIF
    !
  ENDDO
  !
  WRITE( stdout, 9120 ) iter
  FLUSH( stdout )
  !
  RETURN
  !
  ! formats
  !
9062 FORMAT( '     - averaged Fock potential =',0PF17.8,' Ry' )
9063 FORMAT( '     + Fock energy (ACE)       =',0PF17.8,' Ry' )
9064 FORMAT( '     + Fock energy (full)      =',0PF17.8,' Ry' )
9066 FORMAT(/,A2,'   total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' )
9067 FORMAT('     est. exchange err (dexx)  =',0PF17.8,' Ry' )
9068 FORMAT('     est. exchange err (dexx)  =',1PE17.1,' Ry' )
9101 FORMAT(/'     EXX self-consistency reached' )
9120 FORMAT(/'     EXX convergence NOT achieved after ',i3,' iterations: stopping' )
9121 FORMAT(/'     scf convergence threshold =',1PE17.1,' Ry' )

  RETURN
END SUBROUTINE









