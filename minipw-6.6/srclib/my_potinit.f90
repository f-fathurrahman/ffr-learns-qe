!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE my_potinit()
  !----------------------------------------------------------------------------
  !
  ! ... This routine initializes the self consistent potential in the array
  ! ... vr. There are three possible cases:
  !
  ! ... a) the code is restarting from a broken run:
  ! ...    read rho from data stored during the previous run
  ! ... b) the code is performing a non-scf calculation following a scf one:
  ! ...    read rho from the file produced by the scf calculation
  ! ... c) the code starts a new calculation:
  ! ...    calculate rho as a sum of atomic charges
  ! 
  ! ... In all cases the scf potential is recalculated and saved in vr
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : pi
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : alat, omega
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE basis,                ONLY : starting_pot
  USE klist,                ONLY : nelec
  USE lsda_mod,             ONLY : lsda, nspin
  USE lsda_mod,             ONLY : starting_magnetization
  USE fft_base,             ONLY : dfftp, dffts
  USE gvect,                ONLY : ngm, gstart, g, gg, ig_l2g
  USE gvecs,                ONLY : doublegrid
  USE control_flags,        ONLY : lscf, gamma_only
  USE scf,                  ONLY : rho, rho_core, rhog_core, &
                                   vltot, v, vrs, kedtau
  USE funct,                ONLY : dft_is_meta
  USE ener,                 ONLY : ehart, etxc, vtxc, epaw
  USE ldaU,                 ONLY : lda_plus_u, Hubbard_lmax, eth, &
                                   niter_with_fixed_ns, lda_plus_u_kind, &
                                   nsg, nsgnew
  USE noncollin_module,     ONLY : noncolin, report
  USE io_files,             ONLY : restart_dir, input_drho, check_file_exist
  USE spin_orb,             ONLY : domag, lforcet
  USE mp,                   ONLY : mp_sum
  USE mp_bands ,            ONLY : intra_bgrp_comm, root_bgrp
  USE io_global,            ONLY : ionode, ionode_id
  USE io_rho_xml,           ONLY : read_scf
  USE io_base,              ONLY : read_rhog
  USE fft_rho,              ONLY : rho_g2r, rho_r2g
  !
  USE uspp,                 ONLY : becsum
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_init,             ONLY : PAW_atomic_becsum
  USE paw_onecenter,        ONLY : PAW_potential
  !
  IMPLICIT NONE
  !
  REAL(DP)              :: charge           ! the starting charge
  REAL(DP)              :: etotefield       !
  REAL(DP)              :: fact
  INTEGER               :: is
  LOGICAL               :: exst 
  CHARACTER(LEN=320)    :: filename

  write(*,*)
  write(*,*) '<div> ENTER my_potinit'
  write(*,*)


  !
  filename = TRIM (restart_dir( )) // 'charge-density'
  exst     = check_file_exist( TRIM(filename) // '.dat' )

  !
  IF( starting_pot == 'file' .AND. exst ) THEN
    !
    stop 'Disabled in my_potinit'
    !
  ELSE
    !
    ! ... Case c): the potential is built from a superposition 
    ! ... of atomic charges contained in the array rho_at
    !
    WRITE( UNIT = stdout, &
           FMT = '(/5X,"Initial potential from superposition of free atoms")' )
    !
    CALL my_atomic_rho_g( rho%of_g, nspin )

    ! in the DFT+U(+V) case set the initial value of ns (or nsg)
    !
    IF( lda_plus_u ) THEN
      stop 'lda_plus_u is disabled in my_potinit'
    ENDIF

    ! ... in the paw case uses atomic becsum
    IF( okpaw )      CALL PAW_atomic_becsum()
  
    IF( input_drho /= ' ' ) THEN
      write(*,*) 'input_drho = ', input_drho
      stop 'input_drho is disabled in my_potinit'
    ENDIF
    !
  ENDIF
  !
  ! ... check the integral of the starting charge, renormalize if needed
  !
  charge = 0.D0
  IF( gstart == 2 ) THEN
    charge = omega*REAL( rho%of_g(1,1) )
  ENDIF
  CALL mp_sum( charge , intra_bgrp_comm )

  write(*,*) 'my_potinit: charge = ', charge

  !
  IF ( lscf .AND. ABS( charge - nelec ) > ( 1.D-7 * charge ) ) THEN
    !
    IF( charge > 1.D-8 .AND. nat > 0 ) THEN
      WRITE( stdout, '(/,5X,"starting charge ",F10.5, &
                       & ", renormalised to ",F10.5)') charge, nelec
      rho%of_g = rho%of_g / charge * nelec
    ELSE 
      WRITE( stdout, '(/,5X,"Starting from uniform charge")')
      rho%of_g(:,1:nspin) = (0.0_dp,0.0_dp)
      IF( gstart == 2 ) rho%of_g(1,1) = nelec / omega
      if(nspin == 2) then
        rho%of_g(1,2) = sum(starting_magnetization(:))
      endif
    ENDIF
    !
  ELSEIF( .NOT. lscf .AND. ABS( charge - nelec ) > (1.D-3 * charge ) ) THEN
    !
    CALL errore( 'potinit', 'starting and expected charges differ', 1 )
    !
  ENDIF
  !
  ! ... bring starting rho from G- to R-space
  !
  CALL rho_g2r(dfftp, rho%of_g, rho%of_r)

  write(*,*) 'my_potinit: integ rho%of_r(:,1): ', &
    sum(rho%of_r(:,1))*omega/(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)

  !
  IF( dft_is_meta() ) THEN
     IF (starting_pot /= 'file') THEN
        ! ... define a starting (TF) guess for rho%kin_r from rho%of_r
        ! ... to be verified for LSDA: rho is (tot,magn), rho_kin is (up,down)
        fact = (3.d0/5.d0)*(3.d0*pi*pi)**(2.0/3.0)
        DO is = 1, nspin
           rho%kin_r(:,is) = fact * abs(rho%of_r(:,is)*nspin)**(5.0/3.0)/nspin
        END DO
        !if (nspin==2) then
        !     rho%kin_r(:,1) = fact * abs(rho%of_r(:,1)+rho%of_r(:,2))**(5.0/3.0)/2.0
        !     rho%kin_r(:,2) = fact * abs(rho%of_r(:,1)-rho%of_r(:,2))**(5.0/3.0)/2.0
        !endif
        ! ... bring it to g-space
        CALL rho_r2g(dfftp, rho%kin_r, rho%kin_g)
     ELSE
        ! ... rho%kin was read from file in G-space, bring it to R-space
        CALL rho_g2r(dfftp, rho%kin_g, rho%kin_r)
     ENDIF
     !
  ENDIF

  ! plugin contribution to local potential
  CALL plugin_scf_potential(rho, .FALSE., -1.d0, vltot)
  
  ! compute the potential and store it in v
  CALL my_v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )

  IF( okpaw ) then
    write(*,*) 'my_potinit: sum ddd_PAW before PAW_potential = ', sum(ddd_PAW)
    CALL PAW_potential(rho%bec, ddd_PAW, epaw)
    write(*,*) 'my_potinit: sum rho%bec = ', sum(rho%bec)
    write(*,*) 'my_potinit: EHxc (in Ha) = ', 0.5d0*epaw
    write(*,*) 'my_potinit: sum ddd_PAW after PAW_potential (in Ha) = ', 0.5d0*sum(ddd_PAW)
    write(*,*) 'Wrote ddd_PAW to fort.999, shape ddd_PAW: ', shape(ddd_PAW)
  endif

  ! define the total local potential (external+scf)
  CALL my_set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )

  ! vrs is a global variable  

  !
  ! ... write on output the parameters used in the DFT+U(+V) calculation
  !
  IF( lda_plus_u ) THEN
    stop 'lda_plus_u is disabled in my_potinit'
  ENDIF
  !
  IF( report /= 0 .AND. &
      noncolin .AND. domag .AND. lscf ) CALL report_mag()

  write(*,*)
  write(*,*) '</div> EXIT my_potinit'
  write(*,*)

  RETURN
  !
END SUBROUTINE my_potinit

