!----------------------------------------------------------------------------
SUBROUTINE my_v_of_rho( rho, rho_core, rhog_core, &
                        ehart, etxc, vtxc, eth, etotefield, charge, v )
!----------------------------------------------------------------------------
  !! This routine computes the Hartree and Exchange and Correlation
  !! potential and energies which corresponds to a given charge density
  !! The XC potential is computed in real space, while the
  !! Hartree potential is computed in reciprocal space.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  USE ions_base,        ONLY : nat, tau
  USE ldaU,             ONLY : lda_plus_u, lda_plus_u_kind, ldmx_b, &
                               nsg, v_nsg 
  USE funct,            ONLY : dft_is_meta, get_meta
  USE scf,              ONLY : scf_type
  USE cell_base,        ONLY : alat
  USE control_flags,    ONLY : ts_vdw
  USE tsvdw_module,     ONLY : tsvdw_calculate, UtsvdW
  USE cell_base,        ONLY : omega   ! ffr
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  TYPE(scf_type), INTENT(INOUT) :: v
  !! the scf (Hxc) potential 
  !=================> NB: NOTE that in F90 derived data type must be INOUT and 
  !=================> not just OUT because otherwise their allocatable or pointer
  !=================> components are NOT defined 
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
  !! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
  !! the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: vtxc
  !! the integral V_xc * rho
  REAL(DP), INTENT(OUT) :: etxc
  !! the E_xc energy
  REAL(DP), INTENT(OUT) :: ehart
  !! the hartree energy
  REAL(DP), INTENT(OUT) :: eth
  !! the hubbard energy
  REAL(DP), INTENT(OUT) :: charge
  !! the integral of the charge
  REAL(DP) :: eth1
  !! the hubbard energy coming from the background states
  REAL(DP), INTENT(INOUT) :: etotefield
  !! electric field energy - inout due to the screwed logic of add_efield
  !
  INTEGER :: is, ir

  write(*,*)
  write(*,*) '<div> ENTER my_v_of_rho'
  write(*,*)

  write(*,*) 'sum v_of_r before my_v_xc, should be zero (in Ha): ', sum(v%of_r)*0.5d0

  ! calculate exchange-correlation potential
  IF (dft_is_meta()) then
    CALL my_v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v%of_r, v%kin_r )   
    !v%kin_r(:,:) = 0.d0 ! ffr
  ELSE
    CALL my_v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
  ENDIF

  write(*,*) 'sum v_of_r after my_v_xc (in Ha):', sum(v%of_r)*0.5d0

  ! add a magnetic field  (if any)
  CALL add_bfield( v%of_r, rho%of_r )

  ! calculate hartree potential
  CALL my_v_h( rho%of_g(:,1), ehart, charge, v%of_r )

  write(*,*) 'sum v_of_r after my_v_h (in Ha):', sum(v%of_r)*0.5d0

  ! DFT+U(+V): build up (extended) Hubbard potential 
  IF (lda_plus_u) THEN
    stop 'lda_plus_u is disable in my_v_of_rho'
  ENDIF

  write(*,*)
  write(*,*) 'my_v_of_rho: Ehartree (in Ha) = ', ehart*0.5d0
  write(*,*) 'my_v_of_rho: Exc      (in Ha) = ', etxc*0.5d0
  write(*,*) 'my_v_of_rho: Evtxc    (in Ha) = ', vtxc*0.5d0
  write(*,*)

  ! add an electric field
  ! 
  DO is = 1, nspin_lsda
    CALL add_efield(v%of_r(1,is), etotefield, rho%of_r(:,1), .false. )
  ENDDO

  ! add Tkatchenko-Scheffler potential (factor 2: Ha -> Ry)
  IF (ts_vdw) THEN
     CALL tsvdw_calculate(tau*alat,rho%of_r(:,1))
     DO is = 1, nspin_lsda
        DO ir=1,dfftp%nnr
           v%of_r(ir,is) = v%of_r(ir,is) + 2.0d0*UtsvdW(ir)
        END DO
     END DO
  END IF

  write(*,*)
  write(*,*) '</div> EXIT my_v_of_rho'
  write(*,*)

  RETURN
  !
END SUBROUTINE my_v_of_rho
!

