!----------------------------------------------------------------------------
SUBROUTINE my_init_run()
!----------------------------------------------------------------------------
  !
  USE klist,              ONLY : nkstot
  USE symme,              ONLY : sym_rho_init
  USE wvfct,              ONLY : nbnd, et, wg, btype
  USE control_flags,      ONLY : lmd, gamma_only, smallmem, ts_vdw
  USE gvect,              ONLY : g, gg, mill, gcutm, ig_l2g, ngm, ngm_g, &
                                 gshells, gstart ! to be comunicated to the Solvers if gamma_only
  USE gvecs,              ONLY : gcutms, ngms
  USE cell_base,          ONLY : at, bg, set_h_ainv
  USE cellmd,             ONLY : lmovecell
  USE dynamics_module,    ONLY : allocate_dyn_vars
  USE paw_variables,      ONLY : okpaw
  USE paw_init,           ONLY : paw_init_onecenter, allocate_paw_internals
  USE bp,                 ONLY : allocate_bp_efield, bp_global_map
  USE fft_base,           ONLY : dfftp, dffts
  USE funct,              ONLY : dft_is_hybrid
  USE recvec_subs,        ONLY : ggen, ggens
  USE my_dfunct,             ONLY : my_newd
  USE esm,                ONLY : do_comp_esm, esm_init
  USE tsvdw_module,       ONLY : tsvdw_initialize
  USE Coul_cut_2D,        ONLY : do_cutoff_2D, cutoff_fact 
  !
  IMPLICIT NONE

  write(*,*)
  write(*,*) '******************************************************************'
  write(*,*) '                          ENTER init_run                          '
  write(*,*) '******************************************************************'


  !
  ! ... calculate limits of some indices, used in subsequent allocations
  !
  CALL my_pre_init()
  !
  ! ... determine the data structure for fft arrays
  !
  CALL data_structure( gamma_only )
  !
  ! ... print a summary and a memory estimate before starting allocating
  !
  CALL summary()
  CALL memory_report()
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL allocate_fft()
  !
  ! ... generate reciprocal-lattice vectors and fft indices
  !
  IF( smallmem ) THEN
    CALL ggen( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
         g, gg, mill, ig_l2g, gstart, no_global_sort = .TRUE. )
  ELSE
    CALL ggen( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
      g, gg, mill, ig_l2g, gstart )
  ENDIF
  !
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms )
  if (gamma_only) THEN
     ! ... Solvers need to know gstart
     call export_gstart_2_solvers(gstart)
  END IF
  !
  IF(do_comp_esm) CALL esm_init()
  !
  ! ... setup the 2D cutoff factor
  !
  IF(do_cutoff_2D) CALL cutoff_fact()
  !
  CALL gshells( lmovecell )
  !
  ! ... variable initialization for parallel symmetrization
  !
  CALL sym_rho_init( gamma_only )
  !
  ! ... allocate memory for all other arrays (potentials, wavefunctions etc)
  !
  CALL allocate_nlpot()
  IF( okpaw ) THEN
    CALL allocate_paw_internals() ! ddd_paw
    CALL paw_init_onecenter()
  ENDIF
  CALL allocate_locpot()
  CALL allocate_bp_efield()
  CALL bp_global_map()
  !
  call plugin_initbase()
  !
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ), btype( nbnd, nkstot ) )
  !
  et(:,:) = 0.D0
  wg(:,:) = 0.D0
  !
  btype(:,:) = 1
  !
  IF (ts_vdw) THEN
     CALL tsvdw_initialize()
     CALL set_h_ainv()
  END IF
  !
  CALL allocate_wfc_k()
  CALL openfil()
  !
  CALL my_hinit0()
  !
  CALL my_potinit()
  !
  CALL my_newd()
  !
  CALL my_wfcinit()
  !
  IF ( lmd ) CALL allocate_dyn_vars()

  write(*,*)
  write(*,*) '******************************************************************'
  write(*,*) '                          EXIT init_run                          '
  write(*,*) '******************************************************************'

  !
  RETURN
  !
END SUBROUTINE my_init_run


!----------------------------------------------------------------------------
SUBROUTINE my_pre_init()
!----------------------------------------------------------------------------
  !
  USE ions_base,        ONLY : nat, nsp, ityp
  USE uspp_param,       ONLY : upf, lmaxkb, nh, nhm, nbetam
  USE uspp,             ONLY : nkb, nkbus
  IMPLICIT NONE
  INTEGER :: na, nt, nb
  !
  !     calculate the number of beta functions for each atomic type
  !
  lmaxkb = -1
  DO nt = 1, nsp
     !
     nh (nt) = 0
     !
     ! do not add any beta projector if pseudo in 1/r fmt (AF)
     IF ( upf(nt)%tcoulombp ) CYCLE 
     !
     DO nb = 1, upf(nt)%nbeta
        nh (nt) = nh (nt) + 2 * upf(nt)%lll(nb) + 1
        lmaxkb = MAX (lmaxkb, upf(nt)%lll(nb) )
     ENDDO
     !
  ENDDO
  !
  ! calculate the maximum number of beta functions
  !
  nhm = MAXVAL (nh (1:nsp))
  nbetam = MAXVAL (upf(:)%nbeta)
  !
  ! calculate the number of beta functions of the solid
  !
  nkb = 0
  nkbus = 0
  do na = 1, nat
     nt = ityp(na)
     nkb = nkb + nh (nt)
     if (upf(nt)%tvanp) nkbus = nkbus + nh (nt)
  enddo


END SUBROUTINE my_pre_init
