INCLUDE 'prepare_all.f90'


PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL debug_PAW_xc_potential()
END PROGRAM



!----------------------------------
SUBROUTINE debug_PAW_xc_potential()
!----------------------------------
  use kinds, only : DP
  use PAW_variables, only : okpaw
  use scf, only : rho
  use paw_init, only : PAW_atomic_becsum
  USE atom, ONLY : g => rgrid
  USE ions_base, ONLY : nat, ityp
  USE lsda_mod, ONLY : nspin
  USE uspp_param, ONLY : upf
  USE PAW_variables, only : paw_info
  use paw_onecenter, only : paw_rho_lm
  !
  implicit none
  !
  REAL(DP) :: energy !! if present compute E[rho]
  REAL(DP) :: e_cmp(nat,2,2) !! components of the energy
  !
  ! ... local variables
  !
  INTEGER, PARAMETER :: AE = 1, PS = 2,&   ! All-Electron and Pseudo
                        H  = 1, XC = 2     ! Hartree and XC
  REAL(DP), POINTER :: rho_core(:)         ! pointer to AE/PS core charge density 
  TYPE(paw_info) :: i              ! minimal info on atoms
  INTEGER :: i_what                ! counter on AE and PS
  INTEGER :: ia          ! atoms counters and indexes
  INTEGER :: l2, kkbeta, imesh
  !
  REAL(DP), ALLOCATABLE :: v_lm(:,:,:)        ! workspace: potential
  REAL(DP), ALLOCATABLE :: rho_lm(:,:,:)      ! density expanded on Y_lm
  REAL(DP), ALLOCATABLE :: savedv_lm(:,:,:)   ! workspace: potential
  ! fake cross band occupations to select only one pfunc at a time:
  REAL(DP) :: energy_tot
  REAL(DP) :: sgn          ! +1 for AE -1 for PS

  energy_tot = 0.d0

  write(*,*)
  write(*,*) '----------------------------'
  write(*,*) 'ENTER debug_PAW_xc_potential'
  write(*,*) '----------------------------'

  if( .not. okpaw ) then
    write(*,*) 'Not using PAW: early return'
    return
  endif

  ! From potinit
  ! no need to guard agains okpaw here
  CALL PAW_atomic_becsum()
  write(*,*) 'sum rho%bec = ', sum(rho%bec)

  ! Choose atom index and which partial waves to be used (AE or PS)
  ia = 1
  i_what = AE

  i%a = ia   ! atom's index
  i%t = ityp(ia) ! type of atom ia
  i%m = g(i%t)%mesh ! radial mesh size for atom i%t
  i%b = upf(i%t)%nbeta  ! number of beta functions for i%t
  i%l = upf(i%t)%lmax_rho + 1  ! max ang.mom. in augmentation for ia
  l2  = i%l**2
  kkbeta = upf(i%t)%kkbeta
  imesh  = i%m

  IF( .not. upf(i%t)%tpawp ) THEN
    write(*,*) 'Early return'
    return
  endif 

  !
  ! Arrays are allocated inside the cycle to allow reduced
  ! memory usage as different atoms have different meshes
  ALLOCATE( v_lm(i%m,l2,nspin)      )
  ALLOCATE( savedv_lm(i%m,l2,nspin) )
  ALLOCATE( rho_lm(i%m,l2,nspin)    )
  write(*,*) 'size rho_lm = ', shape(rho_lm)

  v_lm(:,:,:) = 0.d0
  savedv_lm(:,:,:) = 0.d0
  rho_lm(:,:,:) = 0.d0

  ! STEP: 1 [ build rho_lm (PAW_rho_lm) ]
  !
  i%ae = i_what
  NULLIFY( rho_core )
  !
  IF( i_what == AE ) THEN
    ! Compute rho spherical harmonics expansion from becsum and pfunc
    CALL PAW_rho_lm( i, rho%bec, upf(i%t)%paw%pfunc, rho_lm )
    write(*,*)
    write(*,*) 'i_what == AE, sum rho_lm = ', sum(rho_lm)
    !
    ! used later for xc potential:
    rho_core => upf(i%t)%paw%ae_rho_atc
    ! sign to sum up the enrgy
    sgn = +1._DP
  ELSE
    !
    CALL PAW_rho_lm( i, rho%bec, upf(i%t)%paw%ptfunc, rho_lm, upf(i%t)%qfuncl )
    !          optional argument for pseudo part (aug. charge) --> ^^^
    write(*,*)
    write(*,*) 'i_what == PS, sum rho_lm = ', sum(rho_lm)
    rho_core => upf(i%t)%rho_atc ! as before
    sgn = -1.0_DP                 ! as before
  ENDIF
             
  ! cleanup auxiliary potentials
  savedv_lm(:,:,:) = 0._DP
   
  !     
  ! Hartree potential stuffs removed
  !
  !................................
             
  !
  ! Then the XC one:
  CALL my_PAW_xc_potential( i, rho_lm, rho_core, v_lm, energy )
  !
  write(*,*) 'X = ', i%a, i_what, sgn*energy
  energy_tot = energy_tot + sgn*energy
  e_cmp(ia, XC, i_what) = sgn*energy
  !
  savedv_lm(:,:,:) = savedv_lm(:,:,:) + v_lm(:,:,:)

  write(*,*)
  write(*,*) '----------------------------'
  write(*,*) 'EXIT debug_PAW_xc_potential'
  write(*,*) '----------------------------'

END SUBROUTINE debug_PAW_xc_potential


!------------------------------------------------------------------
SUBROUTINE my_PAW_xc_potential( i, rho_lm, rho_core, v_lm, energy )
!------------------------------------------------------------------
!  !! Use the density produced by sum_rad_rho to compute xc potential
!  !! and energy, as xc functional is not diagonal on angular momentum
!  !! numerical integration is performed.
!  !
  USE noncollin_module, ONLY : nspin_mag
  use kinds, only : dp
  USE constants, ONLY : e2, eps12
  USE uspp_param, ONLY : upf
  USE lsda_mod, ONLY : nspin
  USE atom, ONLY : g => rgrid
!  USE funct,                  ONLY : dft_is_gradient
  USE xc_lda_lsda, ONLY : xc
  USE paw_variables, ONLY : paw_info, rad
  implicit none
!  !
  TYPE(paw_info), INTENT(IN) :: i !! atom's minimal info
  REAL(DP), INTENT(IN) :: rho_lm(i%m,i%l**2,nspin) !! charge density as lm components
  REAL(DP), INTENT(IN) :: rho_core(i%m) !! core charge, radial and spherical
  REAL(DP), INTENT(OUT) :: v_lm(i%m,i%l**2,nspin) !! potential density as lm components
  REAL(DP), INTENT(OUT) :: energy !! XC energy (if required)
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: rho_loc(:,:)       ! local density (workspace), up and down
  REAL(DP), ALLOCATABLE :: v_rad(:,:,:)    ! radial potential (to be integrated)
  !
  REAL(DP), ALLOCATABLE :: rho_rad(:,:)       ! workspace (only one radial slice of rho)
  !
  REAL(DP), ALLOCATABLE :: e_rad(:)           ! aux, used to store radial slices of energy
  REAL(DP) :: e                               ! aux, used to integrate energy
  !
  INTEGER :: ix,k                             ! counters on directions and radial grid
  INTEGER :: lsd                              ! switch for local spin density
  !
  REAL(DP), ALLOCATABLE :: arho(:,:)
  REAL(DP), ALLOCATABLE :: ex(:), ec(:)
  REAL(DP), ALLOCATABLE :: vx(:,:), vc(:,:)
  REAL(DP), PARAMETER   :: eps = 1.e-30_dp
  !

  ! true if using spin
  lsd = 0
  IF( nspin == 2 )  lsd = 1
  
  ! This will hold the "true" charge density, without r**2 or other factors
  ALLOCATE( rho_loc(i%m,nspin_mag) ) 
  
  rho_loc = 0._DP
  !
  ALLOCATE( rho_rad(i%m,nspin_mag) ) 
  !
  ALLOCATE( arho(i%m,nspin) ) ! XXX: ffr: change second dimension to 2
  ALLOCATE( ex(i%m) )
  ALLOCATE( ec(i%m) )
  ALLOCATE( vx(i%m,2) )
  ALLOCATE( vc(i%m,2) )

  allocate( v_rad(i%m,rad(i%t)%nx,nspin) )

  !
  energy = 0._DP
  ALLOCATE( e_rad(i%m) )
  
  v_rad = 0.0_dp
  DO ix = 1, rad(i%t)%nx
    
    write(*,*)
    write(*,*) 'ix = ', ix

    !
    ! --- LDA (and LSDA) part (no gradient correction) ---
    ! convert _lm density to real density along ix
    !
    CALL my_PAW_lm2rad( i, ix, rho_lm, rho_rad, nspin_mag )
    !
    ! compute the potential along ix
    !
    ! nspin_mag == 4 case is removed
    IF ( nspin == 2 ) THEN
      DO k = 1, i%m
         rho_loc(k,1) = rho_rad(k,1)*g(i%t)%rm2(k)
         rho_loc(k,2) = rho_rad(k,2)*g(i%t)%rm2(k)
      ENDDO
    ELSE
      DO k = 1, i%m
         rho_loc(k,1) = rho_rad(k,1)*g(i%t)%rm2(k)
      ENDDO
    ENDIF

    write(*,*) 'sum rho_loc = ', sum(rho_loc)
    write(*,*) 'sum rho_core = ', sum(rho_core)

    !
    ! Integrate to obtain the energy
    !
    IF( nspin_mag <= 2 ) THEN
      !
      !
      IF( lsd == 0 ) THEN
        !
        arho(:,1) = rho_loc(:,1) + rho_core
        write(*,*) 'sum arho = ', sum(arho)
        !
        CALL xc( i%m, 1, 1, arho(:,1), ex, ec, vx(:,1), vc(:,1) )
        !
        v_rad(:,ix,1) = e2*( vx(:,1) + vc(:,1) )
        e_rad = e2*( ex(:) + ec(:) )
        !
      ELSE
        !
        arho(:,1) = rho_loc(:,1) + rho_loc(:,2) + rho_core(:)
        arho(:,2) = rho_loc(:,1) - rho_loc(:,2)
        !
        CALL xc( i%m, 2, 2, arho, ex, ec, vx, vc )
        !
        v_rad(:,ix,:) = e2*( vx(:,:) + vc(:,:) )
        e_rad(:) = e2*( ex(:) + ec(:) )
        !
      ENDIF
      !
      
      IF (nspin_mag < 2) THEN
        e_rad = e_rad * ( rho_rad(:,1) + rho_core*g(i%t)%r2 )
      ELSEIF (nspin_mag == 2) THEN
        e_rad = e_rad * ( rho_rad(:,1) + rho_rad(:,2) + rho_core*g(i%t)%r2 )
      ENDIF
      !
    ENDIF
    
    write(*,*) 'sum v_rad(:,ix,1) in Ha = ', 0.5d0*sum(v_rad(:,ix,1))
    write(*,*) 'sum e_rad in Ha = ', 0.5d0*sum(e_rad)

    ! Integrate to obtain the energy
    CALL simpson( i%m, e_rad, g(i%t)%rab, e )
    write(*,*) 'integrated energy from ix (in Ha) = ', e*rad(i%t)%ww(ix)*0.5d0
    energy = energy + e * rad(i%t)%ww(ix)
  
  ENDDO

  write(*,*) 'Before: sum v_rad (in Ha) = ', 0.5d0*sum(v_rad)
  write(*,*) 'Before: sum v_lm (in Ha)  = ', 0.5d0*sum(v_lm)

  ! Recompose the sph. harm. expansion
  CALL my_PAW_rad2lm( i, v_rad, v_lm, i%l, nspin_mag )
  
  write(*,*) 'energy (in Ha)    = ', 0.5d0*energy

  write(*,*) 'After: sum v_rad (in Ha) = ', 0.5d0*sum(v_rad)
  write(*,*) 'After: sum v_lm (in Ha)  = ', 0.5d0*sum(v_lm)


  !
  ! Add gradient correction, if necessary
  !
  !IF( dft_is_gradient() ) CALL PAW_gcxc_potential( i, rho_lm, rho_core, v_lm, energy )

  DEALLOCATE( v_rad )
  DEALLOCATE( e_rad )
  !!
  DEALLOCATE( rho_rad )
  DEALLOCATE( rho_loc )
  !!
  DEALLOCATE( arho )
  DEALLOCATE( ex )
  DEALLOCATE( ec )
  DEALLOCATE( vx )
  DEALLOCATE( vc )

  RETURN
  !
END SUBROUTINE


!---------------------------------------------------------------------------------
SUBROUTINE my_PAW_lm2rad( i, ix, F_lm, F_rad, nspin )
!---------------------------------------------------------------------------------
  use kinds, only : dp
  use paw_variables, only : paw_info, rad
  implicit none
  !! Build radial charge distribution from its spherical harmonics expansion.
  !
  TYPE(paw_info), INTENT(IN) :: i
  !! atom's minimal info
  INTEGER :: ix
  !! line of the ylm matrix to use
  !! actually it is one of the nx directions
  INTEGER, INTENT(IN) :: nspin
  !! number of spin components
  REAL(DP), INTENT(IN) :: F_lm(i%m,i%l**2,nspin)
  !! Y_lm expansion of rho
  REAL(DP), INTENT(OUT) :: F_rad(i%m,nspin)
  !! charge density on rad. grid
  !
  ! ... local variables
  !
  INTEGER :: ispin, lm ! counters on angmom and spin
  !
  F_rad(:,:) = 0._DP
  ! cycling on spin is a bit less general...
  DO ispin = 1,nspin
    DO lm = 1, i%l**2
      F_rad(:,ispin) = F_rad(:,ispin) + rad(i%t)%ylm(ix,lm)*F_lm(:,lm,ispin)
    ENDDO ! lm
  ENDDO
  return
  !
END SUBROUTINE

!--------------------------------------------------------------------------------
SUBROUTINE my_PAW_rad2lm( i, F_rad, F_lm, lmax_loc, nspin )
!------------------------------------------------------------------------------
  use kinds, only : dp
  use paw_variables, only : paw_info, rad
  implicit none
  !! Computes:
  !! \[ F_{lm}(r) = \int d \Omega\ F(r,\text{th},\text{ph})\ Y_{lm}(\text{th},
  !! \text{ph}) \]
  !
  TYPE(paw_info), INTENT(IN) :: i
  !! atom's minimal info
  INTEGER, INTENT(IN) :: nspin
  !! spin configuration label
  INTEGER,  INTENT(IN) :: lmax_loc
  !! In some cases I have to keep higher angular components
  !! than the default ones (=lmaxq =the ones present in rho)
  REAL(DP), INTENT(OUT):: F_lm(i%m, lmax_loc**2, nspin)
  !! lm component of F up to lmax_loc
  REAL(DP), INTENT(IN) :: F_rad(i%m, rad(i%t)%nx, nspin)
  !! radial samples of F
  !
  ! ... local variables
  !
  INTEGER :: ix    ! counter for integration
  INTEGER :: lm    ! counter for angmom
  INTEGER :: ispin ! counter for spin
  INTEGER :: j

  write(*,*) 'shape F_lm = ', shape(F_lm)
  write(*,*) 'shape F_rad = ', shape(F_rad)

  write(*,*) 'sum F_lm = ', sum(F_lm)

  DO ispin = 1, nspin
    DO lm = 1, lmax_loc**2
      F_lm(:,lm,ispin) = 0.0_dp
      DO ix = 1, rad(i%t)%nx
        DO j = 1, i%m
           F_lm(j, lm, ispin) = F_lm(j, lm, ispin) + F_rad(j,ix,ispin)* rad(i%t)%wwylm(ix,lm)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  return
END SUBROUTINE
