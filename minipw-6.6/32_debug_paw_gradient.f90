INCLUDE 'prepare_all.f90'


PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL debug_PAW_gcxc_potential()
END PROGRAM



!----------------------------------
SUBROUTINE debug_PAW_gcxc_potential()
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
  write(*,*) '------------------------------'
  write(*,*) 'ENTER debug_PAW_gcxc_potential'
  write(*,*) '------------------------------'

  if( .not. okpaw ) then
    write(*,*) 'Not using PAW: early return'
    return
  endif

  ! From potinit
  ! no need to guard agains okpaw here
  CALL PAW_atomic_becsum()
  write(*,*) 'sum rho%bec = ', sum(rho%bec)

  ! Choose atom index and which partial waves to be used (AE or PS)
  ia = 2
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
  ! Then the GCXC one:
  CALL my_PAW_gcxc_potential( i, rho_lm, rho_core, v_lm, energy )
  

  !
  write(*,*) 'X = ', i%a, i_what, sgn*energy
  energy_tot = energy_tot + sgn*energy
  e_cmp(ia, XC, i_what) = sgn*energy
  !
  savedv_lm(:,:,:) = savedv_lm(:,:,:) + v_lm(:,:,:)

  write(*,*)
  write(*,*) '-----------------------------'
  write(*,*) 'EXIT debug_PAW_gcxc_potential'
  write(*,*) '-----------------------------'

END SUBROUTINE debug_PAW_gcxc_potential




!------------------------------------------------------------------------------------
SUBROUTINE my_PAW_gcxc_potential(i, rho_lm, rho_core, v_lm, energy)
!---------------------------------------------------------------------------------
  use kinds, only : DP
  USE noncollin_module, ONLY : nspin_mag
  USE constants, ONLY : e2, eps12
  USE uspp_param, ONLY : upf
  USE xc_lda_lsda, ONLY : xc
  USE paw_variables, ONLY : paw_info, rad

  use paw_onecenter, only : PAW_lm2rad, PAW_rad2lm3, &
    PAW_divergence, PAW_gradient, PAW_rad2lm, compute_rho_spin_lm

  !! Add gradient correction to v_xc, code mostly adapted from ../atomic/vxcgc.f90
  !! in order to support non-spherical charges (as Y_lm expansion).  
  !! Note that the first derivative in vxcgc becomes a gradient, while the second is
  !! a divergence.  
  !! We also have to temporarily store some additional Y_lm components in order not
  !! to loose precision during the calculation, even if only the ones up to 
  !! lmax_rho (the maximum in the density of charge) matter when computing \int v*rho.
  !
  USE lsda_mod,               ONLY : nspin
  USE noncollin_module,       ONLY : noncolin, nspin_mag, nspin_gga
  USE atom,                   ONLY : g => rgrid
  USE constants,              ONLY : sqrtpi, fpi,pi,e2
  USE funct,                  ONLY : igcc_is_lyp
  USE xc_gga,                 ONLY : xc_gcx
  USE mp,                     ONLY : mp_sum

  implicit none


  !
  ! These variables are originally the arguments
  TYPE(paw_info) :: i !! atom's minimal info (input)
  REAL(DP), intent(in) :: rho_lm(i%m,i%l**2,nspin) !! charge density as lm components
  REAL(DP), intent(in) :: rho_core(i%m) !! core charge, radial and spherical
  REAL(DP), intent(inout) :: v_lm(i%m,i%l**2,nspin) !! potential to be updated
  REAL(DP), intent(inout) :: energy  !! if present, add GC to energy
  !
  ! local variables
  !
  REAL(DP), PARAMETER :: epsr = 1.e-6_DP, epsg = 1.e-10_DP
  ! (as in PW/src/gradcorr.f90)
  !
  REAL(DP), ALLOCATABLE :: rho_rad(:,:) ! charge density sampled
  REAL(DP), ALLOCATABLE :: grad(:,:,:)  ! gradient
  REAL(DP), ALLOCATABLE :: gradx(:,:,:) ! gradient (swapped indexes)
  REAL(DP), ALLOCATABLE :: grad2(:,:)   ! square modulus of gradient
                                        ! (first of charge, than of hamiltonian)
  REAL(DP), ALLOCATABLE :: gc_rad(:,:,:)    ! GC correction to V (radial samples)
  REAL(DP), ALLOCATABLE :: gc_lm(:,:,:)     ! GC correction to V (Y_lm expansion)
  REAL(DP), ALLOCATABLE :: h_rad(:,:,:,:)   ! hamiltonian (vector field)
  REAL(DP), ALLOCATABLE :: h_lm(:,:,:,:)    ! hamiltonian (vector field)
                                  ! ^^^^^^^^^^^^^^^^^^ expanded to higher lm than rho !
  REAL(DP), ALLOCATABLE :: div_h(:,:,:)  ! div(hamiltonian)
  !
  REAL(DP), ALLOCATABLE :: rhoout_lm(:,:,:) ! charge density as lm components
  REAL(DP), ALLOCATABLE :: vout_lm(:,:,:)   ! potential as lm components
  REAL(DP), ALLOCATABLE :: segni_rad(:,:)   ! sign of the magnetization
  !
  REAL(DP), ALLOCATABLE :: arho(:,:), grad2_v(:)
  REAL(DP), ALLOCATABLE :: r_vec(:,:)
  !
  REAL(DP), DIMENSION(i%m,nspin_gga) :: v1x, v2x, v1c, v2c  !workspace
  REAL(DP), DIMENSION(i%m) :: sx, sc
  REAL(DP), ALLOCATABLE :: v2cud(:)
  !
  REAL(DP) :: vnull
  !
  REAL(DP), ALLOCATABLE :: e_rad(:)      ! aux, used to store energy
  REAL(DP) :: e, e_gcxc                  ! aux, used to integrate energy
  !
  INTEGER  :: k, ix, is, lm              ! counters on spin and mesh
  REAL(DP) :: sgn                        ! workspace
  REAL(DP) :: co2                        ! workspace
  !
  INTEGER :: mytid, ntids

  REAL(DP), ALLOCATABLE :: egcxc_of_tid(:)

  
  e_gcxc = 0.0d0

  ALLOCATE( gc_rad(i%m,rad(i%t)%nx,nspin_gga) )! GC correction to V (radial samples)
  ALLOCATE( gc_lm(i%m,i%l**2,nspin_gga)       )! GC correction to V (Y_lm expansion)
  ALLOCATE( h_rad(i%m,3,rad(i%t)%nx,nspin_gga))! hamiltonian (vector field)
  ALLOCATE( h_lm(i%m,3,(i%l+rad(i%t)%ladd)**2,nspin_gga) ) 
                                        ! ^^^^^^^^^^^^^^^^^^ expanded to higher lm than rho !
  ALLOCATE(div_h(i%m,i%l**2,nspin_gga))
  ALLOCATE(rhoout_lm(i%m,i%l**2,nspin_gga)) ! charge density as lm components
  ALLOCATE(vout_lm(i%m,i%l**2,nspin_gga))   ! potential as lm components
  ALLOCATE(segni_rad(i%m,rad(i%t)%nx))      ! charge density as lm components
  
  vout_lm=0.0_DP
  
  !
  IF ( nspin_mag == 2 .OR. nspin_mag == 4 ) THEN
    ! transform the noncollinear case into sigma-GGA case
    IF (noncolin) THEN
      CALL compute_rho_spin_lm(i, rho_lm, rhoout_lm, segni_rad)
    ELSE
      rhoout_lm = rho_lm
    ENDIF
  ENDIF

  mytid = 1
  ntids = 1

  ALLOCATE( rho_rad(i%m,nspin_gga)) ! charge density sampled
  ALLOCATE( grad(i%m,3,nspin_gga) ) ! gradient
  ALLOCATE( grad2(i%m,nspin_gga)  ) ! square modulus of gradient
                                      ! (first of charge, than of hamiltonian)
  gc_rad = 0.0d0
  h_rad  = 0.0d0
  !
  
  !IF (PRESENT(energy)) THEN
  ALLOCATE( egcxc_of_tid(ntids) )
  egcxc_of_tid(mytid) = 0.0_dp
  ALLOCATE( e_rad(i%m) )
  !ENDIF
  
  !
  spin:&
  !
  IF( nspin_mag == 1 ) THEN
    !
    ! GGA case
    !
    ALLOCATE( arho(i%m,1), grad2_v(i%m) )
    ALLOCATE( gradx(3,i%m,1) )
    !
    DO ix = 1, rad(i%t)%nx
      !
      !  WARNING: the next 2 calls are duplicated for spin==2
      CALL PAW_lm2rad( i, ix, rho_lm, rho_rad, nspin_mag )
      CALL PAW_gradient( i, ix, rho_lm, rho_rad, rho_core, grad2, grad )
      !
      DO k = 1, i%m
        arho(k,1) = rho_rad(k,1)*g(i%t)%rm2(k) + rho_core(k)
        arho(k,1) = ABS(arho(k,1))
        gradx(:,k,1) = grad(k,:,1)
      ENDDO
      !
      CALL xc_gcx( i%m, 1, arho, gradx, sx, sc, v1x, v2x, v1c, v2c )
      !
      DO k = 1, i%m
        e_rad(k) = e2 * (sx(k)+sc(k)) * g(i%t)%r2(k)
        gc_rad(k,ix,1)  = (v1x(k,1)+v1c(k,1))  !*g(i%t)%rm2(k)
        h_rad(k,:,ix,1) = (v2x(k,1)+v2c(k,1))*grad(k,:,1)*g(i%t)%r2(k)
      ENDDO
      !
      ! integrate energy (if required)
      !IF( PRESENT(energy) ) THEN
         CALL simpson(i%m, e_rad, g(i%t)%rab, e)
         egcxc_of_tid(mytid) = egcxc_of_tid(mytid) + e*rad(i%t)%ww(ix)
      !ENDIF
      !
    ENDDO
    !
    DEALLOCATE( arho, grad2_v ) 
    DEALLOCATE( gradx )
    !
    !
    ELSEIF ( nspin_mag == 2 .OR. nspin_mag == 4 ) THEN
        !
        ALLOCATE( gradx(3,i%m,2) )
        ALLOCATE( r_vec(i%m,2) )
        ALLOCATE( v2cud(i%m) )
        !
        !   this is the \sigma-GGA case
        !
        DO ix = 1, rad(i%t)%nx
           !
           CALL PAW_lm2rad( i, ix, rhoout_lm, rho_rad, nspin_gga )
           CALL PAW_gradient( i, ix, rhoout_lm, rho_rad, rho_core,grad2, grad )
           !
           DO k = 1, i%m
               !
               ! Prepare the necessary quantities
               ! rho_core is considered half spin up and half spin down:
               co2 = rho_core(k)/2
               ! than I build the real charge dividing by r**2
               r_vec(k,1) = rho_rad(k,1)*g(i%t)%rm2(k) + co2
               r_vec(k,2) = rho_rad(k,2)*g(i%t)%rm2(k) + co2
               !
               !
               gradx(:,k,1) = grad(k,:,1)
               gradx(:,k,2) = grad(k,:,2)
           ENDDO
           !
           CALL xc_gcx( i%m, 2, r_vec, gradx, sx, sc, v1x, v2x, v1c, v2c, v2cud )
           !
           DO k = 1, i%m
              !
              e_rad(k) = e2*(sx(k)+sc(k))*g(i%t)%r2(k)
              !
              ! first term of the gradient correction : D(rho*Exc)/D(rho)
              gc_rad(k,ix,1)  = (v1x(k,1)+v1c(k,1)) !*g(i%t)%rm2(k)
              gc_rad(k,ix,2)  = (v1x(k,2)+v1c(k,2)) !*g(i%t)%rm2(k)
              !
              ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
              ! h_rad(k,:,ix,1) =( (v2xup_vec(k)+v2c)*grad(k,:,1)+v2c*grad(k,:,2) )*g(i%t)%r2(k)
              ! h_rad(k,:,ix,2) =( (v2xdw_vec(k)+v2c)*grad(k,:,2)+v2c*grad(k,:,1) )*g(i%t)%r2(k)
              h_rad(k,:,ix,1) =( (v2x(k,1)+v2c(k,1))*grad(k,:,1) + &
                                  v2cud(k)*grad(k,:,2) )*g(i%t)%r2(k)
              h_rad(k,:,ix,2) =( (v2x(k,2)+v2c(k,2))*grad(k,:,2) + &
                                  v2cud(k)*grad(k,:,1) )*g(i%t)%r2(k)
              !
           ENDDO
           !
           ! integrate energy (if required)
           ! NOTE: this integration is duplicated for every spin, FIXME!
           !IF (PRESENT(energy)) THEN
               CALL simpson( i%m, e_rad, g(i%t)%rab, e )
               egcxc_of_tid(mytid) = egcxc_of_tid(mytid) + e * rad(i%t)%ww(ix)
           !ENDIF
           !
        ENDDO ! ix
        !
        DEALLOCATE( gradx )
        DEALLOCATE( r_vec )
        DEALLOCATE( v2cud )
        !
    ELSE spin
    !
        CALL errore( 'PAW_gcxc_v', 'unknown spin number', 2 )
    ENDIF spin
    
    DEALLOCATE( e_rad )
    !
    DEALLOCATE( rho_rad )
    DEALLOCATE( grad  )
    DEALLOCATE( grad2 )
    energy = energy + e_gcxc

    DEALLOCATE( egcxc_of_tid )
    !
    ! convert the first part of the GC correction back to spherical harmonics
    CALL PAW_rad2lm( i, gc_rad, gc_lm, i%l, nspin_gga )
    !
    ! Note that the expansion into spherical harmonics of the derivative 
    ! with respect to theta of the spherical harmonics, is very slow to
    ! converge and would require a huge angular momentum ladd.
    ! This derivative divided by sin_th is much faster to converge, so
    ! we divide here before calculating h_lm and keep into account for
    ! this factor sin_th in the expression of the divergence.
    !
    ! ADC 30/04/2009.
    ! 
    DO ix = 1, rad(i%t)%nx
       h_rad(1:i%m,3,ix,1:nspin_gga) = h_rad(1:i%m,3,ix,1:nspin_gga) / &
                                       rad(i%t)%sin_th(ix)
    ENDDO
    ! We need the gradient of H to calculate the last part of the exchange
    ! and correlation potential. First we have to convert H to its Y_lm expansion
    CALL PAW_rad2lm3( i, h_rad, h_lm, i%l+rad(i%t)%ladd, nspin_gga )
    !
    ! Compute div(H)
    CALL PAW_divergence( i, h_lm, div_h, i%l+rad(i%t)%ladd, i%l )
    !                       input max lm --^  output max lm-^
    !
    ! Finally sum it back into v_xc
    DO is = 1,nspin_gga
      DO lm = 1,i%l**2
         vout_lm(1:i%m,lm,is) = vout_lm(1:i%m,lm,is) + e2*(gc_lm(1:i%m,lm,is)-div_h(1:i%m,lm,is))
      ENDDO
    ENDDO
  !
  IF (nspin_mag == 4 ) then
    stop 'nspin_mag == 4 is disabled'
  endif

  v_lm(:,:,1:nspin_mag) = v_lm(:,:,1:nspin_mag)+vout_lm(:,:,1:nspin_mag)
  
  !
  DEALLOCATE( gc_rad )
  DEALLOCATE( gc_lm  )
  DEALLOCATE( h_rad  )
  DEALLOCATE( h_lm   )
  DEALLOCATE( div_h  )
  DEALLOCATE( rhoout_lm )
  DEALLOCATE( vout_lm   )
  DEALLOCATE( segni_rad )
  !
  ! if(PRESENT(energy)) write(*,*) "gcxc -->", e_gcxc
  !
END SUBROUTINE