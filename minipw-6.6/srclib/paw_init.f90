!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! Modified by ffr

MODULE paw_init
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE

  PUBLIC :: PAW_atomic_becsum
  PUBLIC :: PAW_init_onecenter

  !
  PUBLIC :: allocate_paw_internals, deallocate_paw_internals
  !
  LOGICAL, PARAMETER :: TIMING = .FALSE.
  !
  !
 CONTAINS


!----------------------------------------------------------------------------
SUBROUTINE allocate_paw_internals
!--------------------------------------------------------------------------
  !! Allocate PAW internal variables require for SCF calculation.
  !
  USE lsda_mod,       ONLY : nspin
  USE ions_base,      ONLY : nat
  USE uspp_param,     ONLY : nhm
  !
  USE paw_variables
  !
  IMPLICIT NONE
  write(*,*) 'allocate_paw_internals: nhm = ', nhm
  !
  ALLOCATE( ddd_paw(nhm*(nhm+1)/2,nat,nspin) )
  !
END SUBROUTINE allocate_paw_internals


!------------------------------------------------------------------------------
SUBROUTINE deallocate_paw_internals
!----------------------------------------------------------------------------
  !! Called by \(\textrm{clean_pw}\).
  !
  USE uspp_param,    ONLY : upf
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE paw_variables
  !
  IMPLICIT NONE
  !
  INTEGER :: nt, na
  !
  IF (ALLOCATED(ddd_paw)) DEALLOCATE (ddd_paw)
  !
  IF (ALLOCATED(rad)) THEN
    !
    DO nt = 1,ntyp  
      IF(ASSOCIATED(rad(nt)%ww))       DEALLOCATE( rad(nt)%ww  )  
      IF(ASSOCIATED(rad(nt)%ylm))      DEALLOCATE( rad(nt)%ylm )  
      IF(ASSOCIATED(rad(nt)%wwylm))    DEALLOCATE( rad(nt)%wwylm )  
      IF(ASSOCIATED(rad(nt)%dylmt))    DEALLOCATE( rad(nt)%dylmt )  
      IF(ASSOCIATED(rad(nt)%dylmp))    DEALLOCATE( rad(nt)%dylmp )  
      IF(ASSOCIATED(rad(nt)%cotg_th))  DEALLOCATE( rad(nt)%cotg_th )  
      IF(ASSOCIATED(rad(nt)%cos_phi))  DEALLOCATE( rad(nt)%cos_phi )  
      IF(ASSOCIATED(rad(nt)%sin_phi))  DEALLOCATE( rad(nt)%sin_phi )  
      IF(ASSOCIATED(rad(nt)%cos_th))   DEALLOCATE( rad(nt)%cos_th )  
      IF(ASSOCIATED(rad(nt)%sin_th))   DEALLOCATE( rad(nt)%sin_th )  
    ENDDO
    !
    DEALLOCATE( rad )
    !
  ENDIF
  !
  IF ( ALLOCATED(vs_rad) )   DEALLOCATE( vs_rad )
  !
  paw_is_init = .FALSE.
  !
  RETURN
  !
END SUBROUTINE deallocate_paw_internals


!-----------------------------------------------------------------------------
SUBROUTINE PAW_atomic_becsum()
!--------------------------------------------------------------------------
  !! Initialize becsum with atomic occupations (for PAW atoms only).  
  !! NOTE: requires exact correspondence chi <--> beta in the atom,
  !! that is that all wavefunctions considered for PAW generation are
  !! counted in chi (otherwise the array "oc" does not correspond to beta).
  !
  USE kinds,                ONLY : DP
  USE uspp,                 ONLY : nhtoj, nhtol, indv, becsum
  USE scf,                  ONLY : rho
  USE uspp_param,           ONLY : upf, nh, nhm
  USE ions_base,            ONLY : nat, ityp
  USE lsda_mod,             ONLY : nspin, starting_magnetization
  USE paw_variables,        ONLY : okpaw
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE random_numbers,       ONLY : randy
  USE basis,                ONLY : starting_wfc
  USE noncollin_module,     ONLY : nspin_mag, angle1, angle2
  !
  IMPLICIT NONE
  !
  !REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin)
  INTEGER :: ispin, na, nt, ijh, ih, jh, nb, mb
  REAL(DP) :: noise = 0._DP
  !
  write(*,*)
  write(*,*) '------------------------------------------------------------'
  write(*,*) 'ENTER PAW_atomic_becsum'
  write(*,*) '------------------------------------------------------------'
  !
  !
  IF (.NOT. okpaw) RETURN
  IF (.NOT. ALLOCATED(becsum))   CALL errore( 'PAW_init_becsum', &
                 'Something bad has happened: becsum is not allocated yet', 1 )
  !
  ! Add a bit of random noise if not starting from atomic or saved wfcs:
  IF ( starting_wfc=='atomic+random') noise = 0.05_DP
  IF ( starting_wfc=='random')        noise = 0.10_DP
  !

  noise = 0.d0
  !if(noise > 0.d0) then
  !  write(*,*) 'PAW_atomic_becsum: Random component will be added to becsum'
  !endif

  becsum = 0.0_DP
  na_loop: DO na = 1, nat
    nt = ityp(na)
    is_paw: IF (upf(nt)%tpawp) THEN
      !
      ijh = 1
      ih_loop: DO ih = 1, nh(nt)
        nb = indv(ih,nt)
        !
        IF (nspin == 1) THEN
          !
          becsum(ijh,na,1) = upf(nt)%paw%oc(nb) / DBLE(2*nhtol(ih,nt)+1)
          !
        ELSEIF (nspin == 2) THEN
          !
          becsum(ijh,na,1) = 0.5_dp*(1._DP+starting_magnetization(nt))* &
                             upf(nt)%paw%oc(nb) / DBLE(2*nhtol(ih,nt)+1)
          becsum(ijh,na,2) = 0.5_dp*(1._DP-starting_magnetization(nt))* &
                             upf(nt)%paw%oc(nb) / DBLE(2*nhtol(ih,nt)+1)
          !
        ELSEIF (nspin == 4) THEN
          !
          becsum(ijh,na,1) = upf(nt)%paw%oc(nb)/DBLE(2*nhtol(ih,nt)+1)
          IF (nspin_mag == 4) THEN
            becsum(ijh,na,2) = becsum(ijh,na,1) *              &
                               starting_magnetization(nt)*     &
                               SIN(angle1(nt))*COS(angle2(nt))
            becsum(ijh,na,3) = becsum(ijh,na,1) *              &
                               starting_magnetization(nt)*     &
                               SIN(angle1(nt))*SIN(angle2(nt))
            becsum(ijh,na,4) = becsum(ijh,na,1) *              &
                               starting_magnetization(nt)*     &
                               COS(angle1(nt))
          ENDIF
          !
        ENDIF
        !
        ijh = ijh + 1
        !
        jh_loop: &
         DO jh = ( ih + 1 ), nh(nt)
          !mb = indv(jh,nt)
          DO ispin = 1, nspin_mag
            IF (noise > 0._DP) &
               becsum(ijh,na,ispin) = becsum(ijh,na,ispin) + noise *2._DP*(.5_DP-randy())
          ENDDO
          !
          ijh = ijh + 1
          !
        ENDDO jh_loop
      ENDDO ih_loop
    ENDIF is_paw
  ENDDO na_loop
  !
  ! ... copy becsum in scf structure and symmetrize it
  rho%bec(:,:,:) = becsum(:,:,:)

  write(*,*) 'PAW_atomic_becsum: sum(becsum) before PAW_symmetrize = ', sum(rho%bec)
  write(*,*) 'becsum(1,1,1) before = ', rho%bec(1,1,1)
  write(*,*) 'becsum(2,1,1) before = ', rho%bec(2,1,1)
  if(nspin == 2) then
    write(*,*) 'becsum(1,1,2) before = ', rho%bec(1,1,2)
    write(*,*) 'becsum(2,1,2) before = ', rho%bec(2,1,2)
  endif
  
  ! Write for comparison
  !write(101,*) becsum

  CALL PAW_symmetrize( rho%bec )

  ! For comparison
  !write(102,*) rho%bec

  write(*,*) 'PAW_atomic_becsum: sum(becsum) after PAW_symmetrize = ', sum(rho%bec)
  write(*,*) 'becsum(1,1,1) after = ', rho%bec(1,1,1)
  write(*,*) 'becsum(2,1,1) after = ', rho%bec(2,1,1)
  if(nspin == 2) then
    write(*,*) 'becsum(1,1,2) after = ', rho%bec(1,1,2)
    write(*,*) 'becsum(2,1,2) after = ', rho%bec(2,1,2)
  endif

  write(*,*)
  write(*,*) '------------------------------------------------------------'
  write(*,*) 'EXIT PAW_atomic_becsum'
  write(*,*) '------------------------------------------------------------'


  !
END SUBROUTINE PAW_atomic_becsum


  
!-------------------------------------------------------------------
SUBROUTINE PAW_init_onecenter()
!-----------------------------------------------------------------
  !! This allocates space to store onecenter potential and 
  !! calls PAW_rad_init to initialize onecenter integration.
  !
  USE ions_base,          ONLY : nat, ityp, ntyp => nsp
  USE paw_variables,      ONLY : xlm, lm_fact, lm_fact_x,  &
                                 rad, paw_is_init, vs_rad, &
                                 total_core_energy, only_paw
  USE atom,               ONLY : g => rgrid
  USE radial_grids,       ONLY : do_mesh
  USE uspp_param,         ONLY : upf
  USE lsda_mod,           ONLY : nspin
  USE spin_orb,           ONLY : domag
  USE noncollin_module,   ONLY : noncolin
  USE funct,              ONLY : dft_is_gradient
  USE mp_images,          ONLY : me_image, nproc_image
  USE mp,                 ONLY : mp_sum
  !
  ! ... local variables
  !
  INTEGER :: nt, lmax_safe, lmax_add, ia, ia_s, ia_e, na, mykey, max_mesh, &
             max_nx
  !
  CHARACTER(LEN=12) :: env='            '
  !

  write(*,*)
  write(*,*) '------------------------------------------------------------'
  write(*,*) 'ENTER PAW_init_onecenter'
  write(*,*) '------------------------------------------------------------'

  IF ( paw_is_init ) THEN
    CALL errore( 'PAW_init_onecenter', 'Already initialized!', 1 )
    RETURN
  ENDIF

  !
  ! Init only for the atoms that it will actually use later.
  ! Parallel: divide among processors for the same image
  CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )
  !
  ! Sum all core energies to get...
  total_core_energy = 0._dp
  only_paw = .TRUE.
  max_nx = 0
  max_mesh = 0
  !
  ! Calculate total_core_energy (will be reported at the end of SCF loop)
  DO na = 1, nat
    only_paw = only_paw .AND. upf(ityp(na))%tpawp
    !
    IF( upf(ityp(na))%tpawp ) &
    total_core_energy = total_core_energy + upf(ityp(na))%paw%core_energy
  ENDDO

  write(*,*) 'PAW_init_onecenter: total_core_energy = ', total_core_energy

  !
  ! initialize for integration on angular momentum and gradient, integrating
  ! up to 2*lmaxq (twice the maximum angular momentum of rho) is enough for
  ! H energy and for XC energy. If I have gradient correction I have to go a bit higher
  ALLOCATE( rad(ntyp) )
  !
  DO nt = 1, ntyp
    NULLIFY( rad(nt)%ww    )
    NULLIFY( rad(nt)%ylm   )
    NULLIFY( rad(nt)%wwylm )
    NULLIFY( rad(nt)%dylmt )
    NULLIFY( rad(nt)%dylmp )
    NULLIFY( rad(nt)%cotg_th )
    NULLIFY( rad(nt)%cos_phi )
    NULLIFY( rad(nt)%sin_phi )
    NULLIFY( rad(nt)%cos_th  )
    NULLIFY( rad(nt)%sin_th  )
  ENDDO
  !
  types : &
  DO nt = 1,ntyp
    IF(.NOT. upf(nt)%tpawp) CYCLE types
    ! only allocate radial grid integrator for atomic species
    ! that are actually present on this parallel node:
    DO ia = ia_s, ia_e
      IF (ityp(ia) == nt ) THEN

        WRITE(*,*)
        write(*,*) 'PAW_init_onecenter: nt = ', nt
        write(*,*) 'PAW_init_onecenter: ia = ', ia
        WRITE(*,*) 'PAW_init_onecenter: upf(nt)%lmax_rho = ', upf(nt)%lmax_rho
        
        IF (upf(nt)%lmax_rho == 0) THEN
          ! no need for more than one direction, when it is spherical!
          lmax_safe = 0
          lmax_add  = 0
          write(*,*) 'It is spherical'
        ELSE
          ! 
          IF ( dft_is_gradient() ) THEN
            ! Integrate up to a higher maximum lm if using gradient
            ! correction check expression for d(y_lm)/d\theta for details
            lmax_safe = lm_fact_x*upf(nt)%lmax_rho
            lmax_add  = xlm
          ELSE
            ! no gradient correction:
            lmax_safe = lm_fact*upf(nt)%lmax_rho
            lmax_add  = 0 
          ENDIF
        ENDIF

        write(*,*) 'lmax_safe = ', lmax_safe
        write(*,*) 'lmax_add  = ', lmax_add

        !
        CALL PAW_rad_init( lmax_safe, lmax_add, rad(nt) )
        !
        max_mesh = MAX( max_mesh, g(nt)%mesh )
        max_nx = MAX( max_nx, rad(nt)%nx )

        write(*,*) 'radial grid: rgrid%mesh = ', g(nt)%mesh
        write(*,*) 'rad(nt)%nx = ', rad(nt)%nx
        write(*,*) 'max_mesh = ', max_mesh
        write(*,*) 'max_nx = ', max_nx
        write(*,*) 'upf(nt)%mesh (should be equal to rgrid%mesh) = ', upf(nt)%mesh
        ! ffr: max_mesh and max_nx are used to allocate memory for vs_rad
        ! ffr: vs_rad is used for noncolinear case.

        !write(*,*)
        !write(*,*) 'Radial integrator: '
        !write(*,*) 'nt = ', nt
        !write(*,*) 'rad(nt)%nx = ', rad(nt)%nx
        !write(*,*) 'rad(nt)%ww = ', rad(nt)%ww

        !
        CYCLE types
      ENDIF
    ENDDO

  ENDDO types
  !
  IF (noncolin .AND. domag)  ALLOCATE( vs_rad(max_mesh,max_nx,nat) )
  !
  paw_is_init = .TRUE.

  write(*,*)
  write(*,*) '------------------------------------------------------------'
  write(*,*) 'EXIT PAW_init_onecenter'
  write(*,*) '------------------------------------------------------------'

END SUBROUTINE PAW_init_onecenter

  
!----------------------------------------------------------------------------------
SUBROUTINE PAW_rad_init( l, ls, rad )
!--------------------------------------------------------------------------------
  !! Initialize several quantities related to radial integration: spherical harmonics and their 
  !! gradients along a few (depending on lmaxq) directions, weights for spherical integration.
  !
  ! IMPORTANT: routine PW/summary.f90 has the initialization parameters hardcoded in it
  !            remember to update it if you change this!
  !
  USE constants,              ONLY : pi, fpi, eps8
  USE funct,                  ONLY : dft_is_gradient
  USE paw_variables,          ONLY : paw_radial_integrator
  !
  INTEGER, INTENT(IN) :: l
  !! max angular momentum component that will be integrated
  !! exactly (to numerical precision).
  INTEGER, INTENT(IN) :: ls
  !! additional max l that will be used when computing gradient
  !! and divergence in speherical coords
  TYPE(paw_radial_integrator), INTENT(OUT) :: rad
  !! containt weights and more info to integrate
  !! on radial grid up to lmax = l
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: x(:)    ! nx versors in smart directions
  REAL(DP), ALLOCATABLE :: w(:)    ! temporary integration weights
  REAL(DP), ALLOCATABLE :: r(:,:)  ! integration directions
  REAL(DP), ALLOCATABLE :: r2(:)   ! square modulus of r
  REAL(DP), ALLOCATABLE :: ath(:), aph(:)
                              ! angles in sph coords for r
  INTEGER :: i, ii, n, nphi   ! counters
  INTEGER :: lm, m            ! indexes for ang.mom
  REAL(DP) :: phi, dphi, rho  ! spherical coordinates
  REAL(DP) :: z               ! cartesian coordinates
  ! for gradient corrections:
  INTEGER :: ipol
  REAL(DP), ALLOCATABLE :: aux(:,:)  ! workspace
  REAL(DP) :: vth(3), vph(3)         !versors for theta and phi
  
  write(*,*)
  write(*,*) '--------------------------------------------------------'
  write(*,*) 'Enter PAW_rad_init'
  write(*,*)
  write(*,*) 'l = ', l
  write(*,*) 'ls = ', ls


  IF (TIMING) CALL start_clock( 'PAW_rad_init' )
  !
  
  ! maximum value of l correctly integrated
  rad%lmax = l+ls
  rad%ladd = ls
  
  ! volume element for angle phi
  nphi = rad%lmax + 1 + MOD(rad%lmax,2)
  dphi = 2._DP*pi/nphi !(rad%lmax+1)
  
  ! number of samples for theta angle
  n = (rad%lmax+2)/2
  !
  ALLOCATE( x(n), w(n) )
  
  ! compute weights for theta integration
  CALL gauss_weights( x, w, n )
  
  !
  ! number of integration directions
  rad%nx = n*nphi !(rad%lmax+1)
  !write(*,*) "paw --> directions",rad%nx," lmax:",rad%lmax
  
  !
  ALLOCATE( r(3,rad%nx), r2(rad%nx), rad%ww(rad%nx), ath(rad%nx), aph(rad%nx) )
  !
  ! compute real weights multiplying theta and phi weights
  ii = 0
  !
  DO i = 1, n
    z = x(i)
    rho = SQRT(1._DP-z**2)
    DO m = 1, nphi !rad%lmax
      ii= ii+1
      phi = dphi*DBLE(m-1)
      r(1,ii) = rho*COS(phi)
      r(2,ii) = rho*SIN(phi)
      r(3,ii) = z
      rad%ww(ii) = w(i)*2._dp*pi/nphi !(rad%lmax+1)
      r2(ii) = r(1,ii)**2 + r(2,ii)**2 + r(3,ii)**2
      ! these will be used later:
      ath(ii) = ACOS(z/SQRT(r2(ii)))
      aph(ii) = phi
    ENDDO
  ENDDO

  ! cleanup
  DEALLOCATE( x, w )
  !
  ! initialize spherical harmonics that will be used
  ! to convert rho_lm to radial grid
  rad%lm_max = (rad%lmax+1)**2
  !
  ALLOCATE( rad%ylm(rad%nx,rad%lm_max) )
  CALL ylmr2( rad%lm_max, rad%nx, r, r2, rad%ylm )
  ! As I will mostly use the product ww*ylm I can 
  ! precompute it here:
  ALLOCATE( rad%wwylm(rad%nx,rad%lm_max) )
  !
  DO i = 1, rad%nx
    DO lm = 1, rad%lm_max
      rad%wwylm(i,lm) = rad%ww(i) * rad%ylm(i,lm)
    ENDDO
  ENDDO

  !
  ALLOCATE( rad%cos_phi(rad%nx) )
  ALLOCATE( rad%sin_phi(rad%nx) )
  ALLOCATE( rad%cos_th(rad%nx) )
  ALLOCATE( rad%sin_th(rad%nx) )
  !
  DO i = 1, rad%nx
    rad%cos_phi(i) = COS(aph(i))
    rad%sin_phi(i) = SIN(aph(i))
    rad%cos_th(i) = COS(ath(i))
    rad%sin_th(i) = SIN(ath(i))
  ENDDO
  !
  ! if gradient corrections will be used than we need
  ! to initialize the gradient of ylm, as we are working in spherical
  ! coordinates the formula involves \hat{theta} and \hat{phi}
  gradient: IF (dft_is_gradient()) THEN

    write(*,*) 'Using gradient correction in PAW_rad_init'
    write(*,*) 'rad%lm_max = ', rad%lm_max

    ALLOCATE( rad%dylmt(rad%nx,rad%lm_max), &
              rad%dylmp(rad%nx,rad%lm_max), &
              aux(rad%nx,rad%lm_max) )
    ALLOCATE( rad%cotg_th(rad%nx) )
    !
    rad%dylmt(:,:) = 0._DP
    rad%dylmp(:,:) = 0._DP
    !
    ! Compute derivative along x, y and z => gradient, then compute the
    ! scalar products with \hat{theta} and \hat{phi} and store them in
    ! dylmt and dylmp respectively.
    !
    DO ipol = 1, 3 !x,y,z
      !
      CALL dylmr2( rad%lm_max, rad%nx, r,r2, aux, ipol )
      !
      DO lm = 1, rad%lm_max
        DO i = 1, rad%nx
          vph = (/-SIN(aph(i)), COS(aph(i)), 0._DP/)
          ! this is the explicit form, but the cross product trick (below) is much faster:
          ! vth = (/COS(aph(i))*COS(ath(i)), SIN(aph(i))*COS(ath(i)), -SIN(ath(i))/)
          vth = (/vph(2)*r(3,i)-vph(3)*r(2,i),&
                  vph(3)*r(1,i)-vph(1)*r(3,i),&
                  vph(1)*r(2,i)-vph(2)*r(1,i)/)
          rad%dylmt(i,lm) = rad%dylmt(i,lm) + aux(i,lm)*vth(ipol)
          ! CHECK: the 1/SIN(th) factor should be correct, but deals wrong result, why?
          rad%dylmp(i,lm) = rad%dylmp(i,lm) + aux(i,lm)*vph(ipol) !/SIN(ath(i))
        ENDDO
      ENDDO
      !
    ENDDO
    !
    DO i = 1, rad%nx
       rad%cotg_th(i) = COS(ath(i))/SIN(ath(i))
    ENDDO
    !
    DEALLOCATE( aux )
    !
  ENDIF gradient
  ! cleanup
  DEALLOCATE( r, r2, ath, aph )
  !
  IF (TIMING) CALL stop_clock( 'PAW_rad_init' )
  

  write(*,*)
  write(*,*) 'EXIT PAW_rad_init'
  write(*,*) '--------------------------------------------------------'


  !
CONTAINS

  !---------------------------------------------------
  SUBROUTINE gauss_weights( x, w, n )
  !--------------------------------------------------
    !! Computes weights for gaussian integrals,
    !! from numerical recipes.
    !
    USE constants, ONLY : pi, eps => eps12
    !
    IMPLICIT NONE
    !
    INTEGER :: n, i, j, m
    REAL(8) :: x(n), w(n), z, z1, p1, p2, p3, pp
    !
    m = (n+1)/2
    !
    DO i = 1, m
      z1 = 2._DP
      z = COS(pi*(i-0.25_DP)/(n+0.5_DP))
      DO WHILE( ABS(z-z1) > eps )
        p1 = 1._DP
        p2 = 0._DP
        DO j = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2._DP*j-1._DP)*z*p2-(j-1._DP)*p3)/j
        ENDDO
        pp = n*(z*p1-p2)/(z*z-1._DP)
        z1 = z
        z  = z1-p1/pp
      ENDDO
      x(i) = -z
      x(n+1-i) = z
      w(i) = 2._DP/((1._DP-z*z)*pp*pp)
      w(n+1-i) = w(i)
    ENDDO
    !
  END SUBROUTINE gauss_weights ! internal subroutine in PAW_rad_init
    !
END SUBROUTINE PAW_rad_init 


END MODULE paw_init
