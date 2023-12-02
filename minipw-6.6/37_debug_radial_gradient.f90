INCLUDE 'prepare_all.f90'

PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL debug_radial_gradient()
END PROGRAM

INCLUDE 'my_radial_gradient.f90'

!---------------------------------
subroutine debug_radial_gradient()
!---------------------------------

  USE kinds,          ONLY : DP
  USE paw_variables,  ONLY : paw_info
  USE constants,         ONLY : sqrtpi
  USE paw_onecenter,     ONLY : paw_rho_lm
  USE atom,              ONLY : g => rgrid
  USE ions_base,         ONLY : ityp
  USE lsda_mod,          ONLY : nspin
  USE uspp_param,        ONLY : upf
  USE scf,               ONLY : scf_type
  USE fft_types,         ONLY : fft_index_to_3d
  USE splinelib,         ONLY : spline, splint
  USE scf, ONLY : rho

  implicit none

  logical :: withcore
  TYPE(paw_info) :: i ! minimal info on atoms
  INTEGER :: ir   ! counter on grid point
  INTEGER :: is   ! spin index
  INTEGER :: lm   ! counters on angmom and radial grid
  INTEGER :: ia
  REAL(DP), ALLOCATABLE :: d1y(:), d2y(:)
  REAL(DP), ALLOCATABLE :: rho_lm(:,:,:), rho_lm_ae(:,:,:), rho_lm_ps(:,:,:)
  REAL(DP) :: first, second

  ! choose atom index
  ia = 1
  withcore = .false.

  !
  i%a = ia                      ! atom's index
  i%t = ityp(ia)                ! type of atom ia
  i%m = g(i%t)%mesh             ! radial mesh size for atom i%t
  i%b = upf(i%t)%nbeta          ! number of beta functions for i%t
  i%l = upf(i%t)%lmax_rho+1 ! max ang.mom. in augmentation for ia

  IF( .NOT. upf(i%t)%tpawp ) THEN
    write(*,*) 'atom ', ia, ' is not using PAW: early return'
    return
  ENDIF
  
  ALLOCATE( rho_lm_ae(i%m,i%l**2,nspin) )
  allocate( rho_lm_ps(i%m,i%l**2,nspin) )
  ALLOCATE( rho_lm(i%m,i%l**2,nspin) )

  write(*,*) 'sum rho%bec = ', sum(rho%bec)

  !
  ! Compute rho spherical harmonics expansion from becsum and pfunc
  CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%pfunc,  rho_lm_ae)
  CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%ptfunc, rho_lm_ps, upf(i%t)%qfuncl)
  !
  DO is=1,nspin
    DO lm = 1,i%l**2
      DO ir = 1, i%m
        rho_lm(ir,lm,is) = ( rho_lm_ae(ir,lm,is) - rho_lm_ps(ir,lm,is) ) * g(i%t)%rm2(ir)
      ENDDO
    ENDDO
    !
    ! add core charge
    !
    IF (withcore) THEN
      DO ir = 1, i%m
        rho_lm(ir,1,is) = rho_lm(ir,1,is) + upf(i%t)%paw%ae_rho_atc(ir) / nspin * (2._DP * sqrtpi)
      ENDDO
    ENDIF
  ENDDO
  write(*,*) 'sum rho_lm = ', sum(rho_lm_ae)
  write(*,*) 'sum rho_ps = ', sum(rho_lm_ps)
  write(*,*) 'sum rho_lm = ', sum(rho_lm)
  ! deallocate asap
  DEALLOCATE(rho_lm_ae, rho_lm_ps)
  !
  !
  ALLOCATE( d1y(upf(i%t)%kkbeta), d2y(upf(i%t )%kkbeta) )
  
  is = 1 ! spin index
  lm = 1 ! max i%l**2

  write(*,*) 'kkbeta = ', upf(i%t)%kkbeta
  write(*,*) 'Nrmesh = ', i%m
  
  ! last argument = 1 means using coarse grid algorithm
  CALL my_radial_gradient(rho_lm(1:upf(i%t)%kkbeta,lm,is), d1y, &
                        g(i%t)%r, upf(i%t)%kkbeta, 1)
  CALL my_radial_gradient(d1y, d2y, g(i%t)%r, upf(i%t)%kkbeta, 1)
  !
  first  = d1y(1) ! first derivative in first point
  second = d2y(1) ! second derivative in first point

  write(*,*) 'first = ', first
  write(*,*) 'second = ', second
  return

end subroutine


