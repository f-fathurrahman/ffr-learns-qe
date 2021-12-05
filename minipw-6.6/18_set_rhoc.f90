INCLUDE 'prepare_all.f90'

PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()

  call my_set_rhoc()


END PROGRAM


! Input: upf, msh, rgrid
! Output: rho_core, rhog_core
!-----------------------------------------------------------------------
SUBROUTINE my_set_rhoc()
!-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : dp
  USE io_global, ONLY : stdout
  USE atom,      ONLY : msh, rgrid
  USE uspp_param,ONLY : upf
  USE ions_base, ONLY : ntyp => nsp
  USE cell_base, ONLY : omega, tpiba2
  USE fft_base,  ONLY : dfftp
  USE fft_rho,   ONLY : rho_g2r
  USE gvect,     ONLY : ngm, ngl, gl, igtongl
  USE vlocal,    ONLY : strf
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE scf,       ONLY : rho_core, rhog_core
  !
  IMPLICIT NONE
  !
  REAL(DP) , ALLOCATABLE ::  rhocg(:)
  ! the radial fourier transform
  REAL(DP) ::  rhoneg
  ! used to check the core charge
  INTEGER :: ir, nt, ng
  ! counter on mesh points
  ! counter on atomic types
  ! counter on g vectors

  rhog_core(:) = 0.0_DP
  rho_core(:)  = 0.0_DP

  IF( ANY(upf(1:ntyp)%nlcc) ) THEN
    !
    ALLOCATE(rhocg(ngl))    
    !
    ! the sum is on atom types
    !
    DO nt = 1, ntyp
      IF( upf(nt)%nlcc ) THEN
        write(*,*) 'NLCC is used for species: ', nt
        !
        ! drhoc computes the radial fourier transform for each shell of g vec
        !
        CALL my_drhoc(ngl, gl, omega, tpiba2, msh(nt), rgrid(nt)%r, &
          rgrid(nt)%rab, upf(nt)%rho_atc, rhocg)
        !
        !     multiply by the structure factor and sum
        !
        DO ng = 1, ngm
          rhog_core(ng) = rhog_core(ng) + strf(ng,nt) * rhocg(igtongl(ng))
        ENDDO
      ENDIF
    ENDDO
    DEALLOCATE(rhocg)
    !
    CALL rho_g2r( dfftp, rhog_core, rho_core )
    !
    ! test on the charge and computation of the core energy
    !
    rhoneg = 0.d0
    DO ir = 1, dfftp%nnr
      rhoneg = rhoneg + min(0.d0, rho_core(ir) )
      ! XXX: Do not et core charge to be positive definite
    ENDDO
    !
    rhoneg = rhoneg / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
    CALL mp_sum( rhoneg, intra_bgrp_comm )
    !

    IF(rhoneg < -1.0d-6) WRITE(stdout,'(/5x,"Check: negative core charge=",2f12.6)') rhoneg

    ! calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
    ! The term was present in previous versions of the code but it shouldn't
  ENDIF
  
  write(*,*) 'integ rho_core = ', sum(rho_core)*omega/dfftp%nnr
  write(*,*) 'sum rhoe_core = ', sum(rho_core)

  do ir = 1,6
    write(*,'(1x,I8,F18.10)') ir, rho_core(ir)
  enddo

  !do ir = 1,dfftp%nnr
  !  write(*,'(1x,I8,F18.10)') ir, rho_core(ir)
  !enddo

  !do ir = dfftp%nnr-5,dfftp%nnr
  !  write(*,'(1x,I8,F18.10)') ir, rho_core(ir)
  !enddo


  RETURN

END SUBROUTINE


!--------------------------------------------------------------------
SUBROUTINE my_drhoc( ngl, gl, omega, tpiba2, mesh, r, rab, rhoc, rhocg )
!--------------------------------------------------------------------
  !! Calculates the Fourier transform of the core charge.
  !
  USE kinds
  USE constants, ONLY: pi, fpi
  !
  IMPLICIT NONE
  !
  INTEGER :: ngl !! input: the number of g shell
  INTEGER :: mesh !! input: the number of radial mesh points
  REAL(DP) :: gl(ngl) !! input: the number of G shells
  REAL(DP) :: r(mesh) !! input: the radial mesh
  REAL(DP) :: rab(mesh) !! input: the derivative of the radial mesh
  REAL(DP) :: rhoc(mesh) !! input: the radial core charge
  REAL(DP) :: omega !! input: the volume of the unit cell
  REAL(DP) :: tpiba2 !! input: 2 times pi / alat
  REAL(DP) :: rhocg(ngl) !! output: the Fourier transform of the core charge
  !
  ! ... local variables
  !
  REAL(DP) :: gx, rhocg1
  ! the modulus of g for a given shell
  ! the Fourier transform
  REAL(DP), ALLOCATABLE :: aux(:)
  ! auxiliary memory for integration
  INTEGER :: ir, igl, igl0
  ! counter on radial mesh points
  ! counter on g shells
  ! lower limit for loop on ngl
  !

  ALLOCATE( aux(mesh) )
  !
  ! G=0 term
  !
  IF( gl(1) < 1.0d-8 ) THEN
    DO ir = 1, mesh
      aux(ir) = r(ir)**2 * rhoc(ir)
    ENDDO
    CALL simpson( mesh, aux, rab, rhocg1 )
    rhocg(1) = fpi * rhocg1 / omega
    igl0 = 2
  ELSE
    igl0 = 1
  ENDIF
  !
  ! G != 0 term
  !
  DO igl = igl0, ngl
    gx = SQRT(gl(igl) * tpiba2)
    CALL sph_bes( mesh, r, gx, 0, aux )
    DO ir = 1, mesh
      aux(ir) = r(ir)**2 * rhoc(ir) * aux(ir)
    ENDDO
    CALL simpson( mesh, aux, rab, rhocg1 )
    rhocg(igl) = fpi * rhocg1 / omega
  ENDDO
  DEALLOCATE( aux )

  RETURN
END SUBROUTINE 

