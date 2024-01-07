!----------------------------------------------------------------------
SUBROUTINE my_force_cc( forcecc )
!----------------------------------------------------------------------
  !! Calculates the NLCC contribution to the force.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : rgrid, msh
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, g, gg, ngl, gl, igtongl
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions,        ONLY : psic
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP) :: forcecc(3,nat)
  !! output: the NLCC forces on atoms
  !
  ! local variables
  !
  INTEGER :: ig, ir, nt, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms
  REAL(DP), ALLOCATABLE :: vxc(:,:), rhocg(:)
  ! exchange-correlation potential
  ! radial fourier transform of rho core
  REAL(DP) ::  arg, fact

  write(*,*)
  write(*,*) '*** ENTER my_force_cc'
  write(*,*)


  !
  forcecc(:,:) = 0.d0
  !
  IF( ANY( upf(1:ntyp)%nlcc ) ) GOTO 15
  ! early return if there is no pseudopot with nlcc
  RETURN
  !
15 CONTINUE
  IF(gamma_only) THEN
    fact = 2.d0
  ELSE
    fact = 1.d0
  ENDIF
  !
  ! recalculate the exchange-correlation potential (also include rhoe)
  !
  ALLOCATE( vxc(dfftp%nnr,nspin) )
  ! ffr: only for rho_core
  write(*,*) "sum rhoe_core = ", sum(rho_core)
  !
  CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc )
  !
  psic = (0.0_DP,0.0_DP)
  IF (nspin == 1 .OR. nspin == 4) THEN
    DO ir = 1, dfftp%nnr
      psic(ir) = vxc(ir,1)
    ENDDO
  ELSE
     DO ir = 1, dfftp%nnr
        psic (ir) = 0.5d0 * (vxc (ir, 1) + vxc (ir, 2) )
     ENDDO
  ENDIF
  !
  DEALLOCATE( vxc )
  !
  write(*,*) 'sum psic before fwfft (in Ha) = ', 0.5d0*sum(psic)
  CALL fwfft( 'Rho', psic, dfftp )
  write(*,*) 'sum psic after fwfft (in Ha) = ', 0.5d0*sum(psic)
  !
  ! psic contains now Vxc(G)
  !
  ALLOCATE( rhocg(ngl) )
  !
  ! core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  DO nt = 1, ntyp
    IF( upf(nt)%nlcc ) THEN
      !
      CALL drhoc( ngl, gl, omega, tpiba2, msh(nt), rgrid(nt)%r, &
                  rgrid(nt)%rab, upf(nt)%rho_atc, rhocg )
      write(*,*) 'sum rhocg = ', sum(rhocg)
      DO na = 1, nat
        IF (nt == ityp (na) ) THEN
          DO ig = gstart, ngm
            arg = (g(1,ig) * tau(1,na) + g (2, ig) * tau (2, na) &
                 + g(3,ig) * tau(3,na) ) * tpi
            forcecc(1:3,na) = forcecc(1:3, na) + tpiba * omega * &
                    rhocg(igtongl(ig)) * CONJG(psic(dfftp%nl(ig) ) ) * &
                    CMPLX( SIN(arg), COS(arg), KIND=DP) * g(1:3,ig) * fact
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  !
  CALL mp_sum( forcecc, intra_bgrp_comm )
  !
  DEALLOCATE( rhocg )

  write(*,*)
  write(*,*) '*** EXIT my_force_cc'
  write(*,*)

  RETURN
  !
END SUBROUTINE my_force_cc
