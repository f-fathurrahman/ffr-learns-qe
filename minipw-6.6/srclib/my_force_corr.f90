!-----------------------------------------------------------------------
SUBROUTINE my_force_corr(forcescc)
!-----------------------------------------------------------------------
  !   This routine calculates the force term vanishing at full
  !     self-consistency. It follows the suggestion of Chan-Bohnen-Ho
  !     (PRB 47, 4771 (1993)). The true charge density is approximated
  !     by means of a free atom superposition.
  !     (alessio f.)
  ! Uses superposition of atomic charges contained in the array rho_at
  ! and read from pseudopotential files
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : msh, rgrid
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, g, ngl, gl, igtongl
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vnew
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions, ONLY : psic
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE 
  !
  REAL(DP) :: forcescc(3,nat)
  !
  REAL(DP), ALLOCATABLE :: rhocgnt(:), aux(:) ! work space
  REAL(DP) ::  gx, arg, fact ! temp factors
  INTEGER :: ir, isup, isdw, ig, nt, na, ndm ! counters
  

  write(*,*)
  write(*,*) 'ENTER my_force_corr'
  write(*,*)


  !
  ! vnew is V_out - V_in, psic is the temp space
  !
  IF( nspin == 1 .or. nspin == 4 ) THEN 
    psic(:) = vnew%of_r(:, 1)
  ELSE 
    isup = 1
    isdw = 2
    psic(:) = ( vnew%of_r(:, isup) + vnew%of_r(:, isdw) ) * 0.5d0
  ENDIF

  write(*,*) 'sum vnew%of_r in Ha = ', 0.5d0*sum(vnew%of_r(:,1))
  write(*,*) 'sum psic before fwfft in Ha = ', 0.5d0*sum(psic)

  CALL fwfft('Rho', psic, dfftp)

  write(*,*) 'sum psic after fwfft (in Ha) = ', sum(psic)*0.5d0

  write(*,*) 'ngm = ', ngm

  IF( gamma_only ) THEN 
    fact = 2.d0
  ELSE 
    fact = 1.d0
  ENDIF

  ndm = MAXVAL( msh(1:ntyp) )
  ALLOCATE( rhocgnt(ngl) )
  ALLOCATE( aux(ndm) )
  !
  DO nt = 1, ntyp
    !
    ! Here we compute the G /= 0 term
    !
    DO ig = gstart, ngl
      gx = sqrt(gl(ig)) * tpiba
      DO ir = 1, msh (nt)
        IF( rgrid(nt)%r(ir) .lt. 1.0d-8 ) THEN 
          aux(ir) = upf(nt)%rho_at(ir)
        ELSE 
          aux(ir) = upf(nt)%rho_at(ir) * sin(gx*rgrid(nt)%r(ir)) / (rgrid(nt)%r(ir)*gx)
        ENDIF 
      ENDDO 
      CALL simpson( msh(nt), aux, rgrid(nt)%rab, rhocgnt(ig) )
    ENDDO 
    !
    write(*,*) 'sum rhocgnt = ', sum(rhocgnt)
    !
    ! sum over atoms
    DO na = 1, nat
      IF( nt == ityp(na) ) THEN 
        forcescc(1:3, na) = 0.0_DP
        ! sum over G-vectors
        DO ig = gstart, ngm
          arg = ( g(1,ig) * tau(1,na) + g(2,ig) * tau(2,na) + g(3,ig) * tau(3,na) ) * tpi
          forcescc(1:3,na) = forcescc(1:3,na) + fact * &
                  rhocgnt(igtongl(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                  g(1:3,ig) * tpiba * CONJG(psic(dfftp%nl(ig)))
        ENDDO 
      ENDIF 
    ENDDO 
    !
  ENDDO

  CALL mp_sum( forcescc, intra_bgrp_comm )
  
  DEALLOCATE( aux )
  DEALLOCATE( rhocgnt )

  write(*,*)
  write(*,*) 'EXIT my_force_corr'
  write(*,*)

  RETURN 

END SUBROUTINE my_force_corr

