!-----------------------------------------------------------------------
SUBROUTINE my_init_at_1()
!-----------------------------------------------------------------------
  !! This routine computes a table with the radial Fourier transform 
  !! of the atomic wavefunctions.
  !
  USE kinds,        ONLY : DP
  USE atom,         ONLY : rgrid, msh
  USE constants,    ONLY : fpi
  USE cell_base,    ONLY : omega
  USE ions_base,    ONLY : ntyp => nsp
  USE us,           ONLY : tab_at, nqx, dq
  USE uspp_param,   ONLY : upf
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER :: nt, nb, iq, ir, l, startq, lastq, ndm
  !
  REAL(DP), ALLOCATABLE :: aux(:), vchi(:)
  REAL(DP) :: vqint, pref, q
  !
  ndm = MAXVAL(msh(1:ntyp))
  
  write(*,*)
  write(*,*) '<div> ENTER my_init_at_1'
  write(*,*)
  
  ALLOCATE( aux(ndm), vchi(ndm) )
  
  write(*,*) 'ndm = ', ndm
  write(*,*) 'nqx = ', nqx
  write(*,*) 'dq = ', dq
  
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  pref = fpi / SQRT(omega)
  ! needed to normalize atomic wfcs (not a bad idea in general and 
  ! necessary to compute correctly lda+U projections)
  CALL divide( intra_bgrp_comm, nqx, startq, lastq )
  !
  tab_at(:,:,:) = 0.0_DP
  !
  DO nt = 1, ntyp
    DO nb = 1, upf(nt)%nwfc
      !
      IF (upf(nt)%oc(nb) >= 0.0_DP) THEN
        l = upf(nt)%lchi (nb)
        !
        DO iq = startq, lastq
          q = dq * (iq - 1)
          CALL sph_bes( msh(nt), rgrid(nt)%r, q, l, aux )
          DO ir = 1, msh(nt)
            vchi(ir) = upf(nt)%chi(ir,nb) * aux(ir) * rgrid(nt)%r(ir)
          ENDDO
          CALL simpson( msh(nt), vchi, rgrid(nt)%rab, vqint )
          tab_at( iq, nb, nt ) = vqint * pref
        ENDDO
        !
      ENDIF
      !
    ENDDO
    write(*,*) 'nt, sum(tab_at) = ', nt, sum(tab_at(:,:,1))
  ENDDO
  !
  CALL mp_sum( tab_at, intra_bgrp_comm )
  !
  DEALLOCATE( aux, vchi )

  write(*,*)
  write(*,*) '</div> EXIT my_init_at_1'
  write(*,*)

  RETURN
  !
END SUBROUTINE

