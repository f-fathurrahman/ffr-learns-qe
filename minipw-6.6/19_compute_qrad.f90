!
! Compute interpolation table qrad(i,nm,l+1,nt) = Q^{(L)}_{nm,nt}(q_i)
! of angular momentum L, for atom of type nt, on grid q_i, where
! nm = combined index for n,m=1,nh(nt)
!
!----------------------------------------------------------------------
SUBROUTINE my_compute_qrad()
!----------------------------------------------------------------------
  USE kinds,        ONLY : dp
  USE constants,    ONLY : fpi
  USE ions_base,    ONLY : ntyp => nsp
  USE cell_base,    ONLY : omega
  USE atom,         ONLY : rgrid
  USE uspp_param,   ONLY : upf, lmaxq, nbetam, nh, nhm, lmaxkb
  USE us,           ONLY : nqxq, dq, qrad
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER :: ndm, startq, lastq, nt, l, nb, mb, ijv, iq, ir
  ! various indices
  REAL(dp) :: prefr
  ! the prefactor of the Q functions
  REAL(dp) :: q
  REAL(dp), ALLOCATABLE :: aux(:), besr(:)
  ! various work space
  !
  prefr = fpi / omega
  ndm = MAXVAL( upf(:)%kkbeta )
  ALLOCATE(aux(ndm))
  ALLOCATE(besr(ndm))

  CALL divide (intra_bgrp_comm, nqxq, startq, lastq)
  !
  qrad(:,:,:,:) = 0.d0
  !
  DO nt = 1, ntyp

    IF( upf(nt)%tvanp ) then

      DO l = 0, upf(nt)%nqlc -1
        ! note that l is the true (combined) angular momentum
        ! and that the arrays have dimensions 0..l (no more 1..l+1)
        DO iq = startq, lastq
          !
          q = (iq - 1) * dq
          ! here we compute the spherical bessel function for each q_i
          CALL sph_bes( upf(nt)%kkbeta, rgrid(nt)%r, q, l, besr)
          !
          DO nb = 1, upf(nt)%nbeta
            !  the Q are symmetric with respect to indices
            DO mb = nb, upf(nt)%nbeta
              ijv = mb * (mb - 1) / 2 + nb
              IF ( ( l >= abs(upf(nt)%lll(nb) - upf(nt)%lll(mb)) ) .AND. &
                   ( l <=     upf(nt)%lll(nb) + upf(nt)%lll(mb)  ) .AND. &
                   (mod(l+upf(nt)%lll(nb)+upf(nt)%lll(mb),2)==0) ) THEN
                DO ir = 1, upf(nt)%kkbeta
                   aux(ir) = besr(ir) * upf(nt)%qfuncl(ir,ijv,l)
                ENDDO
                !
                !   and then we integrate with all the Q functions
                !
                CALL simpson( upf(nt)%kkbeta, aux, rgrid(nt)%rab, qrad(iq,ijv,l+1, nt) )
              ENDIF
            ENDDO
          ENDDO
        ! igl
        ENDDO
        ! l
      ENDDO
      qrad(:, :, :, nt) = qrad (:, :, :, nt)*prefr
      CALL mp_sum( qrad (:, :, :, nt), intra_bgrp_comm )
    ENDIF 
  ENDDO ! ntyp
  !
  DEALLOCATE(besr)
  DEALLOCATE(aux)
  !
END SUBROUTINE 
