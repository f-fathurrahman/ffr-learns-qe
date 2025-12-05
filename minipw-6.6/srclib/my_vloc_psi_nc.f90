!-----------------------------------------------------------------------
SUBROUTINE my_vloc_psi_nc( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - noncollinear.
  !
  USE parallel_include
  USE kinds,                  ONLY : DP
  USE wvfct,                  ONLY : current_k
  USE klist,                  ONLY : igk_k
  USE fft_base,               ONLY : dffts, dfftp
  USE fft_interfaces,         ONLY : fwfft, invfft
  USE spin_orb,               ONLY : domag
  USE noncollin_module,       ONLY : npol
  USE wavefunctions,          ONLY : psic_nc
  USE fft_helper_subroutines, ONLY : fftx_ntgrp, tg_get_nnr, &
                                     tg_get_group_nr3, tg_get_recip_inc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  REAL(DP), INTENT(IN) :: v(dfftp%nnr,4) ! beware dimensions!
  !! the total pot. in real space (smooth grid)
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,npol,m)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j,ipol, incr
  COMPLEX(DP) :: sup, sdwn
  !
  ! Variables for task groups
  LOGICAL :: use_tg

  write(*,*) 'my_vloc_psi_nc is called'

  incr = 1
  use_tg = dffts%has_task_groups 
  IF( use_tg ) THEN
    stop 'use tg is disable in my_vloc_psi_nc'
  ENDIF
  !
  ! the local potential V_Loc psi. First the psi in real space
  !
  DO ibnd = 1, m, incr
    psic_nc = (0.d0,0.d0)
    DO ipol=1,npol
      DO j = 1, n
        psic_nc(dffts%nl(igk_k(j,current_k)),ipol) = psi(j+(ipol-1)*lda,ibnd)
      ENDDO
      CALL invfft('Wave', psic_nc(:,ipol), dffts)
    ENDDO
    !
    !   product with the potential v = (vltot+vr) on the smooth grid
    !
    IF (domag) THEN
      DO j=1, dffts%nnr
         sup = psic_nc(j,1) * (v(j,1)+v(j,4)) + psic_nc(j,2) * (v(j,2)-(0.d0,1.d0)*v(j,3))
         sdwn = psic_nc(j,2) * (v(j,1)-v(j,4)) + psic_nc(j,1) * (v(j,2)+(0.d0,1.d0)*v(j,3))
         psic_nc(j,1) = sup
         psic_nc(j,2) = sdwn
      ENDDO
    ELSE
      DO j=1, dffts%nnr
        psic_nc(j,:) = psic_nc(j,:) * v(j,1)
      ENDDO
    ENDIF
    !
    ! back to reciprocal space
    !
    DO ipol=1,npol
       CALL fwfft ('Wave', psic_nc(:,ipol), dffts)
    ENDDO
    !
    ! addition to the total product
    !
    DO ipol=1,npol
      DO j = 1, n
        hpsi(j,ipol,ibnd) = hpsi(j,ipol,ibnd) + psic_nc(dffts%nl(igk_k(j,current_k)),ipol)
      ENDDO
    ENDDO
  ENDDO
  write(*,*) 'current_k = ', current_k
  write(*,*) 'sum hpsi in my_vloc_psi_nc = ', sum(hpsi)
  !
  RETURN
END SUBROUTINE my_vloc_psi_nc
