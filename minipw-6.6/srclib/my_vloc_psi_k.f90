!-----------------------------------------------------------------------
SUBROUTINE my_vloc_psi_k( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - k-points:
  !
  !! * fft to real space;
  !! * product with the potential v on the smooth grid;
  !! * back to reciprocal space;
  !! * addition to the hpsi.
  !
  USE parallel_include
  USE kinds,                  ONLY : DP
  USE wvfct,                  ONLY : current_k
  USE klist,                  ONLY : igk_k
  USE mp_bands,               ONLY : me_bgrp
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : fwfft, invfft
  USE fft_helper_subroutines, ONLY : fftx_ntgrp, tg_get_nnr, &
                                     tg_get_group_nr3, tg_get_recip_inc
  USE wavefunctions,          ONLY : psic
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,m)
  !! Hamiltonian dot psi
  REAL(DP), INTENT(IN) :: v(dffts%nnr)
  !! the total pot. in real space (smooth grid) for current spin
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, incr
  INTEGER :: i, right_nnr, right_nr3, right_inc
  !
  ! chunking parameters
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
  !
  ! Task Groups
  LOGICAL :: use_tg
  REAL(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx

  use_tg = dffts%has_task_groups 
  IF( use_tg ) THEN
    stop 'use_tg is disabled in my_vloc_psi_k'
  ENDIF

  DO ibnd = 1, m
    !
    CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
    DO j = 1, n
      psic( dffts%nl(igk_k(j,current_k)) ) = psi(j, ibnd)
    ENDDO
    !
    CALL invfft('Wave', psic, dffts)
    !
    DO j = 1, dffts%nnr
      psic(j) = psic(j) * v(j)
    ENDDO
    !
    CALL fwfft('Wave', psic, dffts)
    !
    !   addition to the total product
    !
    DO j = 1, n
       hpsi(j,ibnd) = hpsi(j,ibnd) + psic(dffts%nl(igk_k(j,current_k)))
    ENDDO

  ENDDO

  !write(*,*) 'my_vloc_psi_k is called: sum(hpsi) = sum(hpsi)'

  RETURN
END SUBROUTINE
