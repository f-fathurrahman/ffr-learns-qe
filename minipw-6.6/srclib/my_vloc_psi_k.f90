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
  !
  IF( use_tg ) THEN
     !     v_siz =  dffts%nnr_tg
     !
     ALLOCATE( tg_v   ( v_siz ) )
     ALLOCATE( tg_psic( v_siz ) )
     !
     CALL tg_gather( dffts, v, tg_v )
  ENDIF
  !
  IF( use_tg ) THEN

     CALL tg_get_nnr( dffts, right_nnr )

     ! compute the number of chuncks
     numblock  = (n+blocksize-1)/blocksize

     DO ibnd = 1, m, fftx_ntgrp(dffts)
        !

        CALL threaded_barrier_memset(tg_psic, 0.D0, fftx_ntgrp(dffts)*right_nnr*2)

        DO idx = 0, MIN(fftx_ntgrp(dffts)-1, m-ibnd)
           DO j = 1, numblock
              tg_psic(dffts%nl (igk_k((j-1)*blocksize+1:MIN(j*blocksize, n),current_k))+right_nnr*idx) = &
                 psi((j-1)*blocksize+1:MIN(j*blocksize, n),idx+ibnd)
           ENDDO
        ENDDO

        !
        CALL  invfft ('tgWave', tg_psic, dffts )
        !write (6,*) 'wfc R ' 
        !write (6,99) (tg_psic(i), i=1,400)
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !

        DO j = 1, dffts%nr1x*dffts%nr2x* right_nr3
           tg_psic (j) = tg_psic (j) * tg_v(j)
        ENDDO

        !write (6,*) 'v psi R ' 
        !write (6,99) (tg_psic(i), i=1,400)
        !
        CALL fwfft ('tgWave',  tg_psic, dffts )
        !
        !   addition to the total product
        !
        CALL tg_get_recip_inc( dffts, right_inc )
        !

        DO idx = 0, MIN(fftx_ntgrp(dffts)-1, m-ibnd)
           DO j = 1, numblock
              hpsi ((j-1)*blocksize+1:MIN(j*blocksize, n), ibnd+idx) = &
                 hpsi ((j-1)*blocksize+1:MIN(j*blocksize, n), ibnd+idx) + &
                 tg_psic( dffts%nl(igk_k((j-1)*blocksize+1:MIN(j*blocksize, n),current_k)) + right_inc*idx )
           ENDDO
        ENDDO
        !
     ENDDO
  ELSE
     DO ibnd = 1, m
        !
        CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
        DO j = 1, n
           psic (dffts%nl (igk_k(j,current_k))) = psi(j, ibnd)
        ENDDO
        !write (6,*) 'wfc G ', ibnd
        !write (6,99) (psic(i), i=1,400)
        !
        CALL invfft ('Wave', psic, dffts)
        !write (6,*) 'wfc R ' 
        !write (6,99) (psic(i), i=1,400)
        !
        DO j = 1, dffts%nnr
           psic (j) = psic (j) * v(j)
        ENDDO
        !write (6,*) 'v psi R ' 
        !write (6,99) (psic(i), i=1,400)
        !
        CALL fwfft ('Wave', psic, dffts)
        !
        !   addition to the total product
        !

        DO j = 1, n
           hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (dffts%nl(igk_k(j,current_k)))
        ENDDO
        !write (6,*) 'v psi G ', ibnd
        !write (6,99) (psic(i), i=1,400)
        !
     ENDDO
  ENDIF
  !
  IF( use_tg ) THEN
     !
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
     !
  ENDIF
  !
99 format ( 20 ('(',2f12.9,')') )

  RETURN
END SUBROUTINE
