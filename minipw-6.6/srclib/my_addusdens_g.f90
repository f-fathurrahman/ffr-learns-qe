!----------------------------------------------------------------------
SUBROUTINE my_addusdens_g( rho )
!----------------------------------------------------------------------
!! This routine adds to the charge density \(\text{rho}(G)\) in reciprocal space
!! the part which is due to the US augmentation.

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, gg, g, &
                                   eigts1, eigts2, eigts3, mill
  USE noncollin_module,     ONLY : nspin_mag
  USE uspp,                 ONLY : becsum, okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: rho(dfftp%ngm,nspin_mag)
  

  !
  ! ... local variables
  !
  
  ! starting/ending indices, local number of G-vectors
  INTEGER :: ngm_s, ngm_e, ngm_l, ngm_s_tmp, ngm_e_tmp, ngm_l_tmp
  
  ! counters
  INTEGER :: ig, na, nt, ih, jh, ijh, is, nab, nb, nij
  
  ! \sum_kv <\psi_kv|\beta_l> <beta_m|\psi_kv> for each species of atoms
  REAL(DP), ALLOCATABLE :: tbecsum(:,:,:)

  ! modulus of G, spherical harmonics  
  REAL(DP), ALLOCATABLE :: qmod (:), ylmk0 (:,:)

  ! structure factors, US contribution to rho
  COMPLEX(DP), ALLOCATABLE :: skk(:,:), aux2(:,:)

  ! work space for rho(G,nspin), Fourier transform of q
  COMPLEX(DP), ALLOCATABLE ::  aux (:,:), qgm(:)

  ! Early return if not using Vanderbilt USPP
  IF (.NOT. okvan) RETURN



  write(*,*)
  write(*,*) 'Entering my_addusdens_g'
  write(*,*)

  !
  ALLOCATE( aux(ngm,nspin_mag) )
  aux(:,:) = (0.d0, 0.d0)

  ! With k-point/bgrp parallelization, distribute G-vectors across all processors
  ! ngm_s = index of first G-vector for this processor (in the k-point x bgrp pool)
  ! ngm_e = index of last  G-vector for this processor (in the k-point x bgrp pool)
  ! ngm_l = local number of G-vectors 
  CALL divide( inter_pool_comm, ngm, ngm_s_tmp, ngm_e_tmp )
  ngm_l_tmp = ngm_e_tmp - ngm_s_tmp + 1
  CALL divide( inter_bgrp_comm, ngm_l_tmp, ngm_s, ngm_e )
  ngm_l = ngm_e - ngm_s + 1 
  ngm_s = ngm_s + ngm_s_tmp - 1
  ngm_e = ngm_e + ngm_s_tmp - 1

  ! for the extraordinary unlikely case of more processors than G-vectors
  IF ( ngm_l <= 0 ) GOTO 10
  

  ! Allocate memory
  ALLOCATE( qmod(ngm_l), qgm(ngm_l)   )
  ALLOCATE( ylmk0(ngm_l, lmaxq*lmaxq) )
  write(*,*) 'ngm_l = ', ngm_l
  write(*,*) 'lmaxq = ', lmaxq


  CALL ylmr2( lmaxq*lmaxq, ngm_l, g(1,ngm_s), gg(ngm_s), ylmk0 )
  DO ig = 1, ngm_l
    qmod(ig) = SQRT(gg(ngm_s+ig-1))*tpiba
  ENDDO
  !
  DO nt = 1, ntyp

    IF( upf(nt)%tvanp ) THEN

      ! nij = max number of (ih,jh) pairs per atom type nt
      nij = nh(nt)*(nh(nt)+1)/2

      ! count max number of atoms of type nt
      nab = 0
      DO na = 1, nat
        IF( ityp(na) == nt ) nab = nab + 1
      ENDDO
      
      write(*,*) 'nab = ', nab
      write(*,*) 'nij = ', nij
      
      ALLOCATE( skk(ngm_l,nab), tbecsum(nij,nab,nspin_mag), aux2(ngm_l,nij) )

      nb = 0
      DO na = 1, nat
        IF( ityp(na) == nt ) THEN
          nb = nb + 1
          tbecsum(:,nb,:) = becsum(1:nij,na,1:nspin_mag)
          DO ig = 1, ngm_l
             skk(ig,nb) = eigts1(mill(1,ngm_s+ig-1),na) * &
                          eigts2(mill(2,ngm_s+ig-1),na) * &
                          eigts3(mill(3,ngm_s+ig-1),na)
          ENDDO
        ENDIF
      ENDDO

      DO is = 1, nspin_mag
        ! sum over atoms
        CALL dgemm( 'N', 'T', 2*ngm_l, nij, nab, 1.0_dp, skk, 2*ngm_l, &
                       tbecsum(1,1,is), nij, 0.0_dp, aux2, 2*ngm_l )

        ! sum over lm indices of Q_{lm}
        ijh = 0
        DO ih = 1, nh(nt)
          DO jh = ih, nh(nt)
            ijh = ijh + 1
            ! qgm is complex here
            CALL my_qvan2( ngm_l, ih, jh, nt, qmod, qgm, ylmk0 )
            !
            write(*,*) 'ih = ', ih, ' jh = ', jh, ' ijh = ', ijh
            !
            DO ig = 1, ngm_l
              aux(ngm_s+ig-1,is) = aux(ngm_s+ig-1,is) + aux2(ig,ijh)*qgm(ig)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE( aux2, tbecsum, skk )
    ENDIF
  ENDDO
  !
  DEALLOCATE( ylmk0 )
  DEALLOCATE( qgm, qmod )
  !
  10 CONTINUE
  CALL mp_sum( aux, inter_bgrp_comm )
  CALL mp_sum( aux, inter_pool_comm )
  !
  ! add aux to the charge density in reciprocal space
  !
  rho(:,:) = rho(:,:) + aux(:,:)
  !
  DEALLOCATE( aux )

  RETURN
  !
END SUBROUTINE
