!--------------------------------------------------------------------
SUBROUTINE my_s_psi( lda, n, m, psi, spsi )
!--------------------------------------------------------------------
  !! This routine applies the S matrix to m wavefunctions psi and puts 
  !! the results in spsi.
  !! Requires the products of psi with all beta functions in array 
  !! becp(nkb,m) (calculated in h_psi or by calbec).
  !
  !! \(\textit{Wrapper routine}\): performs bgrp parallelization on 
  !! non-distributed bands if suitable and required, calls old S\psi
  !! routine s_psi_ . See comments in h_psi.f90 about band 
  !! parallelization.
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE funct,            ONLY : exx_is_active
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,               ONLY : mp_allgather, mp_size, &
                               mp_type_create_column_section, mp_type_free
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the m wavefunctions
  COMPLEX(DP), INTENT(OUT)::spsi(lda*npol,m)
  !! S matrix dot wavefunctions psi
  !
  ! ... local variables
  !
  INTEGER :: m_start, m_end
  INTEGER :: column_type
  INTEGER, ALLOCATABLE :: displs(:)
  !
  IF (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) THEN
    stop 'Not supported my_s_psi: 41'
  ELSE
     ! don't use band parallelization here
     CALL my_s_psi_( lda, n, m, psi, spsi )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE

!
!----------------------------------------------------------------------------
SUBROUTINE my_s_psi_( lda, n, m, psi, spsi )
!----------------------------------------------------------------------------
  !! This routine applies the S matrix to m wavefunctions psi and puts 
  !! the results in spsi.
  !! Requires the products of psi with all beta functions in array 
  !! becp(nkb,m) (calculated in h_psi or by calbec).
  !
  USE kinds,            ONLY: DP
  USE becmod,           ONLY: becp
  USE uspp,             ONLY: vkb, nkb, okvan, qq_at, qq_so, indv_ijkb0
  USE spin_orb,         ONLY: lspinorb
  USE uspp_param,       ONLY: upf, nh, nhm
  USE ions_base,        ONLY: nat, nsp, ityp
  USE control_flags,    ONLY: gamma_only 
  USE noncollin_module, ONLY: npol, noncolin
  USE realus,           ONLY: real_space, invfft_orbital_gamma,     &
                              fwfft_orbital_gamma, calbec_rs_gamma, &
                              s_psir_gamma, invfft_orbital_k,       &
                              fwfft_orbital_k, calbec_rs_k, s_psir_k
  USE wavefunctions,    ONLY: psic
  USE fft_base,         ONLY: dffts
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the m wavefunctions
  COMPLEX(DP), INTENT(OUT)::spsi(lda*npol,m)
  !! S matrix dot wavefunctions psi
  !
  ! ... local variables
  !
  INTEGER :: ibnd
  !
  !
  ! ... initialize  spsi
  !
  CALL threaded_memcpy( spsi, psi, lda*npol*m*2 )
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
  !
  ! ... The product with the beta functions
  !
  IF( gamma_only ) THEN
    stop 'Not yet supported: my_s_psi 102'
  ELSEIF( noncolin ) THEN
    stop 'Not yet supported: my_s_psi 104'
  ELSE 
    IF( real_space ) THEN
      stop 'Not yet supported: my_s_psi: 107'
    ELSE
      !
      CALL s_psi_k()
      !
    ENDIF
  ENDIF    


  RETURN

CONTAINS
     
!-----------------------------------------------------------------------
SUBROUTINE s_psi_k()
!-----------------------------------------------------------------------
  !! k-points version of \(\textrm{s_psi}\) routine.
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ierr
  ! counters
  COMPLEX(DP), ALLOCATABLE :: ps(:,:), qqc(:,:)
  ! ps = product vkb and psi ; qqc = complex version of qq
  !
  ALLOCATE( ps(nkb,m), STAT=ierr )
  !
  IF( ierr /= 0 ) CALL errore( ' s_psi_k ', ' cannot allocate memory (ps) ', ABS(ierr) )

  write(*,*) 'in s_psi_k: sum becp_k = ', sum(becp%k)*0.5d0 ! to Ha

  !
  DO nt = 1, nsp
    !
    IF ( upf(nt)%tvanp ) THEN
      ! qq is real:  copy it into a complex variable to perform
      ! a zgemm - simple but sub-optimal solution
      ALLOCATE( qqc(nh(nt),nh(nt)) )
      DO na = 1, nat
        IF( ityp(na) == nt ) THEN
          qqc(:,:) = CMPLX( qq_at(1:nh(nt),1:nh(nt),na), 0.0_DP, KIND=DP )
          CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_DP,0.0_DP), &
               qqc, nh(nt), becp%k(indv_ijkb0(na)+1,1), nkb, &
               (0.0_DP,0.0_DP), ps(indv_ijkb0(na)+1,1), nkb )
        ENDIF
      ENDDO
      DEALLOCATE( qqc )
    ELSE
      IF( nh(nt) > 0 ) THEN
        DO na = 1, nat
          IF( ityp(na) == nt ) THEN
            ps(indv_ijkb0(na)+1:indv_ijkb0(na)+nh(nt),1:m) = (0.0_DP,0.0_DP)
          ENDIF
        ENDDO
      ENDIF
      !
    ENDIF
    !
  ENDDO

  write(*,*) 'in s_psi_k: sum(ps) = ', sum(ps)*0.5d0 ! to Ha

  IF ( m == 1 ) THEN
     !
     CALL ZGEMV( 'N', n, nkb, ( 1.D0, 0.D0 ), vkb, &
                 lda, ps, 1, ( 1.D0, 0.D0 ), spsi, 1 )
     !
  ELSE
     !
     CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ), vkb, &
                 lda, ps, nkb, ( 1.D0, 0.D0 ), spsi, lda )
     !
  ENDIF
  !
  DEALLOCATE( ps )
  RETURN

END SUBROUTINE


END SUBROUTINE
