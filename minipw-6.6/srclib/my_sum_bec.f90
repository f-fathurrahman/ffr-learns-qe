!----------------------------------------------------------------------------
SUBROUTINE my_sum_bec ( ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd ) 
!----------------------------------------------------------------------------
  !! This routine computes the sum over bands:
  !
  !! \[ \sum_i \langle\psi_i|\beta_l\rangle w_i \langle\beta_m|\psi_i\rangle \]
  !
  !! for point "ik" and, for LSDA, spin "current_spin".  
  !! Calls calbec to compute \(\text{"becp"}=\langle \beta_m|\psi_i \rangle\).  
  !! Output is accumulated (unsymmetrized) into "becsum", module "uspp".
  !
  !! Routine used in sum_band (if okvan) and in compute_becsum, called by hinit1 (if okpaw).
  !
  USE kinds,         ONLY : DP
  USE becmod,        ONLY : becp, calbec, allocate_bec_type
  USE control_flags, ONLY : gamma_only, tqr
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE uspp,          ONLY : nkb, vkb, becsum, ebecsum, indv_ijkb0
  USE uspp_param,    ONLY : upf, nh, nhm
  USE wvfct,         ONLY : nbnd, wg, et, current_k
  USE klist,         ONLY : ngk
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions, ONLY : evc
  USE realus,        ONLY : real_space, &
                            invfft_orbital_gamma, calbec_rs_gamma, &
                            invfft_orbital_k, calbec_rs_k
  USE us_exx,        ONLY : store_becxx0
  USE mp_bands,      ONLY : nbgrp,inter_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! the point the sum is computed for
  INTEGER, INTENT(IN) :: current_spin
  !! the input spin
  INTEGER, INTENT(IN) :: ibnd_start
  !! first band of the parallel block
  INTEGER, INTENT(IN) :: ibnd_end
  !! last band of the parallel block
  INTEGER, INTENT(IN) :: this_bgrp_nbnd
  !! number of bands of the parallel block
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE :: auxk1(:,:), auxk2(:,:), aux_nc(:,:)
  REAL(DP), ALLOCATABLE :: auxg(:,:), aux_gk(:,:), aux_egk(:,:)
  INTEGER :: ibnd, ibnd_loc, nbnd_loc, kbnd  ! counters on bands
  INTEGER :: npw, ikb, jkb, ih, jh, ijh, na, np, is, js


  write(*,*)
  write(*,*) '>>>> Calling my_sum_bec'
  write(*,*)

  ! ffr: tqr is disabled

  IF( noncolin ) THEN
    STOP 'noncolin is disabled in my_sum_bec'
  !
  ELSEIF( gamma_only ) THEN
    STOP 'gamma_only is disabled in my_sum_bec'
  ENDIF


  ! counters on beta functions, atoms, atom types, spin
  !
  npw = ngk(ik)
  IF( .NOT. real_space ) THEN
    ! calbec computes becp = <vkb_i|psi_j>
    CALL calbec( npw, vkb, evc(:,ibnd_start:ibnd_end), becp )
  ELSE
    stop 'Disabled in my_sum_bec'
  ENDIF


  ! In the EXX case with ultrasoft or PAW, a copy of becp will be
  ! saved in a global variable to be rotated later
  CALL store_becxx0(ik, becp)

  !
  DO np = 1, ntyp
    !
    IF( upf(np)%tvanp ) THEN
      !
      ! allocate work space used to perform GEMM operations
      IF( gamma_only ) THEN
        nbnd_loc = becp%nbnd_loc
        ALLOCATE( auxg( nbnd_loc, nh(np) ) )
      ELSE
        ALLOCATE( auxk1( ibnd_start:ibnd_end, nh(np)*npol ), &
                  auxk2( ibnd_start:ibnd_end, nh(np)*npol ) )
      ENDIF
      !
      IF( noncolin ) THEN
        ALLOCATE( aux_nc( nh(np)*npol,nh(np)*npol ) ) 
      ELSE
        ALLOCATE( aux_gk( nh(np),nh(np) ) ) 
        if(tqr) ALLOCATE ( aux_egk( nh(np),nh(np) ) ) 
      ENDIF
      !
      ! In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
      ! run from index i=indv_ijkb0(na)+1 to i=indv_ijkb0(na)+nh(nt)
      DO na = 1, nat
        !
        IF( ityp(na)==np ) THEN
          !
          ! sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
          ! copy into aux1, aux2 the needed data to perform a GEMM
          !
          DO ih = 1, nh(np)
            ikb = indv_ijkb0(na) + ih
            DO kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
              ibnd = ibnd_start + kbnd - 1
              auxk1(ibnd,ih) = becp%k(ikb,kbnd) 
              auxk2(ibnd,ih) = wg(ibnd,ik)*becp%k(ikb,kbnd)
            ENDDO
          ENDDO
          !
          ! only the real part is computed
          !
          CALL DGEMM ( 'C', 'N', nh(np), nh(np), 2*this_bgrp_nbnd, &
               1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
               0.0_dp, aux_gk, nh(np) )
          !
          ! copy output from GEMM into desired format
          !
          ijh = 0
          DO ih = 1, nh(np)
            DO jh = ih, nh(np)
              ijh = ijh + 1
              !
              ! nondiagonal terms summed and collapsed into a
              ! single index (matrix is symmetric wrt (ih,jh))
              !
              IF ( jh == ih ) THEN
                becsum(ijh,na,current_spin) = becsum(ijh,na,current_spin) + aux_gk (ih,jh)
              ELSE
                becsum(ijh,na,current_spin) = becsum(ijh,na,current_spin) + aux_gk(ih,jh)*2.0_dp
              ENDIF
            ENDDO
          ENDDO
        ENDIF  ! ityp(na) == np
      ENDDO
      !
      DEALLOCATE( aux_gk  ) 
      DEALLOCATE( auxk2, auxk1 )
    ENDIF
  ENDDO
  !
END SUBROUTINE my_sum_bec
