!----------------------------------------------------------------------------
SUBROUTINE my_sum_bec( ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd ) 
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


  if(tqr) then
    stop 'tqr is disabled in my_sum_bec'
   endif


  ! counters on beta functions, atoms, atom types, spin
  !
  npw = ngk(ik)
  IF ( .NOT. real_space ) THEN
     ! calbec computes becp = <vkb_i|psi_j>
     CALL calbec( npw, vkb, evc(:,ibnd_start:ibnd_end), becp )
  ELSE
     if (gamma_only) then
        do kbnd = 1, this_bgrp_nbnd, 2 !  ibnd_start, ibnd_end, 2
           ibnd = ibnd_start + kbnd - 1
           call invfft_orbital_gamma(evc,ibnd,ibnd_end) 
           call calbec_rs_gamma(kbnd,this_bgrp_nbnd,becp%r)
        enddo
     else
        current_k = ik
        becp%k = (0.d0,0.d0)
        do kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
           ibnd = ibnd_start + kbnd - 1
           call invfft_orbital_k(evc,ibnd,ibnd_end) 
           call calbec_rs_k(kbnd,this_bgrp_nbnd)
        enddo
     endif
  ENDIF
  !
  ! In the EXX case with ultrasoft or PAW, a copy of becp will be
  ! saved in a global variable to be rotated later
  CALL store_becxx0(ik, becp)
  !
  DO np = 1, ntyp
     !
     IF ( upf(np)%tvanp ) THEN
        !
        ! allocate work space used to perform GEMM operations
        !
        IF ( gamma_only ) THEN
           nbnd_loc = becp%nbnd_loc
           ALLOCATE( auxg( nbnd_loc, nh(np) ) )
        ELSE
           ALLOCATE( auxk1( ibnd_start:ibnd_end, nh(np)*npol ), &
                     auxk2( ibnd_start:ibnd_end, nh(np)*npol ) )
        END IF
        IF ( noncolin ) THEN
           ALLOCATE ( aux_nc( nh(np)*npol,nh(np)*npol ) ) 
        ELSE
           ALLOCATE ( aux_gk( nh(np),nh(np) ) ) 
           if (tqr) ALLOCATE ( aux_egk( nh(np),nh(np) ) ) 
        END IF
        !
        !   In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
        !   run from index i=indv_ijkb0(na)+1 to i=indv_ijkb0(na)+nh(nt)
        !
        DO na = 1, nat
           !
           IF (ityp(na)==np) THEN
              !
              ! sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
              ! copy into aux1, aux2 the needed data to perform a GEMM
              !
              IF ( noncolin ) THEN
                 !
!$omp parallel do default(shared), private(is,ih,ikb,ibnd)
                 DO is = 1, npol
                    DO ih = 1, nh(np)
                       ikb = indv_ijkb0(na) + ih
                       DO kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
                          ibnd = ibnd_start + kbnd - 1
                          auxk1(ibnd,ih+(is-1)*nh(np))= becp%nc(ikb,is,kbnd)
                          auxk2(ibnd,ih+(is-1)*nh(np))= wg(ibnd,ik) * &
                                                        becp%nc(ikb,is,kbnd)
                       END DO
                    END DO
                 END DO
!$omp end parallel do
                 !
                 CALL ZGEMM ( 'C', 'N', npol*nh(np), npol*nh(np), this_bgrp_nbnd, &
                      (1.0_dp,0.0_dp), auxk1, this_bgrp_nbnd, auxk2, this_bgrp_nbnd, &
                      (0.0_dp,0.0_dp), aux_nc, npol*nh(np) )
                 !
              ELSE IF ( gamma_only ) THEN
                 !
!$omp parallel do default(shared), private(ih,ikb,ibnd,ibnd_loc)
                 DO ih = 1, nh(np)
                    ikb = indv_ijkb0(na) + ih
                    DO ibnd_loc = 1, nbnd_loc
                       ibnd = (ibnd_start - 1) + ibnd_loc + becp%ibnd_begin - 1
                       auxg(ibnd_loc,ih)= wg(ibnd,ik)*becp%r(ikb,ibnd_loc) 
                    END DO
                 END DO
!$omp end parallel do
                 CALL DGEMM ( 'N', 'N', nh(np), nh(np), nbnd_loc, &
                      1.0_dp, becp%r(indv_ijkb0(na)+1,1), nkb,    &
                      auxg, nbnd_loc, 0.0_dp, aux_gk, nh(np) )
               if (tqr) then
!$omp parallel do default(shared), private(ih,ikb,ibnd,ibnd_loc)
                 DO ih = 1, nh(np)
                    ikb = indv_ijkb0(na) + ih
                    DO ibnd_loc = 1, nbnd_loc
                       ibnd = (ibnd_start - 1) + ibnd_loc + becp%ibnd_begin - 1
                       auxg(ibnd_loc,ih) = et(ibnd,ik) * auxg(ibnd_loc,ih)
                    END DO
                 END DO
!$omp end parallel do
                 CALL DGEMM ( 'N', 'N', nh(np), nh(np), nbnd_loc, &
                      1.0_dp, becp%r(indv_ijkb0(na)+1,1), nkb,    &
                      auxg, nbnd_loc, 0.0_dp, aux_egk, nh(np) )
               end if
                 !
              ELSE
                 !
!$omp parallel do default(shared), private(ih,ikb,ibnd)
                 DO ih = 1, nh(np)
                    ikb = indv_ijkb0(na) + ih
                    DO kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
                       ibnd = ibnd_start + kbnd - 1
                       auxk1(ibnd,ih) = becp%k(ikb,kbnd) 
                       auxk2(ibnd,ih) = wg(ibnd,ik)*becp%k(ikb,kbnd)
                    END DO
                 END DO
!$omp end parallel do
                 !
                 ! only the real part is computed
                 !
                 CALL DGEMM ( 'C', 'N', nh(np), nh(np), 2*this_bgrp_nbnd, &
                      1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
                      0.0_dp, aux_gk, nh(np) )
                 !
               if (tqr) then
!$omp parallel do default(shared), private(ih,ikb,ibnd)
                 DO ih = 1, nh(np)
                    ikb = indv_ijkb0(na) + ih
                    DO kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
                       ibnd = ibnd_start + kbnd - 1
                       auxk2(ibnd,ih) = et(ibnd,ik)*auxk2(ibnd,ih)
                    END DO
                 END DO
!$omp end parallel do
                 CALL DGEMM ( 'C', 'N', nh(np), nh(np), 2*this_bgrp_nbnd, &
                      1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
                      0.0_dp, aux_egk, nh(np) )
               end if

              END IF
              !
              ! copy output from GEMM into desired format
              !
              IF (noncolin .AND. .NOT. upf(np)%has_so) THEN
                 CALL my_add_becsum_nc (na, np, aux_nc, becsum )
              ELSE IF (noncolin .AND. upf(np)%has_so) THEN
                 CALL my_add_becsum_so (na, np, aux_nc,becsum )
              ELSE
                 ijh = 0
                 DO ih = 1, nh(np)
                    DO jh = ih, nh(np)
                       ijh = ijh + 1
                       !
                       ! nondiagonal terms summed and collapsed into a
                       ! single index (matrix is symmetric wrt (ih,jh))
                       !
                       IF ( jh == ih ) THEN
                          becsum(ijh,na,current_spin) = &
                               becsum(ijh,na,current_spin) + aux_gk (ih,jh)
                          if (tqr) ebecsum(ijh,na,current_spin) = &
                               ebecsum(ijh,na,current_spin) + aux_egk (ih,jh)
                       ELSE
                          becsum(ijh,na,current_spin) = &
                               becsum(ijh,na,current_spin) + aux_gk(ih,jh)*2.0_dp
                          if (tqr) ebecsum(ijh,na,current_spin) = &
                               ebecsum(ijh,na,current_spin) + aux_egk(ih,jh)*2.0_dp
                       END IF
                    END DO
                 END DO
                 !
              END IF
           END IF
           !
        END DO
        !
        IF ( noncolin ) THEN
           DEALLOCATE ( aux_nc )
        ELSE
           DEALLOCATE ( aux_gk  ) 
           if (tqr) DEALLOCATE ( aux_egk  ) 
        END IF
        IF ( gamma_only ) THEN
           DEALLOCATE( auxg )
        ELSE
           DEALLOCATE( auxk2, auxk1 )
        END IF
        !
     END IF
     !
  END DO
  !
END SUBROUTINE my_sum_bec




!----------------------------------------------------------------------------
SUBROUTINE my_add_becsum_nc ( na, np, becsum_nc, becsum )
  !----------------------------------------------------------------------------
  !! This routine multiplies \(\text{becsum_nc}\) by the identity and the
  !! Pauli matrices, saves it in \(\text{becsum}\) for the calculation of 
  !! augmentation charge and magnetization.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, nspin_mag
  USE spin_orb,             ONLY : domag
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na, np
  COMPLEX(DP), INTENT(IN) :: becsum_nc(nh(np),npol,nh(np),npol)
  REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin_mag)
  !
  ! local variables
  !
  INTEGER :: ih, jh, ijh
  REAL(dp) :: fac
  !
  ijh=0
  DO ih = 1, nh(np)
     DO jh = ih, nh(np)
        ijh=ijh+1
        IF ( ih == jh ) THEN
           fac = 1.0_dp
        ELSE
           fac = 2.0_dp
        END IF
        becsum(ijh,na,1)= becsum(ijh,na,1) + fac * &
                DBLE( becsum_nc(ih,1,jh,1) + becsum_nc(ih,2,jh,2) )
        IF (domag) THEN
           becsum(ijh,na,2)= becsum(ijh,na,2) + fac *  &
                DBLE( becsum_nc(ih,1,jh,2) + becsum_nc(ih,2,jh,1) )
           becsum(ijh,na,3)= becsum(ijh,na,3) + fac * DBLE( (0.d0,-1.d0)* &
               (becsum_nc(ih,1,jh,2) - becsum_nc(ih,2,jh,1)) )
           becsum(ijh,na,4)= becsum(ijh,na,4) + fac * &
                DBLE( becsum_nc(ih,1,jh,1) - becsum_nc(ih,2,jh,2) )
        END IF
     END DO
  END DO
  
END SUBROUTINE my_add_becsum_nc
!
!----------------------------------------------------------------------------
SUBROUTINE my_add_becsum_so( na, np, becsum_nc, becsum )
  !----------------------------------------------------------------------------
  !! This routine multiplies \(\text{becsum_nc}\) by the identity and the Pauli
  !! matrices, rotates it as appropriate for the spin-orbit case, saves it in 
  !! \(\text{becsum}\) for the calculation of augmentation charge and magnetization.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE uspp,                 ONLY : ijtoh, nhtol, nhtoj, indv
  USE noncollin_module,     ONLY : npol, nspin_mag
  USE spin_orb,             ONLY : fcoef, domag
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na
  !! atom index
  INTEGER, INTENT(IN) :: np
  !! atomic type
  COMPLEX(DP), INTENT(IN) :: becsum_nc(nh(np),npol,nh(np),npol)
  !! becsum - noncolin
  REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin_mag)
  !! sum over bands
  !
  ! local variables
  !
  INTEGER :: ih, jh, lh, kh, ijh, is1, is2
  COMPLEX(DP) :: fac
  !
  DO ih = 1, nh(np)
     DO jh = 1, nh(np)
        ijh=ijtoh(ih,jh,np)
        DO kh = 1, nh(np)
           IF (same_lj(kh,ih,np)) THEN
              DO lh=1,nh(np)
                 IF (same_lj(lh,jh,np)) THEN
                    DO is1=1,npol
                       DO is2=1,npol
                          fac=becsum_nc(kh,is1,lh,is2)
                          becsum(ijh,na,1)=becsum(ijh,na,1) + fac * &
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) + &
                               fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  )
                          IF (domag) THEN
                            becsum(ijh,na,2)=becsum(ijh,na,2)+fac * &
                                (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) +&
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  )
                            becsum(ijh,na,3)=becsum(ijh,na,3)+fac*(0.d0,-1.d0)*&
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  )
                           becsum(ijh,na,4)=becsum(ijh,na,4) + fac * &
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  )
                        END IF
                     END DO
                  END DO
               END IF
            END DO
         END IF
      END DO
   END DO
END DO
!
CONTAINS
   LOGICAL FUNCTION same_lj(ih,jh,np)
   INTEGER :: ih, jh, np
   !
   same_lj = ((nhtol(ih,np)==nhtol(jh,np)).AND. &
             (ABS(nhtoj(ih,np)-nhtoj(jh,np))<1.d8).AND. &
             (indv(ih,np)==indv(jh,np)) )
   !
   END FUNCTION same_lj

END SUBROUTINE my_add_becsum_so
