!----------------------------------------------------------------------
SUBROUTINE my_addusforce( forcenl )
!----------------------------------------------------------------------
  !! Wrapper to \(\texttt{addusforce_g}\) or \(\texttt{addusforce_r}\).
  !
  USE kinds,         ONLY : dp
  USE ions_base,     ONLY : nat
  USE control_flags, ONLY : tqr
  USE realus,        ONLY : addusforce_r
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: forcenl(3,nat)
  !! the non-local contribution to the force
  !
  IF( tqr ) THEN
    stop 'tqr is not supported in my_addusforce'
  ELSE
    CALL my_addusforce_g( forcenl )
  ENDIF
  !
END SUBROUTINE my_addusforce

!
!! This routine computes the contribution to atomic forces due
!! to the dependence of the Q function on the atomic position.
!! \[ F_{j,\text{at}} = \sum_G \sum_{lm} iG_j\ \text{exp}(-iG*R_\text{at})
!!    V^*(G)\ Q_{lm}(G)\ \text{becsum}(lm,\text{at}) \]
!! where:
!! \[ \text{becsum}(lm,\text{at}) = \sum_i \langle \psi_i|\beta_l\rangle
!!    w_i\langle \beta_m|\psi_i\rangle \]
!! On output: the contribution is added to \(\text{forcenl}\).
!
!----------------------------------------------------------------------
SUBROUTINE my_addusforce_g( forcenl )
!----------------------------------------------------------------------
  USE kinds,              ONLY : DP
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE cell_base,          ONLY : omega, tpiba
  USE fft_base,           ONLY : dfftp
  USE gvect,              ONLY : ngm, gg, g, eigts1, eigts2, eigts3, mill
  USE noncollin_module,   ONLY : nspin_mag
  USE scf,                ONLY : v, vltot
  USE uspp,               ONLY : becsum, okvan
  USE uspp_param,         ONLY : upf, lmaxq, nh
  USE mp_bands,           ONLY : intra_bgrp_comm
  USE mp_pools,           ONLY : inter_pool_comm
  USE mp,                 ONLY : mp_sum
  USE control_flags,      ONLY : gamma_only
  USE fft_interfaces,     ONLY : fwfft
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: forcenl(3, nat)
  !! the non-local contribution to the force
  !
  ! ... local variables
  !
  ! starting/ending indices, local number of G-vectors
  INTEGER  :: ig, nt, ih, jh, ijh, nij, ipol, is, na, nb, nab
  REAL(DP) :: fact
  COMPLEX(DP) :: cfac
  ! work space
  COMPLEX(DP), ALLOCATABLE :: aux(:), aux1(:,:,:), vg(:,:), qgm(:,:)
  REAL(DP),    ALLOCATABLE :: ddeeq(:,:,:,:), qmod(:), ylmk0(:,:), forceq(:,:)
  !
  IF(.NOT. okvan) RETURN

  write(*,*)
  write(*,*) 'ENTER my_addusforce_g'
  write(*,*)

  !
  ALLOCATE( forceq(3,nat) ) ! main quantity of interest
  forceq(:,:) = 0.0_dp
  IF ( gamma_only ) THEN
     fact = 2.d0*omega
  ELSE
     fact = 1.d0*omega
  ENDIF
  !
  ! fourier transform of the total effective potential
  !
  ALLOCATE( vg(ngm,nspin_mag) )
  ALLOCATE( aux(dfftp%nnr) )
  DO is = 1, nspin_mag
    IF(nspin_mag == 4 .AND. is /= 1) THEN
      aux(:) = v%of_r(:,is)
    ELSE
      aux(:) = vltot(:) + v%of_r(:,is)
    ENDIF
    write(*,*) 'sum aux before fwfft in (Ha) = ', sum(aux)*0.5d0
    CALL fwfft( 'Rho', aux, dfftp )
    write(*,*) 'sum aux after fwfft in (Ha) = ', sum(aux)*0.5d0
    ! Note the factors -i and 2pi/a *units of G) here in V(G) !
    vg(:,is) = aux(dfftp%nl(:)) * tpiba * (0.d0, -1.d0)
    write(*,*) 'sum vg (in Ha) * 2*pi = ', sum(vg)*0.5d0
  ENDDO
  DEALLOCATE( aux )
  ! Finish calculation -im*V_eff(G)


  ! Parallelization over G is disabled

  ALLOCATE( ylmk0(ngm, lmaxq*lmaxq) )
  CALL ylmr2( lmaxq * lmaxq, ngm, g(:,1:ngm), gg(1:ngm), ylmk0 )
  !
  ALLOCATE( qmod(ngm) )
  DO ig = 1, ngm
    qmod(ig) = SQRT( gg(ig) )*tpiba
  ENDDO
  !
  DO nt = 1, ntyp

    write(*,*)
    write(*,*) 'Begin species loop = ', nt

    IF ( upf(nt)%tvanp ) THEN
      !
      ! nij = max number of (ih,jh) pairs per atom type nt
      ! qgm contains the Q functions in G space
      !
      nij = nh(nt)*(nh(nt)+1)/2
      ALLOCATE( qgm(ngm,nij) )
      ijh = 0
      DO ih = 1, nh(nt)
        DO jh = ih, nh(nt)
          ijh = ijh + 1
          CALL qvan2( ngm, ih, jh, nt, qmod, qgm(:,ijh), ylmk0 )
          !write(*,*) 'ijh = ', ijh, ' sum(qgm) = ', sum(qgm(:,ijh))
        ENDDO
      ENDDO
      !
      !
      ! nab = number of atoms of type nt
      !
      nab = 0
      DO na = 1, nat
        IF( ityp(na) == nt ) nab = nab + 1
      ENDDO
      ALLOCATE( aux1(ngm, nab, 3) )
      ALLOCATE( ddeeq(nij, nab, 3, nspin_mag) )
      !
      DO is = 1, nspin_mag
        nb = 0
        DO na = 1, nat
          IF (ityp(na) == nt) THEN
            nb = nb + 1
            !
            ! aux1 = product of potential, structure factor and iG
            !
            DO ig = 1, ngm
              cfac = vg(ig, is) * &
                   CONJG(eigts1(mill(1,ig),na) * eigts2(mill(2,ig),na) * eigts3(mill(3,ig),na) )
              aux1(ig,nb,1) = g(1,ig) * cfac
              aux1(ig,nb,2) = g(2,ig) * cfac
              aux1(ig,nb,3) = g(3,ig) * cfac
            ENDDO
            !
          ENDIF
        ENDDO
        write(*,*) 'sum aux1 in (Ha) = ', sum(aux1)*0.5d0
        !
        ! ddeeq = dot product of aux1 with the Q functions
        ! No need for special treatment of the G=0 term (is zero)
        !
        write(*,*) 'sum qgm before dgemm: ', sum(qgm)
        DO ipol = 1, 3
          write(*,*)
          write(*,*) 'ipol = ', ipol, ' sum(aux1) ', sum(aux1(:,:,ipol))
          CALL DGEMM( 'C', 'N', nij, nab, 2*ngm, fact, qgm, 2*ngm, &
               aux1(1,1,ipol), 2*ngm, 0.0_dp, ddeeq(1,1,ipol,is), nij )
          write(*,*) 'sum Ddeeq after DGEMM = ', sum(ddeeq(:,:,ipol,is))
        ENDDO
        !
      ENDDO
      !
      DEALLOCATE(aux1)
      DEALLOCATE(qgm)
      !
      DO is = 1, nspin_mag
        nb = 0
        DO na = 1, nat
          IF( ityp(na) == nt ) THEN
            nb = nb + 1
            DO ipol = 1, 3
              DO ijh = 1, nij
                 forceq(ipol,na) = forceq(ipol,na) + ddeeq(ijh,nb,ipol,is) * becsum(ijh,na,is)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE( ddeeq )
    ENDIF
  ENDDO
  !
  !10 CONTINUE ! NOT USED???
  !
  ! Not required for serial calculation
  CALL mp_sum( forceq, inter_pool_comm )
  CALL mp_sum( forceq, intra_bgrp_comm )
  !
  forcenl(:,:) = forcenl(:,:) + forceq(:,:)
  !
  DEALLOCATE(qmod)
  DEALLOCATE(ylmk0)
  DEALLOCATE(vg)
  DEALLOCATE(forceq)
  
  write(*,*)
  write(*,*) 'EXIT my_addusforce_g'
  write(*,*)

  !
  RETURN
  !
END SUBROUTINE my_addusforce_g
