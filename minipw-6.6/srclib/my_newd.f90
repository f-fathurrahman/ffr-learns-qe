MODULE my_dfunct
!
CONTAINS
!
!-------------------------------------------------------------------------
SUBROUTINE my_newq( vr, deeq, skip_vltot )
  !----------------------------------------------------------------------
  !! This routine computes the integral of the perturbed potential with
  !! the Q function
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : omega, tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : g, gg, ngm, gstart, mill, &
                                   eigts1, eigts2, eigts3
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vltot
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : psic
  USE noncollin_module,     ONLY : nspin_mag
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), INTENT(IN)  :: vr( dfftp%nnr, nspin )
  !! Input: potential
  REAL(KIND=DP), INTENT(OUT) :: deeq( nhm, nhm, nat, nspin )
  !! Output: contribution to integral
  LOGICAL, INTENT(IN) :: skip_vltot
  !! If .FALSE. vltot is added to vr when necessary
  !
  ! ... local variables
  !
  ! starting/ending indices, local number of G-vectors
  INTEGER :: ig, nt, ih, jh, na, is, ijh, nij, nb, nab
  ! counters on g vectors, atom type, beta functions x 2,
  !   atoms, spin, aux, aux, beta func x2 (again)
  COMPLEX(DP), ALLOCATABLE :: vaux(:,:), aux(:,:), qgm(:,:)
  ! work space
  REAL(DP), ALLOCATABLE :: ylmk0(:,:), qmod(:), deeaux(:,:)
  ! spherical harmonics, modulus of G
  REAL(DP) :: fact
  !
  IF( gamma_only ) THEN
    fact = 2.0_dp
  ELSE
    fact = 1.0_dp
  ENDIF
  !
  deeq(:,:,:,:) = 0.D0
  

  ! ffr: Parallelization wrt G-vectors is disabled

  ALLOCATE( vaux(ngm,nspin_mag), qmod(ngm), ylmk0(ngm,lmaxq*lmaxq) )
  !
  CALL ylmr2( lmaxq*lmaxq, ngm, g, gg, ylmk0 )
  DO ig = 1, ngm
     qmod(ig) = SQRT(gg(ig))*tpiba
  ENDDO

  !
  ! ... fourier transform of the total effective potential
  !
  DO is = 1, nspin_mag
    IF( (nspin_mag == 4 .AND. is /= 1) .OR. skip_vltot ) THEN 
      DO ig = 1, dfftp%nnr
        psic(ig) = vr(ig,is)
      ENDDO
    ELSE
      DO ig = 1, dfftp%nnr
        psic(ig) = vltot(ig) + vr(ig,is)
      ENDDO
    ENDIF
    !write(*,*) 'my_newq: sum(psic) Veff (in Ha) = ', sum(psic)*0.5d0
    !
    CALL fwfft( 'Rho', psic, dfftp )
    DO ig = 1, ngm
      vaux(ig,is) = psic(dfftp%nl(ig))
    ENDDO
    !write(*,*) 'my_newq: sum vaux = ', sum(vaux)
  ENDDO
  ! vaux is V_eff(G)

  DO nt = 1, ntyp

    IF( upf(nt)%tvanp ) THEN
      !
      ! nij = max number of (ih,jh) pairs per atom type nt
      nij = nh(nt)*(nh(nt)+1)/2
      ALLOCATE( qgm(ngm,nij) )
      !
      ! Compute and store Q(G) for this atomic species 
      ! (without structure factor)
      !
      ijh = 0
      DO ih = 1, nh(nt)
        DO jh = ih, nh(nt)
          ijh = ijh + 1
          CALL qvan2( ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0 )
          !write(*,*) 'ijh = ', ijh, ' sum qgm = ', sum(qgm(:,ijh))
        ENDDO
      ENDDO

      !
      ! count max number of atoms of type nt
      !
      nab = 0
      DO na = 1, nat
        IF( ityp(na) == nt ) nab = nab + 1
      ENDDO
      !
      ALLOCATE( aux(ngm,nab ), deeaux(nij,nab) )
      !
      ! Compute and store V(G) times the structure factor e^(-iG*tau)
      DO is = 1, nspin_mag
        nb = 0
        DO na = 1, nat
          IF( ityp(na) == nt ) THEN
            nb = nb + 1
            DO ig = 1, ngm
              aux(ig, nb) = vaux(ig,is) * CONJG( &
                eigts1(mill(1,ig),na) * &
                eigts2(mill(2,ig),na) * &
                eigts3(mill(3,ig),na) )
            ENDDO
            !write(*,*) 'nb = ', nb, ' sum aux = ', sum(aux(:,nb))
          ENDIF
        ENDDO
        !
        ! here we compute the integral Q*V for all atoms of this kind
        CALL DGEMM( 'C', 'N', nij, nab, 2*ngm, fact, qgm, 2*ngm, aux, &
                 2*ngm, 0.0_dp, deeaux, nij )
        IF ( gamma_only .AND. gstart == 2 ) &
             CALL DGER( nij, nab, -1.0_dp, qgm, 2*ngm, aux, 2*ngm, deeaux, nij )
        !
        nb = 0
        DO na = 1, nat
          IF( ityp(na) == nt ) THEN
            nb = nb + 1
            ijh = 0
            DO ih = 1, nh(nt)
              DO jh = ih, nh(nt)
                ijh = ijh + 1
                deeq(ih,jh,na,is) = omega * deeaux(ijh,nb)
                !write(*,'(4I4,F18.10)') ih, jh, na, is, deeq(ih, jh, na, is)*0.5d0 ! to Ha
                IF (jh > ih) deeq(jh,ih,na,is) = deeq(ih,jh,na,is)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        !
      ENDDO ! loop is: over nspin_mag
      !
      DEALLOCATE( deeaux, aux, qgm )
      !
    ENDIF ! tvanp
    !
  ENDDO ! over ntyp
  

  write(*,*) '**** my_newq: sum(deeq) (in Ha) = ', sum(deeq)*0.5d0

  DEALLOCATE( qmod, ylmk0, vaux )
  
  !CALL mp_sum( deeq( :, :, :, 1:nspin_mag ), inter_pool_comm )
  !CALL mp_sum( deeq( :, :, :, 1:nspin_mag ), intra_bgrp_comm )
  
  !
END SUBROUTINE



!----------------------------------------------------------------------------
SUBROUTINE my_newd()
  !----------------------------------------------------------------------------
  !! This routine computes the integral of the effective potential with
  !! the Q function and adds it to the bare ionic D term which is used
  !! to compute the non-local term in the US scheme.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,             ONLY : nspin
  USE uspp,                 ONLY : deeq, dvan, deeq_nc, dvan_so, okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE spin_orb,             ONLY : lspinorb, domag
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE uspp,                 ONLY : nhtol, nhtolm
  USE scf,                  ONLY : v
  USE realus,               ONLY : newq_r
  USE control_flags,        ONLY : tqr
  USE ldaU,                 ONLY : lda_plus_U, U_projection
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, nt, ih, jh, na, is, nht, nb, mb
  ! counters on g vectors, atom type, beta functions x 2,
  !   atoms, spin, aux, aux, beta func x2 (again)
  !
  
  write(*,*)
  write(*,*) '---------------------------------------------------'
  write(*,*) 'ENTER my_newd'
  write(*,*) '---------------------------------------------------'
  
  ! For debugging
  !DO na = 1, nat
  !  nt  = ityp(na)
  !  nht = nh(nt)
  !  DO is = 1, nspin
  !    deeq(1:nht,1:nht,na,is) = dvan(1:nht,1:nht,nt)
  !  ENDDO
  !enddo
  !write(*,*) 'Early return, set to bare coefs'
  !return


  IF( .NOT. okvan ) THEN
    ! no ultrasoft potentials: use bare coefficients for projectors
    DO na = 1, nat
      nt  = ityp(na)
      nht = nh(nt)
      DO is = 1, nspin
         deeq(1:nht,1:nht,na,is) = dvan(1:nht,1:nht,nt)
      ENDDO
    ENDDO
    ! early return
    RETURN
  ENDIF


  ! Calculate Q integral here
  IF(tqr) THEN
    stop 'tqr in my_newd is disabled'
  ELSE
    CALL my_newq( v%of_r, deeq, .FALSE. )
  ENDIF

  IF( noncolin ) CALL add_paw_to_deeq( deeq )

  write(*,*) 'After newq: '
  write(*,*) 'sum Dvan = ', sum(dvan)
  !write(*,*) 'Some Deeq'
  
  atoms : &
  DO na = 1, nat
    nt  = ityp(na)
    if_noncolin:&
    IF( noncolin ) THEN
      stop 'Disabled 291 in my_newd'
    ELSE if_noncolin
      DO is = 1, nspin
         DO ih = 1, nh(nt)
            DO jh = ih, nh(nt)
              deeq(ih,jh,na,is) = deeq(ih,jh,na,is) + dvan(ih,jh,nt)
              deeq(jh,ih,na,is) = deeq(ih,jh,na,is)
              !write(*,'(4I4,F18.10)') ih, jh, na, is, deeq(ih, jh, na, is)
            ENDDO
         ENDDO
      ENDDO
     ENDIF if_noncolin
  ENDDO atoms

  IF(.NOT. noncolin) CALL add_paw_to_deeq(deeq)

  IF (lda_plus_U .AND. (U_projection == 'pseudo')) CALL add_vhub_to_deeq( deeq )

  write(*,*) 'sum Deeq = ', sum(deeq)


  write(*,*)
  write(*,*) '---------------------------------------------------'
  write(*,*) 'EXIT my_newd'
  write(*,*) '---------------------------------------------------'


  !
  RETURN
  !
!  CONTAINS


END SUBROUTINE my_newd
!
!
END MODULE my_dfunct
