#define ZERO ( 0._dp, 0._dp )

!----------------------------------------------------------------------------
SUBROUTINE my_mix_rho( input_rhout, rhoin, alphamix, dr2, tr2_min, iter, n_iter,&
                    iunmix, conv )
  !----------------------------------------------------------------------------
  !! * Modified Broyden's method for charge density mixing: D.D. Johnson,
  !!   PRB 38, 12807 (1988) ;
  !! * Thomas-Fermi preconditioning described in: Raczkowski, Canning, Wang,
  !!   PRB 64,121101 (2001) ;
  !! * Extended to mix also quantities needed for PAW, meta-GGA, DFT+U(+V) ;
  !! * Electric field (all these are included into \(\text{mix_type}\)) ;
  !! * On output: the mixed density is in \(\text{rhoin}\), 
  !!   \(\text{input_rhout}\) is unchanged.
  !
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat, ntyp => nsp
  USE gvect,          ONLY : ngm
  USE gvecs,          ONLY : ngms
  USE lsda_mod,       ONLY : nspin
  USE control_flags,  ONLY : imix, ngm0, tr2, io_level
  ! ... for PAW:
  USE uspp_param,     ONLY : nhm
  USE scf,            ONLY : scf_type, create_scf_type, destroy_scf_type, &
                             mix_type, create_mix_type, destroy_mix_type, &
                             assign_scf_to_mix_type, assign_mix_to_scf_type, &
                             mix_type_AXPY, davcio_mix_type, rho_ddot, &
                             high_frequency_mixing, nsg_ddot, &
                             mix_type_COPY, mix_type_SCAL
  USE io_global,     ONLY : stdout

  USE io_files,      ONLY : diropn


  !
  IMPLICIT NONE
  !
  ! ... First the I/O variable
  !
  INTEGER, INTENT(IN) :: iter !! counter of the number of iterations
  INTEGER, INTENT(IN) :: n_iter !! number of iterations used in mixing
  INTEGER, INTENT(IN) :: iunmix !! I/O unit where data from previous iterations is stored
  REAL(DP), INTENT(IN) :: alphamix !! mixing factor
  REAL(DP), INTENT(IN) :: tr2_min
  !! estimated error in diagonalization. If the estimated
  !! scf error is smaller than this, exit: a more accurate 
  !! diagonalization is needed
  REAL(DP), INTENT(OUT) :: dr2
  !! the estimated error on the energy
  LOGICAL, INTENT(OUT) :: conv
  !! .TRUE. if the convergence has been reached
  !
  TYPE(scf_type), INTENT(INOUT) :: input_rhout
  TYPE(scf_type), INTENT(INOUT) :: rhoin
  !
  ! ... local variables
  !
  TYPE(mix_type) :: rhout_m, rhoin_m
  INTEGER, PARAMETER :: &
    maxmix = 25     ! max number of iterations for charge mixing
  INTEGER ::       &
    iter_used,     &! actual number of iterations used
    ipos,          &! index of the present iteration
    inext,         &! index of the next iteration
    i, j,          &! counters on number of iterations
    info,          &! flag saying if the exec. of libr. routines was ok
    ldim,          &! 2 * Hubbard_lmax + 1
    iunmix_nsg,    &! the unit for Hubbard mixing within DFT+U+V
    nt              ! index of the atomic type
  REAL(DP),ALLOCATABLE :: betamix(:,:), work(:)
  INTEGER, ALLOCATABLE :: iwork(:)
  COMPLEX(DP), ALLOCATABLE :: nsginsave(:,:,:,:,:),  nsgoutsave(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: deltansg(:,:,:,:,:)
  LOGICAL :: exst
  REAL(DP) :: gamma0

  !
  ! ... saved variables and arrays
  !
  INTEGER, SAVE :: mixrho_iter = 0    ! history of mixing

  ! XXX Uh-oh saved variables ....
  TYPE(mix_type), ALLOCATABLE, SAVE :: &
    df(:),        &! information from preceding iterations
    dv(:)          !     "  "       "     "        "  "

  REAL(DP) :: norm
  INTEGER, PARAMETER :: read_ = -1, write_ = +1
  !
  ! ... external functions
  !
  INTEGER, EXTERNAL :: find_free_unit

  ! LDA+U stuffs are removed


  write(*,*)
  write(*,*) '************* Enter my_mix_rho *************'
  write(*,*)

  !
  ngm0 = ngms
  !
  mixrho_iter = iter
  !
  IF ( n_iter > maxmix ) CALL errore( 'mix_rho', 'n_iter too big', 1 )

  !
  ! define mix_type variables and copy scf_type variables there
  !
  call create_mix_type(rhout_m) 
  call create_mix_type(rhoin_m)
  !
  ! copy several variables from scf_type to mix_type
  call assign_scf_to_mix_type(rhoin, rhoin_m)
  call assign_scf_to_mix_type(input_rhout, rhout_m)
  ! copy?

  ! This is used to test rho_ddot
  ! TODO: move this to separate/dedicated subroutine

!  write(*,*)
!  write(*,*) 'sum rhoin%of_g up = ', sum(rhoin%of_g(:,1))
!  if(nspin == 2) then
!    write(*,*) 'sum rhoin%of_g dn = ', sum(rhoin%of_g(:,2))
!  endif
!  dr2 = rho_ddot( rhoin_m, rhoin_m, ngms )  ! must assign first ??
!  write(*,*) 'test rho_ddot(rhoin_m, rhoin_m) (in Ha) = ', dr2*0.5d0
!
!  write(*,*)
!  write(*,*) 'sum input_rhout%of_g up = ', sum(input_rhout%of_g(:,1))
!  if(nspin == 2) then
!    write(*,*) 'sum input_rhout%of_g dn = ', sum(input_rhout%of_g(:,2))
!  endif
!  dr2 = rho_ddot( rhout_m, rhout_m, ngms )
!  write(*,*) 'test rho_ddot(rhout_m, rhout_m) = ', dr2*0.5d0


  ! compute differences
  call mix_type_AXPY( -1.d0, rhoin_m, rhout_m )  ! rhoout_m <- (-1)*rhoin_m + rhout_m


  ! Y <= A*X + Y
  !
  ! compute the "norm" ?
  dr2 = rho_ddot( rhout_m, rhout_m, ngms )  !!!! this used to be ngm NOT ngms
  !
  IF (dr2 < 0.0_DP) CALL errore('mix_rho', 'negative dr2', 1)
  !
  conv = ( dr2 < tr2 )
  write(*,'(1x,A,2ES18.10)') 'my_mix_rho: dr2, tr2 ', dr2, tr2
  write(*,*) 'my_mix_rho: conv = ', conv

  !
  ! This is the case when convergence is achieved
  !
  IF( conv .OR. dr2 < tr2_min ) THEN
    !
    ! if convergence is achieved or if the self-consistency error (dr2) is
    ! smaller than the estimated error due to diagonalization (tr2_min),
    ! exit and leave rhoin and rhocout unchanged
    !
    IF( ALLOCATED( df ) ) THEN
      DO i=1, n_iter
        call destroy_mix_type(df(i))
      END DO
      DEALLOCATE( df )
    ENDIF
    IF( ALLOCATED( dv ) ) THEN
      DO i=1, n_iter
        call destroy_mix_type(dv(i))
      END DO
      DEALLOCATE( dv )
    ENDIF
    !
    call destroy_mix_type(rhoin_m)
    call destroy_mix_type(rhout_m)
    !
    RETURN
  ENDIF

  !
  ! Not yet converged ....
  !
  ! We allocate vectors of mix_type's: df and dv
  !
  IF( .NOT. ALLOCATED( df ) ) THEN
    ALLOCATE( df( n_iter ) )
    DO i = 1,n_iter
      CALL create_mix_type( df(i) )
    ENDDO
  ENDIF
  
  IF( .NOT. ALLOCATED( dv ) ) THEN
    ALLOCATE( dv( n_iter ) )
    DO i=1,n_iter
      CALL create_mix_type( dv(i) )
    ENDDO
  ENDIF
  !
  ! iter_used = mixrho_iter-1  if  mixrho_iter <= n_iter
  ! iter_used = n_iter         if  mixrho_iter >  n_iter
  !
  iter_used = MIN( ( mixrho_iter - 1 ), n_iter )
  !
  ! ipos is the position in which results from the present iteration
  ! are stored. ipos=mixrho_iter-1 until ipos=n_iter, then back to 1,2,...
  !
  ipos = mixrho_iter - 1 - ( ( mixrho_iter - 2 ) / n_iter ) * n_iter
  !
  write(*,*) 'my_mix_rho: iter        = ', iter
  write(*,*) 'my_mix_rho: n_iter      = ', n_iter
  write(*,*) 'my_mix_rho: mixrho_iter = ', mixrho_iter
  write(*,*) 'my_mix_rho: iter_used   = ', iter_used
  write(*,*) 'my_mix_rho: ipos        = ', ipos
  !
  !
  IF ( mixrho_iter > 1 ) THEN
    ! IO stuffs (using disk?)
    CALL davcio_mix_type( df(ipos), iunmix, 1, read_ ) ! read ?? rhout_m previous
    CALL davcio_mix_type( dv(ipos), iunmix, 2, read_ ) ! read ?? rhoin_m previous
    !write(200+iter,*) df(ipos)%of_g
    !write(300+iter,*) dv(ipos)%of_g
    !
    call mix_type_AXPY( -1.d0, rhout_m, df(ipos) )
    call mix_type_AXPY( -1.d0, rhoin_m, dv(ipos) )
  ENDIF
  !
  ! load data other dat for df and dv at positions other than ipos
  DO i = 1, iter_used
    IF ( i /= ipos ) THEN
      !write(*,'(1x,A,3I4)') 'iter, iter_used, will_read i = ', iter, iter_used, i
      CALL davcio_mix_type( df(i), iunmix, 2*i+1, read_ )
      CALL davcio_mix_type( dv(i), iunmix, 2*i+2, read_ )
    ENDIF
  ENDDO
  !
  ! write data for next iter?
  CALL davcio_mix_type( rhout_m, iunmix, 1, write_ )  ! will be read as df(ipos) ?
  !write(400+iter,*) rhout_m%of_g
  CALL davcio_mix_type( rhoin_m, iunmix, 2, write_ )
  !write(500+iter,*) rhoin_m%of_g
  !
  ! write for next iter?
  IF ( mixrho_iter > 1 ) THEN
    CALL davcio_mix_type( df(ipos), iunmix, 2*ipos+1, write_ )
    CALL davcio_mix_type( dv(ipos), iunmix, 2*ipos+2, write_ )
  ENDIF

  !
  ! Nothing else to do on first iteration
  ! FIXME: Is this label used?
  skip_on_first: &
  IF (iter_used > 0) THEN
    !
    ALLOCATE(betamix(iter_used, iter_used)) !iter_used))
    betamix = 0._dp
    !
    DO i = 1, iter_used
      DO j = i, iter_used
        betamix(i,j) = rho_ddot( df(j), df(i), ngm0 )
        betamix(j,i) = betamix(i,j)
      ENDDO
    ENDDO
    !write(100+iter,*) betamix(1:iter_used,1:iter_used)
    !write(*,*) 'betamix written to files'
    !
    ! allocate(e(iter_used), v(iter_used, iter_used))
    ! CALL rdiagh(iter_used, betamix, iter_used, e, v)
    ! write(*,'(1e11.3)') e(:)
    ! write(*,*)
    ! deallocate(e,v)
    allocate(work(iter_used), iwork(iter_used))
    !write(*,*) betamix(:,:)
    CALL DSYTRF('U', iter_used, betamix, iter_used, iwork, work, iter_used, info )
    CALL errore('broyden', 'factorization', abs(info) )
    !
    CALL DSYTRI('U', iter_used, betamix, iter_used, iwork, work, info )
    CALL errore('broyden', 'DSYTRI', abs(info) )    !
    deallocate(iwork)
    !
    FORALL( i = 1:iter_used, &
            j = 1:iter_used, j > i ) betamix(j,i) = betamix(i,j)
    !
    DO i = 1, iter_used
      work(i) = rho_ddot( df(i), rhout_m, ngm0 )
    ENDDO
    !
    DO i = 1, iter_used
      !
      gamma0 = DOT_PRODUCT( betamix(1:iter_used,i), work(1:iter_used) )
      !write(*,*) 'gamma0 = ', gamma0
      !
      call mix_type_AXPY( -gamma0, dv(i), rhoin_m )
      call mix_type_AXPY( -gamma0, df(i), rhout_m )
    ENDDO
    DEALLOCATE(betamix, work)
    !
    ! auxiliary vectors dv and df not needed anymore
    !
  ENDIF skip_on_first
  !
  IF ( ALLOCATED( df ) ) THEN
    DO i=1, n_iter
      call destroy_mix_type(df(i))
    END DO
    DEALLOCATE( df )
  ENDIF
  IF( ALLOCATED( dv ) ) THEN
    DO i=1, n_iter
      call destroy_mix_type(dv(i))
    END DO
    DEALLOCATE( dv )
  ENDIF


  !
  ! preconditioning the new search direction
  !
  write(*,*) 'imix = ', imix
  IF ( imix == 1 ) THEN
    CALL my_approx_screening( rhout_m )
  ELSEIF ( imix == 2 ) THEN
    CALL my_approx_screening2( rhout_m, rhoin_m )
  ENDIF

  !
  ! set new trial density
  !
  ! rhoin_m <- alphamix*rhout_m + rhoin_m
  call mix_type_AXPY( alphamix, rhout_m, rhoin_m )

  ! simple mixing for high_frequencies (and set to zero the smooth ones)
  ! rhoin <- input_rhout
  call high_frequency_mixing( rhoin, input_rhout, alphamix )

  ! add the mixed rho for the smooth frequencies
  ! rhoin_m ->  rhoin
  call assign_mix_to_scf_type(rhoin_m, rhoin)
  !
  call destroy_mix_type(rhout_m)
  call destroy_mix_type(rhoin_m)

  ! Final output is rhoin

  write(*,*)
  write(*,*) '*********** Exit my_mix_rho ************'
  write(*,*)

  RETURN
  !
END SUBROUTINE my_mix_rho



!----------------------------------------------------------------------------
SUBROUTINE my_approx_screening( drho )
!----------------------------------------------------------------------------
  !! Apply an average TF preconditioning to \(\text{drho}\).
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, pi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, ngm
  USE klist,         ONLY : nelec
  USE control_flags, ONLY : ngm0
  USE scf,           ONLY : mix_type
  USE wavefunctions, ONLY : psic
  !
  IMPLICIT NONE
  !
  type (mix_type), intent(INOUT) :: drho ! (in/out)
  !
  REAL(DP) :: rs, agg0

  write(*,*) 'Calling my_approx_screening'
  !
  rs = ( 3.D0 * omega / fpi / nelec )**( 1.D0 / 3.D0 )
  !
  agg0 = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / rs
  !
  drho%of_g(:ngm0,1) =  drho%of_g(:ngm0,1) * gg(:ngm0) / (gg(:ngm0)+agg0)
  !
  RETURN
  !
END SUBROUTINE my_approx_screening



!----------------------------------------------------------------------------
SUBROUTINE my_approx_screening2( drho, rhobest )
!----------------------------------------------------------------------------
  !! Apply a local-density dependent TF preconditioning to \(\text{drho}\).
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : e2, pi, tpi, fpi, eps8, eps32
  USE cell_base,            ONLY : omega, tpiba2
  USE gvect,                ONLY : gg, ngm
  USE wavefunctions, ONLY : psic
  USE klist,                ONLY : nelec
  USE control_flags,        ONLY : ngm0, gamma_only
  USE scf,                  ONLY : mix_type, local_tf_ddot
  USE mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  type(mix_type), intent(inout) :: drho
  type(mix_type), intent(in) :: rhobest
  !
  INTEGER, PARAMETER :: mmx = 12
  !
  INTEGER :: &
    iwork(mmx), i, j, m, info, is
  REAL(DP) :: &
    avg_rsm1, target, dr2_best
  REAL(DP) :: &
    aa(mmx,mmx), invaa(mmx,mmx), bb(mmx), work(mmx), vec(mmx), agg0
  COMPLEX(DP), ALLOCATABLE :: &
    v(:,:),     &! v(ngm0,mmx)
    w(:,:),     &! w(ngm0,mmx)
    dv(:),      &! dv(ngm0)
    vbest(:),   &! vbest(ngm0)
    wbest(:)     ! wbest(ngm0)
  REAL(DP), ALLOCATABLE :: &
    alpha(:)     ! alpha(dffts%nnr)
  !
  INTEGER             :: ir, ig
  REAL(DP), PARAMETER :: one_third = 1.D0 / 3.D0
  !
  target = 0.D0
  !
  IF ( gg(1) < eps8 ) drho%of_g(1,1) = ZERO
  !
  ALLOCATE( alpha( dffts%nnr ) )
  ALLOCATE( v( ngm0, mmx ), &
            w( ngm0, mmx ), dv( ngm0 ), vbest( ngm0 ), wbest( ngm0 ) )
  !
  !$omp parallel
     !
     CALL threaded_barrier_memset(psic, 0.0_DP, dffts%nnr*2)
     !$omp do
     DO ig = 1, ngm0
        psic(dffts%nl(ig)) = rhobest%of_g(ig,1)
     ENDDO
     !$omp end do nowait
     !
  !$omp end parallel
  !
  ! ... calculate alpha from density
  !
  CALL invfft ('Rho', psic, dffts)
  !
  avg_rsm1 = 0.D0
  !
  !$omp parallel do reduction(+:avg_rsm1)
  DO ir = 1, dffts%nnr
     alpha(ir) = ABS( REAL( psic(ir) ) )
     !
     IF ( alpha(ir) > eps32 ) THEN
        !
        alpha(ir) = ( 3.D0 / fpi / alpha(ir) )**one_third
        avg_rsm1  = avg_rsm1 + 1.D0 / alpha(ir)
        !
     END IF   
     !
     alpha(ir) = 3.D0 * ( tpi / 3.D0 )**( 5.D0 / 3.D0 ) * alpha(ir)
     !
  END DO
  !$omp end parallel do
  !
  CALL mp_sum( avg_rsm1 , intra_bgrp_comm )
  avg_rsm1 = ( dffts%nr1*dffts%nr2*dffts%nr3 ) / avg_rsm1
  agg0     = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / avg_rsm1
  !
  ! ... calculate deltaV and the first correction vector
  !
  !$omp parallel
     CALL threaded_barrier_memset(psic, 0.0_DP, dffts%nnr*2)
     !$omp do
     DO ig = 1, ngm0
        psic(dffts%nl(ig)) = drho%of_g(ig,1)
     ENDDO
     !$omp end do nowait
     !
     IF ( gamma_only ) THEN
        !$omp do
        DO ig = 1, ngm0
           psic(dffts%nlm(ig)) = CONJG( psic(dffts%nl(ig)) )
        ENDDO
        !$omp end do nowait
     ENDIF
  !$omp end parallel
  !
  CALL invfft ('Rho', psic, dffts)
  !
  !$omp parallel do
  DO ir = 1, dffts%nnr
     psic(ir) = psic(ir) * alpha(ir)
  ENDDO
  !$omp end parallel do
  !
  CALL fwfft ('Rho', psic, dffts)
  !
  !$omp parallel do
  DO ig = 1, ngm0
     dv(ig) = psic(dffts%nl(ig)) * gg(ig) * tpiba2
     v(ig,1)= psic(dffts%nl(ig)) * gg(ig) / ( gg(ig) + agg0 )
  ENDDO
  !$omp end parallel do
  !
  m       = 1
  aa(:,:) = 0.D0
  bb(:)   = 0.D0
  !
  repeat_loop: DO
     !
     ! ... generate the vector w
     !     
     !$omp parallel
        CALL threaded_barrier_memset(psic, 0.0_DP, dffts%nnr*2)
        !$omp do
        DO ig = 1, ngm0
           !
           w(ig,m) = fpi * e2 * v(ig,m)
           !
           psic(dffts%nl(ig)) = v(ig,m)
        ENDDO
        !$omp end do nowait
        !
        IF ( gamma_only ) THEN
           !$omp do
           DO ig = 1, ngm0
              psic(dffts%nlm(ig)) = CONJG( psic(dffts%nl(ig)) )
           ENDDO
           !$omp end do nowait
        ENDIF
     !$omp end parallel
     !
     CALL invfft ('Rho', psic, dffts)
     !
     !$omp parallel do
     DO ir = 1, dffts%nnr
        psic(ir) = psic(ir) * alpha(ir)
     ENDDO
     !$omp end parallel do
     !
     CALL fwfft ('Rho', psic, dffts)
     !
     !$omp parallel do
     DO ig = 1, ngm0
        w(ig,m) = w(ig,m) + gg(ig) * tpiba2 * psic(dffts%nl(ig))
     ENDDO
     !$omp end parallel do
     !
     ! ... build the linear system
     !
     DO i = 1, m
        !
        aa(i,m) = local_tf_ddot( w(1,i), w(1,m), ngm0)
        !
        aa(m,i) = aa(i,m)
        !
     END DO
     !
     bb(m) = local_tf_ddot( w(1,m), dv, ngm0)
     !
     ! ... solve it -> vec
     !
     invaa = aa
     !
     CALL DSYTRF( 'U', m, invaa, mmx, iwork, work, mmx, info )
     CALL errore( 'broyden', 'factorization', info )
     !
     CALL DSYTRI( 'U', m, invaa, mmx, iwork, work, info )
     CALL errore( 'broyden', 'DSYTRI', info )
     !     
     FORALL( i = 1:m, j = 1:m, j > i ) invaa(j,i) = invaa(i,j)
     !
     FORALL( i = 1:m ) vec(i) = SUM( invaa(i,:)*bb(:) )
     !
     !$omp parallel
        !$omp do
        DO ig = 1, ngm0
           vbest(ig) = ZERO
           wbest(ig) = dv(ig)
        ENDDO
        !$omp end do nowait
        !
        DO i = 1, m
           !$omp do
           DO ig = 1, ngm0
              vbest(ig) = vbest(ig) + vec(i) * v(ig,i)
              wbest(ig) = wbest(ig) - vec(i) * w(ig,i)
           ENDDO
           !$omp end do nowait
        END DO
     !$omp end parallel
     !
     dr2_best = local_tf_ddot( wbest, wbest, ngm0 )
     !
     IF ( target == 0.D0 ) target = MAX( 1.D-12, 1.D-6*dr2_best )
     !
     IF ( dr2_best < target ) THEN
        !
        !$omp parallel
           !$omp do
           DO ig = 1, ngm0
              drho%of_g(ig,1) = vbest(ig)
           ENDDO
           !$omp end do nowait
           !
        !$omp end parallel
        !
        DEALLOCATE( alpha, v, w, dv, vbest, wbest )
        !
        EXIT repeat_loop
        !
     ELSE IF ( m >= mmx ) THEN
        !
        m = 1
        !
        !$omp parallel do
        DO ig = 1, ngm0
           v(ig,m)  = vbest(ig)
        ENDDO
        !$omp end parallel do
        aa(:,:) = 0.D0
        bb(:)   = 0.D0
        !
        CYCLE repeat_loop
        !
     END IF
     !
     m = m + 1
     !
     !$omp parallel do
     DO ig = 1, ngm0
        v(ig,m) = wbest(ig) / ( gg(ig) + agg0 )
     ENDDO
     !$omp end parallel do
     !
  END DO repeat_loop
  !
  RETURN
  !
END SUBROUTINE my_approx_screening2
