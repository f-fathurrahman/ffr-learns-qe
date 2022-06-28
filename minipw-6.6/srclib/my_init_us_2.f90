!----------------------------------------------------------------------
SUBROUTINE my_init_us_2( npw_, igk_, q_, vkb_ )
!----------------------------------------------------------------------
  !! Calculates beta functions (Kleinman-Bylander projectors), with
  !! structure factor, for all atoms, in reciprocal space.
  !
  USE kinds,        ONLY : DP
  USE ions_base,    ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,    ONLY : tpiba, omega
  USE constants,    ONLY : tpi
  USE gvect,        ONLY : eigts1, eigts2, eigts3, mill, g
  USE wvfct,        ONLY : npwx
  USE us,           ONLY : dq, tab, spline_ps
  USE m_gth,        ONLY : mk_ffnl_gth
  USE splinelib
  USE uspp,         ONLY : nkb, nhtol, nhtolm, indv
  USE uspp_param,   ONLY : upf, lmaxkb, nhm, nh
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw_
  !! number of PWs 
  INTEGER, INTENT(IN) :: igk_(npw_)
  !! indices of G in the list of q+G vectors
  REAL(DP), INTENT(IN) :: q_(3)
  !! q vector (2pi/a units)
  COMPLEX(DP), INTENT(OUT) :: vkb_(npwx,nkb)
  !! beta functions (npw_ <= npwx)
  !
  ! ... Local variables
  !
  INTEGER :: i0, i1, i2, i3, ig, lm, na, nt, nb, ih, jkb
  REAL(DP) :: px, ux, vx, wx, arg
  REAL(DP), ALLOCATABLE :: gk(:,:), qg(:), vq(:), ylm(:,:), vkb1(:,:)
  COMPLEX(DP) :: phase, pref
  COMPLEX(DP), ALLOCATABLE :: sk(:)

  !write(*,*)
  !write(*,*) '*** Enter my_init_us_2'
  !write(*,*)

  IF (lmaxkb < 0) RETURN

  !
  IF (spline_ps) THEN
    stop 'spline_ps is not supported'
  ENDIF

  !
  ALLOCATE( vkb1(npw_,nhm) )
  ALLOCATE( sk(npw_) )
  ALLOCATE( qg(npw_) )
  ALLOCATE( vq(npw_) )
  ALLOCATE( ylm(npw_,(lmaxkb+1)**2) )
  ALLOCATE( gk(3,npw_) )

  DO ig = 1, npw_
    gk(1,ig) = q_(1) + g(1,igk_(ig) )
    gk(2,ig) = q_(2) + g(2,igk_(ig) )
    gk(3,ig) = q_(3) + g(3,igk_(ig) )
    qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  ENDDO

  !CALL ylmr2( (lmaxkb+1)**2, npw_, gk, qg, ylm(1:npw_,:) )
  CALL ylmr2( (lmaxkb+1)**2, npw_, gk, qg, ylm )
  
  !
  ! set now qg=|q+G| in atomic units
  !
  DO ig = 1, npw_
     qg(ig) = SQRT(qg(ig))*tpiba
  ENDDO
  
  !
  ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
  jkb = 0
  DO nt = 1, ntyp
    ! ... calculate beta in G-space using an interpolation table:
    !     f_l(q)=\int _0 ^\infty dr r^2 f_l(r) j_l(q.r)
    DO nb = 1, upf(nt)%nbeta
      !
      IF ( upf(nt)%is_gth ) THEN
        CALL mk_ffnl_gth( nt, nb, npw_, omega, qg, vq )
      ELSE
        DO ig = 1, npw_
          ! qg: (G+k)^2
          !write(*,*)
          !write(*,*) 'Pass here in init_us_2'
          px = qg(ig) / dq - INT( qg(ig)/dq )
          ux = 1.d0 - px
          vx = 2.d0 - px
          wx = 3.d0 - px
          i0 = INT( qg(ig)/dq ) + 1
          i1 = i0 + 1
          i2 = i0 + 2
          i3 = i0 + 3
          vq(ig) = tab(i0,nb,nt) * ux * vx * wx / 6.d0 + &
                   tab(i1,nb,nt) * px * vx * wx / 2.d0 - &
                   tab(i2,nb,nt) * px * ux * wx / 2.d0 + &
                   tab(i3,nb,nt) * px * ux * vx / 6.d0
        ENDDO
      ENDIF
      
      ! add spherical harmonic part  (Y_lm(q)*f_l(q)) 
      DO ih = 1, nh(nt)
        IF( nb == indv(ih,nt) ) THEN
          !l = nhtol(ih, nt)
          lm = nhtolm(ih,nt)
          DO ig = 1, npw_
            vkb1(ig,ih) = ylm(ig,lm) * vq(ig)
          ENDDO
        ENDIF
      ENDDO

    ENDDO ! loop over nb

    !
    ! vkb1 contains all betas including angular part for type nt
    ! now add the structure factor and factor (-i)^l
    !
    DO na = 1, nat
      !
      ! ordering: first all betas for atoms of type 1
      !           then  all betas for atoms of type 2  and so on
      !
      IF(ityp(na) == nt) THEN
        !
        arg = ( q_(1) * tau(1,na) + &
                q_(2) * tau(2,na) + &
                q_(3) * tau(3,na) ) * tpi
        phase = CMPLX( COS(arg), -SIN(arg) ,KIND=DP)
        !
        DO ig = 1, npw_
          sk(ig) = eigts1(mill(1,igk_(ig)), na) * &
                   eigts2(mill(2,igk_(ig)), na) * &
                   eigts3(mill(3,igk_(ig)), na)
        ENDDO
        !
        DO ih = 1, nh(nt)
          ! XXXX idx of KB projectors increment here ...
          jkb = jkb + 1
          !
          pref = (0.d0, -1.d0)**nhtol(ih, nt) * phase
          DO ig = 1, npw_
            vkb_(ig, jkb) = vkb1(ig,ih) * sk(ig) * pref
          ENDDO
          ! clean up garbage in the last block
          DO ig = npw_+1, npwx
            vkb_(ig, jkb) = (0.0_DP, 0.0_DP)
          ENDDO
        ENDDO
          !
      ENDIF
       !
    ENDDO

  ENDDO

  DEALLOCATE( gk )
  DEALLOCATE( ylm )
  DEALLOCATE( vq )
  DEALLOCATE( qg )
  DEALLOCATE( sk )
  DEALLOCATE( vkb1 )

  !stop 'ffr my_init_us_2 185'

  !
  RETURN
  !
END SUBROUTINE

