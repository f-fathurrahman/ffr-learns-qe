!-----------------------------------------------------------------------
  FUNCTION my_exx_divergence()
!-----------------------------------------------------------------------
  use kinds, only: dp
  USE constants, ONLY : fpi, e2, pi
  USE cell_base, ONLY : bg, at, alat, omega
  USE gvect, ONLY : ngm, g
  USE gvecw, ONLY : gcutw
  USE control_flags, ONLY : gamma_only
  USE mp, ONLY : mp_sum
  use gvecw, only: ecutwfc
  use exx_base, only: use_regularization, x_gamma_extrapolation, nq1, nq2, nq3, eps, &
                      erfc_scrlen, erf_scrlen, yukawa, on_double_grid, nqs, grid_factor
  !
  IMPLICIT NONE
  !
  REAL(DP) :: my_exx_divergence
  !
  ! local variables
  !
  INTEGER :: iq1,iq2,iq3, ig
  REAL(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, &
              tpiba2, alpha, x, q(3)
  INTEGER :: nqq, iq
  REAL(DP) :: aa, dq

  write(*,*)
  write(*,*) '<div> ENTER my_exx_divergence'
  write(*,*)

  !
  tpiba2 = (fpi / 2.d0 / alat)**2
  !
  alpha  = 10._dp / gcutw
  write(*,*) 'ecutwfc = ', ecutwfc
  write(*,*) 'alat = ', alat
  write(*,*) 'tpiba2 = ', tpiba2
  write(*,*) 'gcutw = ', gcutw
  write(*,*) 'alpha = ', alpha
  write(*,*) 'alat calc = ', SQRT( at(1,1)**2 + at(2,1)**2 + at(3,1)**2 )

  !
  IF ( .NOT. use_regularization ) THEN
    my_exx_divergence = 0._dp
    RETURN
  ENDIF
  !
  dq1 = 1._dp / DBLE(nq1)
  dq2 = 1._dp / DBLE(nq2) 
  dq3 = 1._dp / DBLE(nq3) 
  !
  div = 0._dp
  !
  DO iq1 = 1, nq1
    DO iq2 = 1, nq2
      DO iq3 = 1, nq3
        xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                bg(:,2) * (iq2-1) * dq2 + &
                bg(:,3) * (iq3-1) * dq3
        !
        DO ig = 1, ngm
          !
          q(1) = xq(1) + g(1,ig)
          q(2) = xq(2) + g(2,ig)
          q(3) = xq(3) + g(3,ig)
          qq = ( q(1)**2 + q(2)**2 + q(3)**2 )
          !
          IF (x_gamma_extrapolation) THEN
            on_double_grid = .TRUE.
            x = 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
            on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
            x = 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
            on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
            x = 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
            on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
          ENDIF
          !
          IF (.NOT.on_double_grid) THEN
            IF ( qq > 1.d-8 ) THEN
              IF ( erfc_scrlen > 0 ) THEN
                div = div + EXP( -alpha * qq) / qq * (1._dp-EXP(-qq*tpiba2/4.d0/erfc_scrlen**2)) * grid_factor
              ELSEIF ( erf_scrlen >0 ) THEN
                div = div + EXP( -alpha * qq) / qq * (EXP(-qq*tpiba2/4.d0/erf_scrlen**2)) * grid_factor
              ELSE
                div = div + EXP( -alpha * qq) / (qq + yukawa/tpiba2) * grid_factor
              ENDIF
              !write(*,*) ig, qq, div
            ENDIF
          ENDIF
        ENDDO ! ig
      ENDDO ! 
    ENDDO
  ENDDO
  write(*,*) 'Line 94: div = ', div
  !
  IF (gamma_only) THEN
    div = 2.d0 * div
  ENDIF
  !
  IF ( .NOT. x_gamma_extrapolation ) THEN
    IF ( yukawa > 0._dp) THEN
      div = div + tpiba2/yukawa
    ELSEIF( erfc_scrlen > 0._dp ) THEN
      div = div + tpiba2/4.d0/erfc_scrlen**2
    ELSE
      div = div - alpha
    ENDIF
  ENDIF
  write(*,*) 'Line 114 div = ', div
  !
  div = div * e2 * fpi / tpiba2 / nqs ! 4*pi/(4pi^2 / alat)
  write(*,*) 'e2 = ', e2
  write(*,*) 'fpi / tpiba2 = ', fpi / tpiba2
  write(*,*) 'alat^2/pi = ', alat**2/pi
  write(*,*) '(in Ry) div = ', div ! in 
  !
  alpha = alpha / tpiba2
  !
  nqq = 100000
  dq = 5.0d0 / SQRT(alpha) / nqq
  aa = 0._dp
  !
  DO iq = 0, nqq
    q_ = dq * (iq+0.5d0)
    qq = q_ * q_
    IF ( erfc_scrlen > 0 ) THEN
        aa = aa  -EXP( -alpha * qq) * EXP(-qq/4.d0/erfc_scrlen**2)*dq
    ELSEIF ( erf_scrlen > 0 ) THEN
        aa = 0._dp
    ELSE
        aa = aa - EXP( -alpha * qq) * yukawa / (qq + yukawa)*dq
    ENDIF
  ENDDO
  write(*,*) 'aa = ', aa
  !
  aa = aa * 8.d0/fpi
  aa = aa + 1._dp/SQRT(alpha*0.25d0*fpi)
  write(*,*) 'aa = ', aa
  IF ( erf_scrlen > 0) aa = 1._dp/SQRT((alpha+1._dp/4.d0/erf_scrlen**2)*0.25d0*fpi)
  div = div - e2*omega * aa
  !
  my_exx_divergence = div * nqs

  write(*,*)
  write(*,*) '</div> EXIT my_exx_divergence'
  write(*,*)
  !
  RETURN
END FUNCTION
