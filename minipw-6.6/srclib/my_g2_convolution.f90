
!-----------------------------------------------------------------------
SUBROUTINE my_g2_convolution( ngm, g, xk, xkq, fac )
!----------------------------------------------------------------------
  !! This routine calculates the 1/|r-r'| part of the exact exchange
  !! expression in reciprocal space (the G^-2 factor).
  !! It then regularizes it according to the specified recipe.
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : tpiba, at, tpiba2
  USE constants, ONLY : fpi, e2, pi
  !
  use exx_base, only: gau_scrlen, use_coulomb_vcut_ws, vcut, use_coulomb_vcut_spheric, &
                    & nq1, nq2, nq3, &
                    & x_gamma_extrapolation, yukawa, &
                    & exxdiv, erfc_scrlen, grid_factor, erf_scrlen, &
                    & eps, eps_qdiv
  use coulomb_vcut_module, ONLY : vcut_init, vcut_type, vcut_info, &
                                & vcut_get, vcut_spheric_get
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: ngm !! Number of G vectors
  REAL(DP), INTENT(IN) :: g(3,ngm) !! Cartesian components of G vectors
  REAL(DP), INTENT(IN) :: xk(3) !! current k vector
  REAL(DP), INTENT(IN) :: xkq(3) !! current q vector
  REAL(DP), INTENT(INOUT) :: fac(ngm) !! Calculated convolution
  !
  ! local variables
  !
  INTEGER :: ig !Counters
  REAL(DP) :: q(3), qq, x
  REAL(DP) :: grid_factor_track(ngm), qq_track(ngm)
  REAL(DP) :: nqhalf_dble(3)
  LOGICAL :: odg(3)
  !
  ! First the types of Coulomb potential that need q(3) and an external call
  IF( use_coulomb_vcut_ws ) THEN
    DO ig = 1, ngm
      q(:)= ( xk(:) - xkq(:) + g(:,ig) ) * tpiba
      fac(ig) = vcut_get(vcut,q)
    ENDDO
    RETURN
  ENDIF
  !
  IF( use_coulomb_vcut_spheric ) THEN
    DO ig = 1, ngm
      q(:)= ( xk(:) - xkq(:) + g(:,ig) ) * tpiba
      fac(ig) = vcut_spheric_get(vcut,q)
    ENDDO
    RETURN
  ENDIF
  !
  ! Now the Coulomb potential that are computed on the fly
  !
  nqhalf_dble(1:3) = (/ DBLE(nq1)*0.5_DP, DBLE(nq2)*0.5_DP, DBLE(nq3)*0.5_DP /)
  !
  ! Set the grid_factor_track and qq_track
  !
  IF ( x_gamma_extrapolation ) THEN
    DO ig = 1, ngm
      q(:)= xk(:) - xkq(:) + g(:,ig)
      qq_track(ig) = SUM(q(:)**2) * tpiba2
      x = (q(1)*at(1,1) + q(2)*at(2,1) + q(3)*at(3,1))*nqhalf_dble(1)
      odg(1) = ABS(x-NINT(x)) < eps
      x = (q(1)*at(1,2) + q(2)*at(2,2) + q(3)*at(3,2))*nqhalf_dble(2)
      odg(2) = ABS(x-NINT(x)) < eps
      x = (q(1)*at(1,3) + q(2)*at(2,3) + q(3)*at(3,3))*nqhalf_dble(3)
      odg(3) = ABS(x-NINT(x)) < eps
      IF( ALL( odg(:) ) ) THEN
        grid_factor_track(ig) = 0._DP ! on double grid
      ELSE
        grid_factor_track(ig) = grid_factor ! not on double grid
      ENDIF
    ENDDO
  ELSE
    !ffr: This is the usual one, loop over G-vectors
    DO ig = 1, ngm
      q(:) = xk(:) - xkq(:) + g(:,ig) !ffr: Is this G+q ?
      qq_track(ig) = SUM(q(:)**2) * tpiba2
    ENDDO
    grid_factor_track = 1._DP
  ENDIF
  !
  ! The big loop
  !
  DO ig = 1, ngm
    qq = qq_track(ig)
    IF(gau_scrlen > 0) THEN
      fac(ig) = e2*((pi/gau_scrlen)**(1.5_DP))*EXP(-qq/4._DP/gau_scrlen) * grid_factor_track(ig)
    ELSEIF( qq > eps_qdiv ) THEN
      IF( erfc_scrlen > 0  ) THEN
        fac(ig) = e2*fpi/qq*(1._DP-EXP(-qq/4._DP/erfc_scrlen**2)) * grid_factor_track(ig)
      ELSEIF( erf_scrlen > 0 ) THEN
        fac(ig) = e2*fpi/qq*(EXP(-qq/4._DP/erf_scrlen**2)) * grid_factor_track(ig)
      ELSE
        fac(ig) = e2*fpi/( qq + yukawa ) * grid_factor_track(ig) ! as HARTREE
      ENDIF
      !
    ELSE
      !
      fac(ig) = -exxdiv ! or rather something else (see F.Gygi)
      IF ( yukawa > 0._DP .AND. .NOT. x_gamma_extrapolation) THEN
        fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
      ENDIF
      !
      IF( erfc_scrlen > 0._DP .AND. .NOT. x_gamma_extrapolation) THEN
        fac(ig) = fac(ig) + e2*pi/(erfc_scrlen**2)
      ENDIF
    ENDIF
  ENDDO
END SUBROUTINE my_g2_convolution


!-----------------------------------------------------------------------
SUBROUTINE my_g2_convolution_all( ngm, g, xk, xkq, iq, current_k )
!-----------------------------------------------------------------------
  !! Wrapper for g2_convolution.
  !
  USE kinds, ONLY : DP
  USE klist, ONLY : nks
  !
  use exx_base, only: coulomb_done, nqs, coulomb_fac
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: ngm
  !! Number of G vectors
  REAL(DP), INTENT(IN) :: g(3,ngm)
  !! Cartesian components of G vectors
  REAL(DP), INTENT(IN) :: xk(3)
  !! current k vector
  REAL(DP), INTENT(IN) :: xkq(3)
  !! current q vector
  INTEGER, INTENT(IN) :: current_k
  !! current k-point index
  INTEGER, INTENT(IN) :: iq
  !! q-grid point index
  !
  ! Check if coulomb_fac has been allocated
  IF( .NOT. ALLOCATED( coulomb_fac ) ) ALLOCATE( coulomb_fac(ngm,nqs,nks) )
  !
  ! Check if coulomb_done has been allocated
  IF( .NOT. ALLOCATED( coulomb_done) ) THEN
    ALLOCATE( coulomb_done(nqs,nks) )
    coulomb_done = .FALSE.
  ENDIF
  !
  ! return if this k and k' already computed, otherwise compute it
  IF ( coulomb_done(iq,current_k) ) RETURN
  !
  CALL my_g2_convolution( ngm, g, xk, xkq, coulomb_fac(:,iq,current_k) )
  !
  coulomb_done(iq,current_k) = .TRUE.
  !
END SUBROUTINE my_g2_convolution_all