! Set only one k-point, i.e. (0, 0, 0) kpoint for the calculation.
! Note that, this is different from `gamma_only` where gamma-point
! trick is used.
!------------------------------------------------------------------------------
SUBROUTINE setup_kpoints_000()
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE parameters, ONLY : npk
  USE cell_base, ONLY : bg
  USE klist, ONLY : nks, wk, xk, nkstot
  IMPLICIT NONE
  !
  INTEGER :: ik, ii

  ! k-points for band structure calculation
  ! TODO: read from file
  nks = 1
  nkstot = nks
  wk(1:nks) = 1.0_DP  ! it is not really used in band structure calculation
  xk(:,1:nks) = reshape( (/ 0.00000000, 0.00000000, 0.00000000 /), (/ 3,nks /) )

  CALL cryst_to_cart( nkstot, xk, bg, 1 )
  IF(ionode) THEN
    WRITE(stdout,*) 'npk, nks = ', npk, nks
    WRITE(stdout, '(23x,"cart. coord. in units 2pi/alat")')
    DO ik=1,nkstot
      WRITE(stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') ik, &
        ( xk(ii,ik), ii=1,3), wk(ik)
    ENDDO
  ENDIF

END SUBROUTINE
