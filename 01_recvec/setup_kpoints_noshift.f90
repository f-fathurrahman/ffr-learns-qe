!------------------------------------------------------------------------------
SUBROUTINE setup_kpoints_noshift( nk1_, nk2_, nk3_ )
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE parameters, ONLY : npk
  USE cell_base, ONLY : at, bg
  USE klist, ONLY : nks, wk, xk, nkstot
  USE symm_base, ONLY : nrot, s, t_rev
  ! This is actually can be bypassed. However, we include it and set the
  ! relevant variables in anticipation that they will be accessed by
  ! other subroutines in PWSCF
  USE start_k, ONLY : nk1,nk2,nk3, k1,k2,k3
  IMPLICIT NONE
  !
  INTEGER :: nk1_, nk2_, nk3_
  LOGICAL :: time_reversal, skip_equivalence
  INTEGER :: ik, ii

  time_reversal    = .TRUE.
  skip_equivalence = .FALSE. ! set to FALSE to reduce number of k-points

  ! k1,k2,k3 are the offset, the actual grid sampling are given by nk1,nk2,nk3
  k1 = 0; nk1 = nk1_
  k2 = 0; nk2 = nk2_
  k3 = 0; nk3 = nk3_

  CALL kpoint_grid( 1, time_reversal, skip_equivalence, s, t_rev, bg, &
    npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)

  nkstot = nks

  IF(ionode) THEN
    WRITE(stdout,*) 'npk, nks = ', npk, nks
    WRITE(stdout, '(23x,"cart. coord. in units 2pi/alat")')
    DO ik=1,nkstot
      WRITE(stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') ik, &
        ( xk(ii,ik), ii=1,3), wk(ik)
    ENDDO
  ENDIF

END SUBROUTINE
