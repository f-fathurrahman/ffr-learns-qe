! this will be called by my_exx_base_exx_grid_init
!----------------------------------------------------------------------
SUBROUTINE my_exx_base_exx_qgrid_init(temp_nkqs, xk_collect, temp_xkq, nkqs, temp_index_ikq, dxk)
!------------------------------------------------------------------------
  use kinds, only : dp
  !! Generate q-point mesh compatible with the k-point mesh
  !
  USE klist,     ONLY : nkstot
  USE cell_base, ONLY : at
  USE symm_base, ONLY : nsym
  !
  use exx_base, only: index_xkq, nq1, nq2, nq3, eps, nqs
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: temp_nkqs
  REAL(DP), INTENT(IN) :: xk_collect(:,:), temp_xkq(:,:)
  !
  REAL(DP), INTENT(OUT) :: dxk(:)
  INTEGER, INTENT(OUT) :: temp_index_ikq(:)
  INTEGER, INTENT(OUT) :: nkqs
  !
  INTEGER :: ik, ikq, iq, iq1, iq2, iq3, j, max_nk
  INTEGER, ALLOCATABLE :: new_ikq(:)
  REAL(DP) :: sxk(3), xk_cryst(3), dq1, dq2, dq3
  LOGICAL :: xk_not_found
  !
  max_nk = nkstot * MIN(48, 2 * nsym)
  ALLOCATE( new_ikq(max_nk) )

  IF ( ALLOCATED(index_xkq) ) DEALLOCATE( index_xkq )
  ALLOCATE( index_xkq(nkstot, nqs) )
  !
  nkqs = 0
  new_ikq(:) = 0
  !
  ! define the q-mesh step-sizes
  !
  dq1 = 1._dp / DBLE(nq1)
  dq2 = 1._dp / DBLE(nq2)
  dq3 = 1._dp / DBLE(nq3)
  !
  DO ik = 1, nkstot
    ! go to crystalline coordinates
    xk_cryst(:) = xk_collect(:,ik)
    CALL cryst_to_cart( 1, xk_cryst, at, -1 )
    !
    iq = 0
    !
    DO iq1 = 1, nq1
      sxk(1) = xk_cryst(1) + (iq1-1) * dq1
      DO iq2 = 1, nq2
        sxk(2) = xk_cryst(2) + (iq2-1) * dq2
        DO iq3 = 1, nq3
            sxk(3) = xk_cryst(3) + (iq3-1) * dq3
            iq = iq + 1
            xk_not_found = .TRUE.
            !
            DO ikq = 1, temp_nkqs
              IF ( xk_not_found ) THEN
                  dxk(:) = sxk(:)-temp_xkq(:,ikq) - NINT(sxk(:)-temp_xkq(:,ikq))
                  IF ( ALL(ABS(dxk) < eps ) ) THEN
                      xk_not_found = .FALSE.
                      IF ( new_ikq(ikq) == 0) THEN
                          nkqs = nkqs + 1
                          temp_index_ikq(nkqs) = ikq
                          new_ikq(ikq) = nkqs
                      ENDIF
                      index_xkq(ik,iq) = new_ikq(ikq)
                  ENDIF
              ENDIF
            ENDDO ! ikq
            !
            IF (xk_not_found) THEN
              DEALLOCATE( new_ikq )
              RETURN
            ENDIF
            !
        ENDDO
      ENDDO
    ENDDO
    !
  ENDDO
  !
  DEALLOCATE( new_ikq )
END SUBROUTINE my_exx_base_exx_qgrid_init