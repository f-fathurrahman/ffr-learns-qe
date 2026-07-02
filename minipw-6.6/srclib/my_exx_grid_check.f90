!------------------------------------------------------------------------
SUBROUTINE my_exx_grid_check( xk_collect, ld2 )
!------------------------------------------------------------------------
  !
  use kinds, only: dp
  !
  USE symm_base, ONLY : s
  USE cell_base, ONLY : at
  USE klist, ONLY : nkstot
  !
  use exx_base, only: nq1, nq2, nq3, index_xk, index_xkq, index_sym, eps
  !
  IMPLICIT NONE
  !
  integer :: ld2
  REAL(DP) :: xk_collect(3,ld2)
  !!
  !
  ! local variables
  !
  REAL(DP) :: sxk(3), dxk(3), xk_cryst(3), xkk_cryst(3)
  INTEGER :: iq1, iq2, iq3, isym, ik, ikk, ikq, iq
  REAL(DP) :: dq1, dq2, dq3
  
  write(*,*)
  write(*,*) '<div> ENTER my_exx_grid_check'
  write(*,*)
  !
  dq1 = 1.0_dp/DBLE(nq1)
  dq2 = 1.0_dp/DBLE(nq2)
  dq3 = 1.0_dp/DBLE(nq3)
  write(*,*) 'dq1,dq2,dq3 = ', dq1, dq2, dq3
  !
  DO ik = 1, nkstot
    write(*,*)
    write(*,*) 'Begin ik = ', ik
    xk_cryst(:) = xk_collect(:,ik)
    CALL cryst_to_cart( 1, xk_cryst, at, -1 )
    write(*,*) 'xk_cryst = ', xk_cryst
    !
    iq = 0
    DO iq1 = 1, nq1
      sxk(1) = xk_cryst(1) + (iq1-1) * dq1
      DO iq2 = 1, nq2
        sxk(2) = xk_cryst(2) + (iq2-1) * dq2
        DO iq3 = 1, nq3
          sxk(3) = xk_cryst(3) + (iq3-1) * dq3
          iq = iq + 1
          !
          ikq  = index_xkq(ik,iq)
          ikk  = index_xk(ikq)
          isym = index_sym(ikq)
          !
          xkk_cryst(:) = at(1,:)*xk_collect(1,ikk) + &
                         at(2,:)*xk_collect(2,ikk) + &
                         at(3,:)*xk_collect(3,ikk)
          write(*,*) 'xkk_cryst = ', xkk_cryst
          IF(isym < 0 ) THEN
            xkk_cryst(:) = -xkk_cryst(:)
          ENDIF
          isym = ABS(isym)
          dxk(:) = s(:,1,isym)*xkk_cryst(1) + &
                   s(:,2,isym)*xkk_cryst(2) + &
                   s(:,3,isym)*xkk_cryst(3) - sxk(:)
          write(*,'(1x,A,3I4,3F18.10)') 'Before nint: ', iq1, iq2, iq3, dxk
          dxk(:) = dxk(:) - NINT(dxk(:))
          write(*,'(1x,A,3I4,3F18.10)') 'After nint: ', iq1, iq2, iq3, dxk
          IF ( .NOT. ( ABS(dxk(1)) <= eps .AND. &
                       ABS(dxk(2)) <= eps .AND. &
                       ABS(dxk(3)) <= eps )   ) THEN
            WRITE(*,*) ik,iq
            WRITE(*,*) ikq,ikk,isym
            WRITE(*,*) dxk(:)
            CALL errore( 'exx_grid_check', 'something wrong', 1 )
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  write(*,*)
  write(*,*) '</div> EXIT my_exx_grid_check'
  write(*,*)

  RETURN
END SUBROUTINE

