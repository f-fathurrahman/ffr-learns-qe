!-----------------------------------------------------------------------
SUBROUTINE my_exx_set_symm( nr1, nr2, nr3, nr1x, nr2x, nr3x )
  !-----------------------------------------------------------------------
  !! Uses \(\text{nkqs}\) and \(\text{index_sym}\) from module \(\texttt{exx}\),
  !! computes \(\text{rir}\).
  !
  use kinds, only: dp
  USE symm_base,  ONLY : nsym, s, ft
  !
  use exx_base, only: rir
  !
  IMPLICIT NONE
  !
  INTEGER :: nr1, nr2, nr3, nr1x, nr2x, nr3x 
  !
  ! local variables
  !
  INTEGER :: isym, i,j,k, ri,rj,rk, ir, nxxs
  INTEGER, allocatable :: ftau(:,:), s_scaled(:,:,:)
  !

  write(*,*) 
  write(*,*) '<div> ENTER my_exx_set_symm'
  write(*,*)
  write(*,*) 'nr1,nr2,nr3 = ', nr1, nr2, nr3
  write(*,*) 'nr1x,nr2x,nr3x = ', nr1x, nr2x, nr3x

  nxxs = nr1x*nr2x*nr3x
  !
  IF (.NOT. ALLOCATED(rir)) THEN   
    ALLOCATE( rir(nxxs,nsym) )
  ELSEIF ((SIZE(rir,1) /= nxxs) ) THEN 
    DEALLOCATE( rir )
    ALLOCATE( rir(nxxs,nsym) )
  ENDIF
  !
  rir = 0 ! rir is an integer array
  ALLOCATE( ftau(3,nsym), s_scaled(3,3,nsym) )
  CALL scale_sym_ops(nsym, s, ft, nr1, nr2, nr3, s_scaled, ftau)
  DO isym = 1, nsym
    DO k = 1, nr3
      DO j = 1, nr2
        DO i = 1, nr1
          CALL rotate_grid_point( s_scaled(1,1,isym), ftau(1,isym), &
              i, j, k, nr1, nr2, nr3, ri, rj, rk )
          ir = i + (j-1)*nr1x + (k-1)*nr1x*nr2x
          rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  DEALLOCATE ( s_scaled, ftau )

  write(*,*) 'Some rir:'
  write(*,*) 'rir(1:4,1) = ', rir(1:4,1)
  write(*,*) 'rir(100:104,2) = ', rir(100:104,2)

  !
  write(*,*) 
  write(*,*) '</div> EXIT my_exx_set_symm'
  write(*,*)
  flush(0) ! flush stderr
  !stop 'Early stop 63 in my_exx_set_symm'
  !
  return
  !
END SUBROUTINE