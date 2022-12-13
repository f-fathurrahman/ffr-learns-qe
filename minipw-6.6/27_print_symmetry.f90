INCLUDE 'prepare_all.f90'

PROGRAM main
  USE kinds, ONLY: DP
  IMPLICIT NONE 

  CALL prepare_all()
  
  CALL print_symmetry()

END PROGRAM main


!--------------------------
SUBROUTINE print_symmetry()
!--------------------------
  USE symm_base, ONLY: sr, Nsym, d1, d2, d3
  IMPLICIT NONE
  integer :: i
  
  WRITE(*,*) 'Nsym = ', Nsym

  WRITE(*,*)
  WRITE(*,*) 'sr matrices'
  DO i = 1,Nsym
    WRITE(*,*)
    WRITE(*,*) 'isym = ', i
    CALL print_matrix_r8(sr(:,:,i), 3, 3)
  ENDDO

  WRITE(*,*)
  WRITE(*,*) 'dy1, dy2, dy3 matrices'
  DO i = 1,Nsym
    WRITE(*,*)
    WRITE(*,*) 'isym = ', i
    WRITE(*,*)
    WRITE(*,*) 'dy1 = '
    CALL print_matrix_r8(d1(:,:,i), 3, 3)
    WRITE(*,*)
    WRITE(*,*) 'dy2 = ' 
    CALL print_matrix_r8(d2(:,:,i), 5, 5)
    WRITE(*,*)
    WRITE(*,*) 'dy3 = '
    CALL print_matrix_r8(d3(:,:,i), 7, 7)
  ENDDO
  
END SUBROUTINE print_symmetry


!----------------------------------------
SUBROUTINE print_matrix_r8(A, Nrows, Ncols)
!----------------------------------------
  IMPLICIT NONE
  INTEGER :: Nrows, Ncols
  REAL(8) :: A(Nrows, Ncols)
  INTEGER :: i, j

  DO i = 1,Nrows
    DO j = 1,Ncols
      WRITE(*,'(1x,F10.5)',advance='no') A(i,j)
    ENDDO
    WRITE(*,*)
  ENDDO

  RETURN
  
END SUBROUTINE print_matrix_r8


