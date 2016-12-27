SUBROUTINE info_structure()
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg, alat
  IMPLICIT NONE 
  INTEGER :: i
  REAL(DP) :: ab(3,3)

  WRITE(*,*)
  WRITE(*,*) 'INFO_STRUCTURE'
  WRITE(*,*) '--------------'
  WRITE(*,*)

  WRITE(*,'(1x,A,F18.10)') 'alat = ', alat

  WRITE(*,*) 'at = '
  DO i = 1,3
    WRITE(*,'(1x,3F18.10)') at(i,1:3)
  ENDDO

  WRITE(*,*) 'bg = '
  DO i = 1,3
    WRITE(*,'(1x,3F18.10)') bg(i,1:3)
  ENDDO

  ab = matmul( at, transpose(bg) )
  WRITE(*,*) 'at * bg'' = '
  DO i = 1,3
    WRITE(*,'(1x,3F18.10)') ab(i,1:3)
  ENDDO

END SUBROUTINE 

