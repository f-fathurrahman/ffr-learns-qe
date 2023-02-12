INCLUDE 'prepare_all.f90'


PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL test_PAW_symmetrize()
END PROGRAM



!----------------------------------
SUBROUTINE test_PAW_symmetrize()
!----------------------------------
  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat, nsp
  USE uspp_param, ONLY : nh, nhm
  USE lsda_mod, ONLY : nspin
  USE paw_symmetry, only : PAW_symmetrize
  IMPLICIT NONE
  INTEGER :: i,dim1
  REAL(8), ALLOCATABLE :: becsum(:,:,:)

  WRITE(*,*) 'nat = ', nat
  WRITE(*,*) 'nsp = ', nsp
  WRITE(*,*) 'nhm = ', nhm

  dim1 = nhm*(nhm+1)/2
  ALLOCATE( becsum(nhm*(nhm+1)/2, nat, nspin) )
  becsum(:,:,:) = 0.d0
  
  WRITE(*,*) 'shape becsum = ', SHAPE(becsum)

  becsum(1:5,1,1) = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0 /)
  becsum(1:5,2,1) = (/ 8.d0, 8.d0, 8.d0, 8.d0, 8.d0 /)

  WRITE(*,*)
  WRITE(*,*) 'Before PAW_symmetrize: sum(becsum) = ', SUM(becsum)
  DO i = 1,5
    WRITE(*,'(1x,A,F18.10)') 'becsum(i,1,1) = ', becsum(i,1,1)
  ENDDO
  WRITE(*,*)
  DO i = 1,5
    WRITE(*,'(1x,A,F18.10)') 'becsum(i,2,1) = ', becsum(i,2,1)
  ENDDO
  
  CALL PAW_symmetrize(becsum)

  WRITE(*,*)
  WRITE(*,*) 'After PAW_symmetrize: sum(becsum) = ', SUM(becsum)  
 
  DO i = 1,5
    WRITE(*,'(1x,A,F18.10)') 'becsum(i,1,1) = ', becsum(i,1,1)
  ENDDO
  WRITE(*,*)
  DO i = 1,5
    WRITE(*,'(1x,A,F18.10)') 'becsum(i,2,1) = ', becsum(i,2,1)
  ENDDO

  DEALLOCATE( becsum )
  
END SUBROUTINE test_PAW_symmetrize

