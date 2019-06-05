include "my_symme.f90"


PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL test_symmetry()
END PROGRAM 

SUBROUTINE my_symvector(nat, vect)
  USE kinds, ONLY: DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : s, nsym, irt


     !-----------------------------------------------------------------------
     ! Symmetrize a function f(i,na), i=cartesian component, na=atom index
     ! e.g. : forces (in cartesian axis) 
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     REAL(DP), intent(INOUT) :: vect(3,nat)
     !
     INTEGER :: na, isym, nar
     REAL(DP) :: dv(3)
     REAL(DP), ALLOCATABLE :: work (:,:)
     !
     IF (nsym == 1) RETURN
     !
     ALLOCATE (work(3,nat))
     !
     ! bring vector to crystal axis
     !
     DO na = 1, nat
        work(:,na) = vect(1,na)*at(1,:) + &
                     vect(2,na)*at(2,:) + &
                     vect(3,na)*at(3,:)
     END DO
     !
     ! symmetrize in crystal axis
     !
     vect (:,:) = 0.0_dp
     DO na = 1, nat
        WRITE(*,*)
        WRITE(*,*) 'before:', work(:,na)
        DO isym = 1, nsym
           nar = irt (isym, na)
           dv(:) = s(:,1,isym)*work(1,nar) + &
                   s(:,2,isym)*work(2,nar) + &
                   s(:,3,isym)*work(3,nar)
           vect(:, na) = vect(:,na) + dv(:)
           WRITE(*,*) 'dv = ', dv(:)
        END DO
        WRITE(*,*) 'after:', vect(:,na)
     END DO
     work (:,:) = vect (:,:) / DBLE(nsym)
     !
     ! bring vector back to cartesian axis
     !
     DO na = 1, nat
        vect(:,na) = work(1,na)*bg(:,1) + &
                     work(2,na)*bg(:,2) + &
                     work(3,na)*bg(:,3)
     END DO
     !
     DEALLOCATE (work)
     !
END SUBROUTINE



SUBROUTINE test_symmetry()
  
  USE symm_base, ONLY : nrot, fft_fact, ft, ftau, s, nsym, irt
  USE symme, ONLY: symvector
  USE gvect, ONLY : ngm, g
  USE my_symme, ONLY: Ngs, my_sym_rho_init_shells
  USE cell_base, ONLY: at, alat, tpiba

  IMPLICIT NONE 
  INTEGER :: isym, ia
  LOGICAL :: non_symmorphic(48)
  REAL(8) :: ft_(3,48), sg(3), arg
  REAL(8) :: v(3,2)

  WRITE(*,*)
  WRITE(*,*) 'nrot = ', nrot

  WRITE(*,*)
  WRITE(*,*) 'fft_fact(1) = ', fft_fact(1)
  WRITE(*,*) 'fft_fact(2) = ', fft_fact(2)
  WRITE(*,*) 'fft_fact(3) = ', fft_fact(3)

  DO isym = 1,nsym
    WRITE(*,'(1x,A,I2,A,3F18.10)') 'ft(:,', isym, ') = ', ft(:,isym)
  ENDDO 

  DO isym = 1,nsym
    WRITE(*,'(1x,A,I2,A,3I5)') 'ftau(:,', isym, ') = ', ftau(:,isym)
  ENDDO 

  DO isym = 1,nsym
    WRITE(*,*) 'irot = ', isym
    WRITE(*,*) s(1,:,isym)
    WRITE(*,*) s(2,:,isym)
    WRITE(*,*) s(3,:,isym)
  ENDDO 


  ft_ = 0.d0

  WRITE(*,*) 'nsym = ', nsym
  WRITE(*,*)
  DO isym = 1,nsym
    non_symmorphic(isym) = ( ft(1,isym) /= 0.d0 .OR. &
                             ft(2,isym) /= 0.d0 .OR. &
                             ft(3,isym) /= 0.d0 )
    WRITE(*,*) 'ft = ', ft(:,isym)
    WRITE(*,*) 'non_symmorphic = ', non_symmorphic(isym)
    WRITE(*,*)
    IF(non_symmorphic(isym)) THEN
      ft_(:,isym) = at(:,1)*ft(1,isym) + at(:,2)*ft(2,isym) + at(:,3)*ft(3,isym)
    ENDIF 
  ENDDO 

  WRITE(*,*) 'No of non_symmorphic = ', count(non_symmorphic)


  CALL my_sym_rho_init_shells(ngm, g)
  WRITE(*,*) 'Ngs = ', Ngs

  DO isym = 1,nsym
    WRITE(*,'(1x,I5,3F18.10)') isym, ft_(:,isym)
  ENDDO 

  WRITE(*,*) at(1,1), at(1,2), at(1,3)
  WRITE(*,*) at(2,1), at(2,2), at(2,3)
  WRITE(*,*) at(3,1), at(3,2), at(3,3)
  WRITE(*,*) 'alat = ', alat

  WRITE(*,*) 'tpiba = ', tpiba
  WRITE(*,*) g(:,2)*tpiba

  DO isym = 1,nsym
    DO ia = 1,size(irt,2)
      WRITE(*,'(1x,I4,x)',advance='no') irt(isym,ia)
    ENDDO 
    WRITE(*,*)
  ENDDO 

  v(:,1) = (/ -0.1d0, 1.d0, 0.d0 /)
  v(:,2) = (/  0.0d0, 0.1d0, 1.d0 /)

  CALL my_symvector(2, v)
  WRITE(*,*) v(:,1)
  WRITE(*,*) v(:,2)


END SUBROUTINE 
