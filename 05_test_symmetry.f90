include "my_symme.f90"


PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL test_symmetry()
END PROGRAM 


SUBROUTINE test_symmetry()
  
  USE symm_base, ONLY : nrot, fft_fact, ft, ftau, s, nsym
  USE gvect, ONLY : ngm, g
  USE my_symme, ONLY: Ngs, my_sym_rho_init_shells
  USE cell_base, ONLY: at, alat, tpiba

  IMPLICIT NONE 
  INTEGER :: isym
  LOGICAL :: non_symmorphic(48)
  REAL(8) :: ft_(3,48), sg(3), arg

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
END SUBROUTINE 
