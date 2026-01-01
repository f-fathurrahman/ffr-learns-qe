include "prepare_all.f90"
include "my_symm_base.f90"


SUBROUTINE test_my_symm_base()
  
  USE ions_base,          ONLY : nat, tau, ityp
  USE noncollin_module,   ONLY : noncolin, m_loc
  USE spin_orb, ONLY: domag
  USE control_flags, ONLY: noinv
  USE extfield, ONLY: gate

  use my_symm_base, ONLY: nsym, nrot, irt, ft, s, sname, fft_fact, t_rev

  IMPLICIT NONE 
  LOGICAL :: magnetic_sym, time_reversal
  INTEGER :: i


  CALL set_sym_bl()
  
  magnetic_sym = noncolin .AND. domag 
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym

  write(*,*) 'gate = ', gate
  CALL find_sym( nat, tau, ityp, magnetic_sym, m_loc, gate )

  DO i = 1,nsym
    WRITE(*,*)
    WRITE(*,*) 'i = ', i
    WRITE(*,*) 'sname = ', sname(i)
    WRITE(*,'(1x,A,3F18.10)') 'ft = ', ft(:,i)
  ENDDO
  WRITE(*,*) 'fft_fact = ', fft_fact

  !WRITE(*,*) 't_rev = ', t_rev ! this is only needed for noncollinear magnetism

  WRITE(*,*) 'nsym = ', nsym
  WRITE(*,*) 'nrot = ', nrot

END SUBROUTINE 


PROGRAM main

  IMPLICIT NONE 
  
  CALL prepare_all()

  CALL test_my_symm_base()

END PROGRAM 



