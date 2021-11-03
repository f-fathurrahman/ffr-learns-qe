INCLUDE 'prepare_all.f90'

PROGRAM main
  !
  USE kinds, ONLY: DP
  USE fft_base, ONLY : dfftp
  !
  IMPLICIT NONE 

  REAL(DP), ALLOCATABLE :: rhoa(:,:)
  INTEGER :: nspina

  CALL prepare_all()

  nspina = 1
  ALLOCATE( rhoa(dfftp%nnr,nspina) )

  CALL my_atomic_rho(rhoa, nspina)

  DEALLOCATE( rhoa )

END PROGRAM 



!-----------------------------------------------------------------------
SUBROUTINE my_atomic_rho( rhoa, nspina )
  !-----------------------------------------------------------------------
  !! Same as \(\texttt{atomic_rho_g}\), with real-space output charge
  !! \(\text{rhoa}(:,\text{nspina})\).
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  USE lsda_mod,             ONLY : lsda
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_g2r
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nspina
  !! number of spin components to be calculated. It may differ from
  !! nspin because in some cases the total charge only is needed, 
  !! even in a LSDA calculation.
  REAL(DP), INTENT(OUT) :: rhoa(dfftp%nnr,nspina)
  !! contains R-space components of the superposition of atomic charges.
  !
  ! ... local variables
  !
  REAL(DP) :: rhoneg
  COMPLEX(DP), allocatable :: rhocg(:,:)
  INTEGER :: ir, is
  !
  ! allocate work space 
  !
  ALLOCATE (rhocg(dfftp%ngm, nspina))
  !
  CALL atomic_rho_g(rhocg, nspina)
  !
  ! bring to real space
  !
  rhoa(:,:) = 0.d0
  CALL rho_g2r( dfftp, rhocg, rhoa )


  write(*,*)
  write(*,*) 'my_atomic_rhoa: integ rhoa = ', sum(rhoa)*omega/dfftp%nnr
  write(*,*)

  DEALLOCATE(rhocg)

  !
  ! check on negative charge
  !
  ! it is useless to set negative terms to zero in real space: 
  ! negative charge will re-appear when Fourier-transformed back and forth
  !
  DO is = 1, nspina
    rhoneg = 0.0_dp
    DO ir = 1, dfftp%nnr
      rhoneg = rhoneg + MIN( 0.0_dp,  DBLE(rhoa(ir,is)) )
    ENDDO
    rhoneg = omega * rhoneg / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
    !
    CALL mp_sum(rhoneg, intra_bgrp_comm)
    !
    IF( (is == 1) .OR. lsda ) THEN
      !
      IF( (rhoneg < -1.0d-4) ) THEN
        IF( lsda ) THEN 
          WRITE( stdout,'(5x,"Check: negative starting charge=", &
                &"(component",i1,"):",f12.6)') is, rhoneg
        ELSE
          WRITE( stdout,'(5x,"Check: negative starting charge=",f12.6)') rhoneg
        ENDIF
      ENDIF
    ENDIF
  ENDDO

END SUBROUTINE

