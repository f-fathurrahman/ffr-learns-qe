! Originally inner subroutine in electrons_scf()

!-----------------------------------------------------------------------
SUBROUTINE my_compute_magnetization()
!-----------------------------------------------------------------------
  use kinds, only : dp
  use scf, only : rho
  use fft_base, only : dfftp
  use lsda_mod, only : lsda, magtot, absmag
  use noncollin_module, only : noncolin, magtot_nc, bfield
  use mp_bands, only : intra_bgrp_comm
  use klist, only : two_fermi_energies, lgauss
  use cell_base, only: omega
  use ener, only: ef_up, ef_dw
  IMPLICIT NONE
  ! local
  INTEGER :: ir, i
  REAL(DP) :: mag !! local magnetization
  !
  IF( lsda ) THEN
    magtot = 0.D0
    absmag = 0.D0
    DO ir = 1, dfftp%nnr
      mag = rho%of_r(ir,2)
      magtot = magtot + mag
      absmag = absmag + ABS( mag )
    ENDDO

    magtot = magtot * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
    absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
    
    CALL mp_sum( magtot, intra_bgrp_comm )
    CALL mp_sum( absmag, intra_bgrp_comm )

    IF( two_fermi_energies .and. lgauss ) bfield(3) = 0.5D0*(ef_up - ef_dw)

  ELSEIF( noncolin ) THEN

    magtot_nc = 0.D0
    absmag    = 0.D0
    DO ir = 1,dfftp%nnr
      mag = SQRT( rho%of_r(ir,2)**2 + &
                  rho%of_r(ir,3)**2 + &
                  rho%of_r(ir,4)**2 )
      DO i = 1, 3
        magtot_nc(i) = magtot_nc(i) + rho%of_r(ir,i+1)
      ENDDO
      absmag = absmag + ABS( mag )
    ENDDO
    CALL mp_sum( magtot_nc, intra_bgrp_comm )
    CALL mp_sum( absmag, intra_bgrp_comm )

    DO i = 1, 3
      magtot_nc(i) = magtot_nc(i) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
    ENDDO
    absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  ENDIF
  !
  RETURN
END SUBROUTINE

