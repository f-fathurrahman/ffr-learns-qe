! Originally inner function in electrons_scf()

!-----------------------------------------------------------------------
FUNCTION my_delta_e() result(res)
  !-----------------------------------------------------------------------
  ! This function computes delta_e, where:
  !
  ! ... delta_e =  - \int rho%of_r(r)  v%of_r(r)
  !                - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
  !                - \sum rho%ns       v%ns       [for LDA+U]
  !                - \sum becsum       D1_Hxc     [for PAW]
  !
  use kinds, only : dp
  USE funct,  ONLY : dft_is_meta
  use scf, only : rho, v
  use fft_base, only : dfftp
  use lsda_mod, only : nspin
  use cell_base, only : omega
  use mp_bands, only : intra_bgrp_comm
  use paw_variables, only : okpaw, ddd_paw
  USE mp, only : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP) :: res
  INTEGER  :: ir

  res = 0._DP
  IF ( nspin==2 ) THEN
    DO ir = 1,dfftp%nnr
      res = res - ( rho%of_r(ir,1) + rho%of_r(ir,2) ) * v%of_r(ir,1) &  ! up
                - ( rho%of_r(ir,1) - rho%of_r(ir,2) ) * v%of_r(ir,2)    ! dw
    ENDDO 
    res = 0.5_DP*res
  ELSE
    res = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
  ENDIF


  IF( dft_is_meta() ) res = res - SUM( rho%kin_r(:,:)*v%kin_r(:,:) )

  res = omega * res / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )

  CALL mp_sum( res, intra_bgrp_comm )
  
  IF(okpaw) res = res - SUM( ddd_paw(:,:,:)*rho%bec(:,:,:) )
  
  RETURN

END FUNCTION
