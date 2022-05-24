! originally an inner function in electrons_scf

!-----------------------------------------------------------------------
FUNCTION my_delta_escf(rhoin, rho) result(res)
!-----------------------------------------------------------------------
  ! This function calculates the difference between the Hartree and XC energy
  ! at first order in the charge density difference delta_rho(r).
  !
  ! ... delta_escf = - \int \delta rho%of_r(r)  v%of_r(r)
  !                  - \int \delta rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
  !                  - \sum \delta rho%ns       v%ns       [for LDA+U]
  !                  - \sum \delta becsum       D1         [for PAW] 
  !
  USE funct,  ONLY : dft_is_meta
  use kinds, only : dp
  USE funct,  ONLY : dft_is_meta
  use scf, only : v, scf_type
  use fft_base, only : dfftp
  use lsda_mod, only : nspin
  use cell_base, only : omega
  use mp_bands, only : intra_bgrp_comm
  use paw_variables, only : okpaw, ddd_paw
  USE mp, only : mp_sum
  !
  IMPLICIT NONE
  type(scf_type) :: rhoin, rho
  !
  REAL(DP) :: res
  real(dp) :: rho_dif(2)
  INTEGER  :: ir
  !
  res = 0._dp
  IF( nspin==2 ) THEN
    DO ir=1, dfftp%nnr
       rho_dif = rhoin%of_r(ir,:) - rho%of_r(ir,:)
       res = res - ( rho_dif(1) + rho_dif(2) ) * v%of_r(ir,1) &  !up
                 - ( rho_dif(1) - rho_dif(2) ) * v%of_r(ir,2)    !dw
    ENDDO
    res = 0.5_dp*res
    !
  ELSE
    ! non-spin polarized case
    res = -SUM( ( rhoin%of_r(:,:)-rho%of_r(:,:) )*v%of_r(:,:) )
  ENDIF
  

  IF( dft_is_meta() ) then
    res = res - SUM( (rhoin%kin_r(:,:)-rho%kin_r(:,:) )*v%kin_r(:,:))
  endif
  res = omega * res / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  CALL mp_sum( res, intra_bgrp_comm )


  IF( okpaw ) then
    res = res - SUM(ddd_paw(:,:,:)*(rhoin%bec(:,:,:)-rho%bec(:,:,:)))
  endif

  RETURN
END FUNCTION