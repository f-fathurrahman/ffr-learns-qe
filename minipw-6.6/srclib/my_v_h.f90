!----------------------------------------------------------------------------
SUBROUTINE my_v_h( rhog, ehart, charge, v )
  !----------------------------------------------------------------------------
  !! Hartree potential VH(r) from n(G)
  !
  USE constants,         ONLY : fpi, e2
  USE kinds,             ONLY : DP
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : invfft
  USE gvect,             ONLY : ngm, gg, gstart
  USE lsda_mod,          ONLY : nspin
  USE cell_base,         ONLY : omega, tpiba2
  USE control_flags,     ONLY : gamma_only
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_sum
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  USE esm,               ONLY : do_comp_esm, esm_hartree, esm_bc
  USE Coul_cut_2D,       ONLY : do_cutoff_2D, cutoff_2D, cutoff_hartree  
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: rhog(ngm)
  !! the charge density in reciprocal space
  REAL(DP), INTENT(INOUT) :: v(dfftp%nnr,nspin)
  !! Hartree potential
  REAL(DP), INTENT(OUT) :: ehart
  !! Hartree energy
  REAL(DP), INTENT(OUT) :: charge
  !
  !  ... local variables
  !
  REAL(DP)              :: fac
  REAL(DP), ALLOCATABLE :: aux1(:,:)
  REAL(DP)              :: rgtot_re, rgtot_im, eh_corr
  INTEGER               :: is, ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), rgtot(:), vaux(:)
  INTEGER               :: nt
  !
  ALLOCATE( aux( dfftp%nnr ), aux1( 2, ngm ) )
  charge = 0.D0
  IF( gstart == 2 ) THEN
    charge = omega*REAL( rhog(1) )
  ENDIF
  CALL mp_sum( charge , intra_bgrp_comm )
  
  ! calculate hartree potential in G-space (NB: V(G=0)=0 )
  IF( do_comp_esm .and. ( esm_bc .ne. 'pbc' ) ) THEN
    ! calculate modified Hartree potential for ESM
    CALL esm_hartree (rhog, ehart, aux)
  ELSE
    ehart     = 0.D0
    aux1(:,:) = 0.D0
    !
    IF (do_cutoff_2D) THEN  !TS
      CALL cutoff_hartree(rhog(:), aux1, ehart)
    ELSE
      DO ig = gstart, ngm
        fac = 1.D0 / gg(ig) 
        rgtot_re = REAL( rhog(ig) )
        rgtot_im = AIMAG( rhog(ig) )
        ehart = ehart + ( rgtot_re**2 + rgtot_im**2 ) * fac
        aux1(1,ig) = rgtot_re * fac
        aux1(2,ig) = rgtot_im * fac
      ENDDO
    ENDIF
    
    fac = e2 * fpi / tpiba2
    ehart = ehart * fac
    aux1 = aux1 * fac

    IF( gamma_only ) THEN
      ehart = ehart * omega
    ELSE 
      ehart = ehart * 0.5D0 * omega
    ENDIF

    write(*,*)
    write(*,*) 'my_v_h: ehart = ', ehart*0.5d0
    write(*,*)

    if( do_comp_mt ) then
      ALLOCATE( vaux( ngm ), rgtot(ngm) )
      rgtot(:) = rhog(:)
      CALL wg_corr_h( omega, ngm, rgtot, vaux, eh_corr )
      aux1(1,1:ngm) = aux1(1,1:ngm) + REAL( vaux(1:ngm))
      aux1(2,1:ngm) = aux1(2,1:ngm) + AIMAG(vaux(1:ngm))
      ehart = ehart + eh_corr
      DEALLOCATE( rgtot, vaux )
    endif

    CALL mp_sum( ehart , intra_bgrp_comm )

    aux(:) = 0.D0
    aux(dfftp%nl(1:ngm)) = CMPLX( aux1(1,1:ngm), aux1(2,1:ngm), KIND=dp )

    IF( gamma_only ) THEN
      aux(dfftp%nlm(1:ngm)) = CMPLX( aux1(1,1:ngm), -aux1(2,1:ngm), KIND=dp )
    ENDIF

  ENDIF

  ! transform hartree potential to real space
  CALL invfft('Rho', aux, dfftp)

  write(*,*)
  write(*,*) 'my_v_h: sum abs aux after invfft (in Ha): ', sum(abs(aux))*0.5d0
  write(*,*)

  ! add hartree potential to the xc potential
  IF( nspin == 4 ) THEN
    v(:,1) = v(:,1) + DBLE(aux(:))
  ELSE
    DO is = 1, nspin
      v(:,is) = v(:,is) + DBLE(aux(:))
    ENDDO
  ENDIF

  DEALLOCATE( aux, aux1 )
  RETURN

END SUBROUTINE my_v_h

