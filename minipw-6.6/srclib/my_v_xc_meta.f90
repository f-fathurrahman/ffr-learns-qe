!----------------------------------------------------------------------------
SUBROUTINE my_v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v, kedtaur )
!----------------------------------------------------------------------------
  !! Exchange-Correlation potential (meta) Vxc(r) from n(r)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : g, ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE funct,            ONLY : get_meta, dft_is_nonlocc, nlc, is_libxc
  USE xc_mgga,          ONLY : xc_metagcx
  USE scf,              ONLY : scf_type, rhoz_or_updw
  USE mp,               ONLY : mp_sum
  USE mp_bands,         ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
  !! the core charge in real space
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
  !! the core charge in reciprocal space
  REAL(DP), INTENT(INOUT) :: v(dfftp%nnr,nspin)
  !! V_xc potential
  REAL(DP), INTENT(INOUT) :: kedtaur(dfftp%nnr,nspin)
  !! local K energy density                     
  REAL(DP), INTENT(INOUT) :: vtxc
  !! integral V_xc * rho
  REAL(DP), INTENT(INOUT) :: etxc
  !! E_xc energy
  !
  ! ... local variables
  !
  REAL(DP) :: zeta, rh, sgn(2)
  INTEGER  :: k, ipol, is, np
  !
  REAL(DP), ALLOCATABLE :: ex(:), ec(:)
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:), v3x(:,:)
  REAL(DP), ALLOCATABLE :: v1c(:,:), v2c(:,:,:), v3c(:,:)
  !
  REAL(DP) :: fac
       
  REAL(DP), DIMENSION(2) :: grho2, rhoneg
  REAL(DP), DIMENSION(3) :: grhoup, grhodw
  !
  REAL(DP), ALLOCATABLE :: grho(:,:,:), h(:,:,:), dh(:)
  REAL(DP), ALLOCATABLE :: rhoout(:)
  COMPLEX(DP), ALLOCATABLE :: rhogsum(:)
  REAL(DP), PARAMETER :: eps12 = 1.0d-12, zero=0._dp
  !
  CALL start_clock( 'v_xc_meta' )
  !
  etxc = zero
  vtxc = zero
  v(:,:) = zero
  rhoneg(:) = zero
  sgn(1) = 1._dp  ;   sgn(2) = -1._dp
  fac = 1.D0 / DBLE( nspin )
  np = 1
  IF (nspin==2) np=3
  !
  ALLOCATE( grho(3,dfftp%nnr,nspin) )
  ALLOCATE( h(3,dfftp%nnr,nspin) )
  ALLOCATE( rhogsum(ngm) )
  !
  ALLOCATE( ex(dfftp%nnr), ec(dfftp%nnr) )
  ALLOCATE( v1x(dfftp%nnr,nspin), v2x(dfftp%nnr,nspin)   , v3x(dfftp%nnr,nspin) )
  ALLOCATE( v1c(dfftp%nnr,nspin), v2c(np,dfftp%nnr,nspin), v3c(dfftp%nnr,nspin) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  ! ... in LSDA case rhogsum is in (up,down) format
  !
  DO is = 1, nspin
     !
     rhogsum(:) = fac*rhog_core(:) + ( rho%of_g(:,1) + sgn(is)*rho%of_g(:,nspin) )*0.5D0
     !
     CALL fft_gradient_g2r( dfftp, rhogsum, g, grho(1,1,is) )
     !
  ENDDO
  DEALLOCATE(rhogsum)
  !
  !
  IF (nspin == 1) THEN
    !
    CALL xc_metagcx( dfftp%nnr, 1, np, rho%of_r, grho, rho%kin_r/e2, ex, ec, &
                      v1x, v2x, v3x, v1c, v2c, v3c )
    !
    DO k = 1, dfftp%nnr
       !
       v(k,1) = (v1x(k,1)+v1c(k,1)) * e2
       !
       ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
       h(:,k,1) = (v2x(k,1)+v2c(1,k,1)) * grho(:,k,1) * e2 
       !
       kedtaur(k,1) = (v3x(k,1)+v3c(k,1)) * 0.5d0 * e2
       !
       etxc = etxc + (ex(k)+ec(k)) * e2
       vtxc = vtxc + (v1x(k,1)+v1c(k,1)) * e2 * ABS(rho%of_r(k,1))
       !
       IF (rho%of_r(k,1) < zero) rhoneg(1) = rhoneg(1)-rho%of_r(k,1)
       !
    ENDDO
    !
  ELSE
    !
    CALL rhoz_or_updw( rho, 'only_r', '->updw' )
    !
    CALL xc_metagcx( dfftp%nnr, 2, np, rho%of_r, grho, rho%kin_r/e2, ex, ec, &
                     v1x, v2x, v3x, v1c, v2c, v3c )
    !
    ! first term of the gradient correction : D(rho*Exc)/D(rho)
    !
    DO k = 1, dfftp%nnr
       !
       v(k,1) = (v1x(k,1) + v1c(k,1)) * e2
       v(k,2) = (v1x(k,2) + v1c(k,2)) * e2
       !
       ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
       !
       h(:,k,1) = (v2x(k,1) * grho(:,k,1) + v2c(:,k,1)) * e2
       h(:,k,2) = (v2x(k,2) * grho(:,k,2) + v2c(:,k,2)) * e2
       !
       kedtaur(k,1) = (v3x(k,1) + v3c(k,1)) * 0.5d0 * e2
       kedtaur(k,2) = (v3x(k,2) + v3c(k,2)) * 0.5d0 * e2
       !
       etxc = etxc + (ex(k)+ec(k)) * e2
       vtxc = vtxc + (v1x(k,1)+v1c(k,1)) * ABS(rho%of_r(k,1)) * e2
       vtxc = vtxc + (v1x(k,2)+v1c(k,2)) * ABS(rho%of_r(k,2)) * e2
       !
       IF ( rho%of_r(k,1) < 0.d0 ) rhoneg(1) = rhoneg(1) - rho%of_r(k,1)
       IF ( rho%of_r(k,2) < 0.d0 ) rhoneg(2) = rhoneg(2) - rho%of_r(k,2)
       !
    ENDDO
    !
    CALL rhoz_or_updw( rho, 'only_r', '->rhoz' )
    !
  ENDIF
  !
  DEALLOCATE( ex, ec )
  DEALLOCATE( v1x, v2x, v3x )
  DEALLOCATE( v1c, v2c, v3c )
  !
  !
  ALLOCATE( dh( dfftp%nnr ) )    
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  ALLOCATE (rhoout(dfftp%nnr))
  DO is = 1, nspin
     !
     CALL fft_graddot( dfftp, h(1,1,is), g, dh )
     !
     v(:,is) = v(:,is) - dh(:)
     !
     ! ... rhoout is in (up,down) format 
     !
     rhoout(:) = ( rho%of_r(:,1) + sgn(is)*rho%of_r(:,nspin) )*0.5D0
     vtxc = vtxc - SUM( dh(:) * rhoout(:) )
     !
  END DO
  DEALLOCATE(rhoout)
  DEALLOCATE(dh)
  !
  CALL mp_sum( rhoneg, intra_bgrp_comm )
  !
  rhoneg(:) = rhoneg(:) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ((rhoneg(1) > eps8) .OR. (rhoneg(2) > eps8)) THEN
    write (stdout, '(/,5x, "negative rho (up,down): ", 2es10.3)') rhoneg(:)
  ENDIF
  !
  vtxc = omega * vtxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 ) 
  etxc = omega * etxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ( dft_is_nonlocc() ) CALL nlc( rho%of_r, rho_core, nspin, etxc, vtxc, v )
  !
  CALL mp_sum(  vtxc , intra_bgrp_comm )
  CALL mp_sum(  etxc , intra_bgrp_comm )
  !
  DEALLOCATE(grho)
  DEALLOCATE(h)
  !
  CALL stop_clock( 'v_xc_meta' )
  !
  RETURN
  !
END SUBROUTINE my_v_xc_meta