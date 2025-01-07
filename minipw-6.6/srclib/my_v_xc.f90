!----------------------------------------------------------------------------
SUBROUTINE my_v_xc( rho, rho_core, rhog_core, etxc, vtxc, v )
!----------------------------------------------------------------------------
  !! Exchange-Correlation potential Vxc(r) from n(r)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2, eps8
  USE io_global,        ONLY : stdout
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE lsda_mod,         ONLY : nspin
  USE cell_base,        ONLY : omega
  USE spin_orb,         ONLY : domag
  USE funct,            ONLY : nlc, dft_is_nonlocc
  USE xc_lda_lsda,      ONLY : xc
  USE scf,              ONLY : scf_type
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE (scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
  !! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
  !! the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: v(dfftp%nnr,nspin)
  !! V_xc potential
  REAL(DP), INTENT(OUT) :: vtxc
  !! integral V_xc * rho
  REAL(DP), INTENT(OUT) :: etxc
  !! E_xc energy
  !
  ! ... local variables
  !
  REAL(DP) :: rhoneg(2), vs
  !
  !REAL(DP), ALLOCATABLE :: arhox(:), amag(:), zeta(:)
  REAL(DP) :: arho, amag
  REAL(DP) :: rhoup2, rhodw2
  REAL(DP), ALLOCATABLE :: ex(:), ec(:)
  REAL(DP), ALLOCATABLE :: vx(:,:), vc(:,:)
  ! In order:
    ! the absolute value of the total charge
    ! the absolute value of the magnetization
    ! zeta = amag / arhox
    ! local exchange energy
    ! local correlation energy
    ! local exchange potential
    ! local correlation potential
  INTEGER :: ir, ipol
    ! counter on mesh points
    ! counter on nspin
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-10, &
                         vanishing_mag    = 1.D-20

  ALLOCATE( ex(dfftp%nnr) )
  ALLOCATE( ec(dfftp%nnr) )
  ALLOCATE( vx(dfftp%nnr,nspin) )
  ALLOCATE( vc(dfftp%nnr,nspin) )
  !
  etxc   = 0.D0
  vtxc   = 0.D0
  v(:,:) = 0.D0
  rhoneg = 0.D0
  !
  !
  rho%of_r(:,1) = rho%of_r(:,1) + rho_core(:) ! XXX add to up only?
  
  write(*,*) 'my_v_xc: sum rho_core: ', sum(rho_core)

  !write(*,*) 'some rho_core: '
  !write(*,*) rho_core(1:3)
  !write(*,*) rho_core(4:6)
  !write(*,*) 'integ rho%of_r(:,1) = ', sum(rho%of_r(:,1))*omega/dfftp%nnr
  !
  IF( nspin == 1 .OR. ( nspin == 4 .AND. .NOT. domag ) ) THEN
    ! ... spin-unpolarized case
    
    write(*,*) 'my_v_xc: sum rho%of_r = ', sum(rho%of_r(:,1))

    write(*,*) 'my_v_xc: before sum abs Vxc (in Ha) = ', sum(abs(v(:,1)))*0.5d0

    CALL xc( dfftp%nnr, 1, 1, rho%of_r, ex, ec, vx, vc )

    write(*,*) 'shape vx = ', shape(vx)
    write(*,*) 'sum vx + vc = ', sum(vx + vc)

    !
    DO ir = 1, dfftp%nnr
      v(ir,1) = e2*( vx(ir,1) + vc(ir,1) )
      etxc = etxc + e2*( ex(ir) + ec(ir) )*rho%of_r(ir,1)
      rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
      vtxc = vtxc + v(ir,1)*rho%of_r(ir,1)
      IF (rho%of_r(ir,1) < 0.D0) rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
    ENDDO

    write(*,*) 'my_v_xc: sum epsxc   (in Ha) = ', sum(abs(ex + ec))*0.5d0
    write(*,*) 'my_v_xc: sum abs Vxc (in Ha) = ', sum(abs(v(:,1)))*0.5d0
    !
  ELSEIF ( nspin == 2 ) THEN
     ! ... spin-polarized case
     !
     CALL xc( dfftp%nnr, 2, 2, rho%of_r, ex, ec, vx, vc )
     ! note that rho%of_r contains rho_core contribution?
     !
     DO ir = 1, dfftp%nnr   !OMP ?
        v(ir,:) = e2*( vx(ir,:) + vc(ir,:) ) ! convert to Ry?
        etxc = etxc + e2*( (ex(ir) + ec(ir))*rho%of_r(ir,1) )
        !
        rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)  ! recover original rhoe
        !
        vtxc = vtxc + ( ( v(ir,1) + v(ir,2) )*rho%of_r(ir,1) + &
                        ( v(ir,1) - v(ir,2) )*rho%of_r(ir,2) )
        !
        rhoup2 = rho%of_r(ir,1) + rho%of_r(ir,2)
        rhodw2 = rho%of_r(ir,1) - rho%of_r(ir,2)
        IF (rhoup2 < 0.d0) rhoneg(1) = rhoneg(1) + rhoup2
        IF (rhodw2 < 0.d0) rhoneg(2) = rhoneg(2) + rhodw2
     ENDDO
     !
     vtxc   = 0.5d0 * vtxc
     rhoneg = 0.5d0 * rhoneg
     !
     !
  ELSE IF ( nspin == 4 ) THEN
     ! ... noncolinear case
     !
     CALL xc( dfftp%nnr, 4, 2, rho%of_r, ex, ec, vx, vc )
     !
     DO ir = 1, dfftp%nnr  !OMP ?
        arho = ABS( rho%of_r(ir,1) )
        IF ( arho < vanishing_charge ) CYCLE
        vs = 0.5D0*( vx(ir,1) + vc(ir,1) - vx(ir,2) - vc(ir,2) )
        v(ir,1) = e2*( 0.5D0*( vx(ir,1) + vc(ir,1) + vx(ir,2) + vc(ir,2) ) )
        !
        amag = SQRT( SUM( rho%of_r(ir,2:4)**2 ) )
        IF ( amag > vanishing_mag ) THEN
           v(ir,2:4) = e2 * vs * rho%of_r(ir,2:4) / amag
           vtxc = vtxc + SUM( v(ir,2:4) * rho%of_r(ir,2:4) )
        ENDIF
        etxc = etxc + e2*( ex(ir) + ec(ir) ) * arho
        !
        rho%of_r(ir,1) = rho%of_r(ir,1) - rho_core(ir)
        IF ( rho%of_r(ir,1) < 0.D0 )  rhoneg(1) = rhoneg(1) - rho%of_r(ir,1)
        IF ( amag / arho > 1.D0 )  rhoneg(2) = rhoneg(2) + 1.D0/omega
        vtxc = vtxc + v(ir,1) * rho%of_r(ir,1)
     ENDDO
     !
     !
  ENDIF
  !
  DEALLOCATE( ex, ec )
  DEALLOCATE( vx, vc )
  !
  CALL mp_sum(  rhoneg , intra_bgrp_comm )
  !
  rhoneg(:) = rhoneg(:) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  IF ( rhoneg(1) > eps8 .OR. rhoneg(2) > eps8 ) &
     WRITE( stdout,'(/,5X,"negative rho (up, down): ",2ES10.3)') rhoneg
  !
  ! ... energy terms, local-density contribution
  !
  vtxc = omega * vtxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  etxc = omega * etxc / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )

  !
  ! ... add gradient corrections (if any)
  !
  CALL gradcorr( rho%of_r, rho%of_g, rho_core, rhog_core, etxc, vtxc, v )
  !
  ! ... add non local corrections (if any)
  !
  IF ( dft_is_nonlocc() ) CALL nlc( rho%of_r, rho_core, nspin, etxc, vtxc, v )
  !
  CALL mp_sum(  vtxc , intra_bgrp_comm )
  CALL mp_sum(  etxc , intra_bgrp_comm )

  RETURN

END SUBROUTINE my_v_xc