! similar to ld1x_scf, but with many debugging prints
!
!---------------------------------------------------------------
SUBROUTINE ld1x_my_scf(ic)
  !---------------------------------------------------------------
  !
  ! this routine performs the atomic self-consistent procedure
  ! self-interaction-correction allowed
  !
  USE kinds, ONLY : dp
  USE funct, ONLY : dft_is_meta
  USE radial_grids, ONLY : ndmx
  USE constants, ONLY: e2
  USE ld1inc, ONLY : grid, zed, psi, isic, vpot, vh, vxt, rho, iter, &
                     lsd, rel, latt, enne, beta, nspin, tr2, eps0, &
                     nwf, nn, ll, jj, enl, oc, isw, core_state, frozen_core, &
                     tau, vtau, vsic, vsicnew, vhn1, egc, relpert, noscf
  IMPLICIT NONE
 
  INTEGER, INTENT(in) :: ic

  LOGICAL:: meta, conv
  INTEGER:: nerr, nstop, n, i, is, id, nin
  REAL(DP) ::  vnew(ndmx,2), vtaunew(ndmx), rhoc1(ndmx), ze2
  INTEGER, PARAMETER :: maxter=200
  REAL(DP), PARAMETER :: thresh=1.0e-10_dp
  INTEGER :: ii
  REAL(DP) :: integRho

  meta = dft_is_meta()
  ze2 = -zed * e2
  rhoc1 = 0.0_dp
  
  IF(.not. frozen_core .or. ic == 1) psi = 0.0_dp
  
  DO iter = 1,maxter

    WRITE(*,*)
    WRITE(*,*) '------------------------------------'
    WRITE(*,*) 'Begin iterSCF = ', iter
    WRITE(*,*) '------------------------------------'

    !
    nerr = 0
    vnew = vpot
    vtaunew = vtau
    !
    ! loop over all states (number of wavefunctions)
    DO n = 1,nwf
      !
      ! Solve one-particle equation (Schroedinger or Dirac)
      !
      ! only solve when occupation is not zero
      IF( oc(n) >= 0.0_dp ) THEN
        
        IF( ic==1 .or. .not. frozen_core .or. .not. core_state(n) ) THEN
          !
          is = isw(n) ! spin index
          !
          IF( isic /= 0 .and. iter > 1 ) THEN
            vnew(:,is) = vpot(:,is) - vsic(:,n)
          ENDIF 
          !    
          IF( rel == 0 ) THEN
            !
            ! nonrelativistic calculation
            !
            IF( meta ) THEN
              !
              ! Meta-GGA version of lschps
              !
              CALL lschps_meta( 2, zed, thresh, grid, nin, nn(n), ll(n),&
                         enl(n), vnew(:,is), vtaunew, psi(:,:,n), nstop )
            ELSE
              !
              ! Non meta-GGA
              !
              ! This is the "normal" case
              CALL my_ascheq( nn(n), ll(n), enl(n), grid%mesh, grid, vnew(:,is), & ! potential
                    &  ze2, thresh, psi(:,:,n), nstop )
               
            ENDIF ! meta
            !
          ELSEIF( rel == 1 ) THEN
            !
            ! relativistic scalar calculation
            !
            IF( meta ) THEN
              CALL lschps_meta( 1, zed, thresh, grid, nin, nn(n), ll(n), &
                             &  enl(n), vnew(:,is), vtaunew, psi(:,:,n), nstop)
            ELSE
              CALL my_lschps( 1, zed, thresh, grid, nin, nn(n), ll(n), &
                        &  enl(n), vnew(:,is), psi(:,:,n), nstop)
            ENDIF
            IF( nstop > 0 .and. oc(n) < 1.e-10_DP) nstop=0
            !
          ELSEIF( rel == 2 ) THEN
            !
            ! Dirac equation
            !
            CALL dirsol( ndmx, grid%mesh, nn(n), ll(n), jj(n), iter, enl(n), &
                    &    thresh, grid, psi(1,1,n), vnew(1,is), nstop )
          ELSE
            CALL errore('scf', 'relativistic not programmed', 1)
          
          ENDIF
          !      write(6,*) nn(n),ll(n),enl(n)
          ! if (nstop /= 0) write(6,'(4i6)') iter,nn(n),ll(n),nstop
          nerr = nerr + nstop
        ENDIF
        !  
      ELSE
        !
        ! Case oc(n) is negative, zero out energies and psi
        !
        enl(n) = 0.0_dp
        psi(:,:,n) = 0.0_dp
      
      ENDIF
    
    ENDDO ! loop over Nwf

    WRITE(*,*)
    WRITE(*,*) 'Energy eigenvalues (in Ha): '
    DO ii = 1,nwf
      WRITE(*,'(1x,A,F18.10)') 'enl = ', enl(ii)/2
      !write(*,*) 'psi1 = ', psi(:,:,ii)
    ENDDO

    !
    ! calculate charge density (spherical approximation)
    !
    rho = 0.0_dp
    IF( noscf ) GOTO 500
    DO n = 1,nwf
      rho(1:grid%mesh,isw(n)) = rho(1:grid%mesh,isw(n)) + &
        & oc(n)*( psi(1:grid%mesh,1,n)**2 + psi(1:grid%mesh,2,n)**2 )
    ENDDO
    CALL simpson(grid%mesh, rho(1:grid%mesh,1), grid%rab, integRho) 
    ! XXX: Only integrate up spin?
    WRITE(*,*)
    WRITE(*,*) 'SCF: integrated rho = ', integRho
    !
    ! Calculate kinetic energy density (spherical approximation) if needed
    !
    IF( meta ) then
      CALL kin_e_density(ndmx, grid%mesh, nwf, ll, oc, psi, grid%r, grid%r2, grid%dx, tau)
    endif
    
    !
    ! calculate new potential
    !
    CALL my_new_potential( ndmx, grid%mesh, grid, zed, vxt, &
       &                lsd, .false., latt, enne, rhoc1, rho, vh, vnew, 1 )
    
    !
    ! calculate SIC correction potential (if present)
    !
    IF (isic /= 0) THEN
      DO n=1,nwf
        IF (oc(n) >= 0.0_dp) THEN
          is=isw(n)
          CALL sic_correction(n,vhn1,vsicnew,egc)
          !
          ! use simple mixing for SIC correction
          !
          vsic(:,n) = (1.0_dp - beta)*vsic(:,n)+beta*vsicnew(:)
        ENDIF
      ENDDO
    ENDIF

    !
    ! mix old and new potential
    !
    id = 3
    IF( isic /= 0 .and. relpert )  id=1
    !
    CALL vpack(grid%mesh, ndmx, nspin, vnew, vpot, 1)
    CALL dmixp(grid%mesh*nspin, vnew, vpot, beta, tr2, iter, id, eps0, conv, maxter)
    CALL vpack(grid%mesh, ndmx, nspin, vnew, vpot, -1)
    
    WRITE(*,'(1x,A,I5,ES18.10)') 'SCF iter ', iter, eps0
    
    !
    ! mix old and new metaGGA potential - use simple mixing
    !
    IF( meta ) vtau(:) = (1.0_dp-beta)*vtaunew(:)+beta*vtau(:)



500 IF( noscf ) THEN
      conv = .true.
      eps0 = 0.0_DP
    ENDIF

    IF( conv ) THEN
      IF (nerr /= 0) CALL infomsg ('scf','warning: at least one error in KS equations')
      EXIT ! exit cycle
    ENDIF
  ENDDO
  
  IF( .not. conv ) CALL infomsg('scf','warning: convergence not achieved')

END SUBROUTINE ld1x_my_scf

