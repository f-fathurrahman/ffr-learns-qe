!---------------------------------------------------------------
SUBROUTINE ld1x_debug_v01()
!---------------------------------------------------------------
  ! modeled after ld1x_all_electron


  !
  !  this routine is a driver to an all-electron calculation
  !  with the parameters given in input
  !
  !
  USE kinds, ONLY : DP
  USE radial_grids, ONLY: ndmx
  USE ld1inc, ONLY: isic, grid, rho, enne, vpot, vxt, enl, &
                 & deld, encl, etot, ecxc, evxt, ehrt, ekin, &
                 & vh, nspin, vdw, nn, ll, oc, nwf, &
                 & zed, zval, vxc, exc, excgga, v0, verbosity, &
                 & relpert, evel, edar, eso, vsic, vsicnew, vhn1, egc, el, &
                 & isw
  IMPLICIT NONE

  ! originally arguments
  INTEGER :: ic
  LOGICAL :: ild
  !
  integer :: iwf

  ic = 1
  ild = .false.

  WRITE(*,*)
  WRITE(*,*) '***** ENTER ld1x_debug_v01 *****'
  WRITE(*,*)

  ! print out some variables
  WRITE(*,*) 'ndmx = ', ndmx
  WRITE(*,*) 'zed = ', zed, ''
  WRITE(*,*) 'zval = ', zval
  WRITE(*,*) 'nwfx = ', size(oc)
  WRITE(*,*) 'nwf = ', nwf
  write(*,*) 'Configurations: iwf, label (el), isw (spin), nn, ll, oc'
  do iwf = 1,nwf
    write(*,'(1x,I4,A3,3I4,F8.3)') iwf, el(iwf), isw(iwf), nn(iwf), ll(iwf), oc(iwf)
  enddo

  CALL starting_potential( ndmx, grid%mesh, zval, zed, nwf, oc, nn, ll,&
                           grid%r, enl, v0, vxt, vpot, enne, nspin )

  !
  ! isic stuffs are removed
  !

  ! loop over all states (number of wavefunctions)
  DO n = 1,nwf
    !
    ! Solve one-particle equation (Schroedinger or Dirac)
    !
    ! only solve when occupation is not zero
    IF( oc(n) >= 0.0_dp ) THEN
    
        IF( ic==1 .or. .not. frozen_core .or. .not. core_state(iwf) ) THEN
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
  WRITE(*,*) '***** EXIT ld1x_debug_v01 *****'
  WRITE(*,*)

  RETURN
  !
END SUBROUTINE 

