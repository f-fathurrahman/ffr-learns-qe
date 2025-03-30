!---------------------------------------------------------------
SUBROUTINE ld1x_debug_solve_sch()
!---------------------------------------------------------------
  ! modeled after ld1x_all_electron
  USE kinds, ONLY : DP
  USE radial_grids, ONLY: ndmx
  USE constants, ONLY: e2  
  USE ld1inc, ONLY: grid, enne, vpot, vxt, enl, &
                 & nspin, nn, ll, oc, nwf, &
                 & zed, zval, v0, el, &
                 & isw, core_state, rel, frozen_core, psi
  IMPLICIT NONE
  !
  ! automatic arrays
  REAL(DP) :: vnew(ndmx,2)
  ! originally arguments
  INTEGER :: ic
  LOGICAL :: ild
  !
  integer :: iwf, ispin
  integer :: nstop, nerr
  real(dp) :: ze2
  INTEGER, PARAMETER :: maxter=200
  REAL(DP), PARAMETER :: thresh=1.0e-10_dp
  integer :: nin

  ! what? convert to Ry and change the sign
  ze2 = -zed * e2
  ! This is used in ascheq

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
  ! XXX: isic stuffs are removed
  ! XXX: metaGGA version is removed


  vnew = vpot


  ! loop over all states (number of wavefunctions)
  DO iwf = 1,nwf
    !
    ! Solve one-particle equation (Schroedinger or Dirac)
    !
    ! only solve when occupation is not zero
    IF( oc(iwf) >= 0.0_dp ) THEN
      ! 
      IF( ic==1 .or. .not. frozen_core .or. .not. core_state(iwf) ) THEN
        !
        ispin = isw(iwf) ! spin index
        !    
        IF( rel == 0 ) THEN
          !
          ! nonrelativistic calculation
          !
          CALL my_ascheq( nn(iwf), ll(iwf), enl(iwf), grid%mesh, grid, vnew(:,ispin), & ! potential
                      &  ze2, thresh, psi(:,:,iwf), nstop )
          !
        ELSEIF( rel == 1 ) THEN
          !
          ! scalar relativistic
          !
          CALL my_lschps( 1, zed, thresh, grid, nin, nn(iwf), ll(iwf), &
                        &  enl(iwf), vnew(:,ispin), psi(:,:,iwf), nstop)
          ! mode = 1 find energy and wavefunction of bound states,
          !          scalar-relativistic (all-electron)
          ! XXX what's this?
          IF( nstop > 0 .and. oc(iwf) < 1.e-10_DP) then
            nstop = 0
          endif
          !
        ELSE
          !
          CALL errore('ld1x_debug_v01', 'relativistic not programmed', 1)
          !  
        ENDIF
          ! write(6,*) nn(n),ll(n),enl(n)
          ! if (nstop /= 0) write(6,'(4i6)') iter,nn(n),ll(n),nstop
          nerr = nerr + nstop
      ENDIF
        !  
    ELSE
      !
      ! Case oc(n) is negative, zero out energies and psi
      !
      enl(iwf) = 0.d0
      psi(:,:,iwf) = 0.d0
    ENDIF
    
  ENDDO ! loop over Nwf

  WRITE(*,*)
  WRITE(*,*) 'Energy eigenvalues (in Ha): '
  DO iwf = 1,nwf
    WRITE(*,'(1x,A,F18.10)') 'enl = ', enl(iwf)/2
  ENDDO


  WRITE(*,*)
  WRITE(*,*) '***** EXIT ld1x_debug_solve_sch *****'
  WRITE(*,*)

  RETURN
  !
END SUBROUTINE 

