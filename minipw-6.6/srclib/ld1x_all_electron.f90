!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
SUBROUTINE all_electron(ild,ic)
  !---------------------------------------------------------------
  !
  !  this routine is a driver to an all-electron calculation
  !  with the parameters given in input
  !
  !
  USE kinds, ONLY : DP
  USE radial_grids, ONLY: ndmx
  USE ld1inc, ONLY: isic, grid, zeta, rho, enne, vpot, vxt, enl, &
                     deld, encl, etot, ecxc, evxt, ehrt, epseu, ekin, &
                     vnl, vh, lsd, nspin, nlcc, vdw, nn, ll, oc, nwf, &
                     zed, zval, vxc, exc, excgga, v0, verbosity, &
                     relpert, evel, edar, eso, vsic, vsicnew, vhn1, egc, el
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ic   ! counter on configurations
  LOGICAL :: ild    ! if true compute log der
  INTEGER :: i

  WRITE(*,*)
  WRITE(*,*) 'Enter all_electron'
  WRITE(*,*)

  WRITE(*,*) 'ndmx = ', ndmx
  WRITE(*,*) 'zed = ', zed, ''
  WRITE(*,*) 'zval = ', zval
  WRITE(*,*) 'nwfx = ', size(oc)
  WRITE(*,*) 'nwf = ', nwf
  do i = 1,nwf
    write(*,*)
    write(*,'(1x,A,A)') 'el = ', el(i)
    write(*,'(1x,A,I5)') 'nn = ', nn(i)
    write(*,'(1x,A,I5)') 'll = ', ll(i)
    write(*,'(1x,A,F18.10)') 'oc = ', oc(i)
  enddo
  !write(*,*) 'enne = ', enne
  !write(*,*) 'enl = ', enl(1:nwf)

  write(*,*) 'nn = ', nn(1:nwf)
  write(*,*) 'll = ', ll(1:nwf)
  write(*,*) 'oc = ', oc(1:nwf)

  write(*,*)

  flush(6)

  !stop ! ffr


  !
  ! Compute an initial estimate of the potential
  !
  write(*,*) 'Before starting_potential:'
  write(*,*) 'v0   = ', v0(1:2)
  write(*,*) 'vxt  = ', vxt(1:2)
  write(*,*) 'vpot1 = ', vpot(1:2,1)
  write(*,*) 'vpot2 = ', vpot(1:2,2)

  write(*,*) 'nrmesh = ', grid%mesh

  call starting_potential( ndmx, grid%mesh, zval, zed, nwf, oc, nn, ll,&
                           grid%r, enl, v0, vxt, vpot, enne, nspin )
  
  write(*,*) 'After starting_potential: (Ha unit)'
  write(*,*) 'v0   = ', v0(1:2)*0.5
  write(*,*) 'vxt  = ', vxt(1:2)*0.5
  write(*,*) 'vpot1 = ', vpot(1:2,1)*0.5
  write(*,*) 'vpot2 = ', vpot(1:2,2)*0.5
  write(*,*) 'enl = ', enl(1:nwf)*0.5


  !
  ! allocate variables for SIC, if needed
  !
  IF( isic /= 0 ) THEN
    WRITE(*,*) 'Using SIC'
    ALLOCATE( vsic(ndmx,nwf), vsicnew(ndmx), vhn1(ndmx), egc(ndmx) )
    vsic = 0.0_dp
  ENDIF
  WRITE(*,*) 'ic = ', ic

  !stop ! ffr
  
  !
  ! solve the eigenvalue self-consistent equation
  !
  !call scf(ic)
  call ld1x_my_scf(ic)

  !
  ! compute relativistic corrections to the eigenvalues
  !
  IF( relpert ) call compute_relpert( evel, edar, eso )
  !
  ! compute total energy
  !
  CALL elsd(zed, grid, rho, vxt, vh, vxc, exc, excgga, nwf, nspin, enl, oc,    &
            etot, ekin, encl, ehrt, ecxc, evxt)
  !
  IF( verbosity=='high' ) call elsd_highv(ic)
  !
  ! add sic correction if needed
  !
  IF( isic /= 0 ) call esic()
  !
  !   print results
  !
  call write_results()
  !
  !  compute logarithmic derivative
  !
  IF( deld > 0.0_DP .and. ild ) call lderiv()
  !
  ! compute C6 coefficient if required
  !
  IF (vdw) THEN
    call c6_tfvw( grid%mesh, zed, grid, rho(1,1) )
    call c6_dft( grid%mesh, zed, grid )
  ENDIF
  !
  IF (isic /= 0) THEN
    DEALLOCATE(egc, vhn1, vsicnew, vsic)
  ENDIF
  !

  write(*,*)
  write(*,*) 'Exit all_electron'
  write(*,*)

  RETURN
  !
END SUBROUTINE all_electron
