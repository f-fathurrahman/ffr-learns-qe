!
! Copyright (C) 2004-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
PROGRAM ld1
  !---------------------------------------------------------------
  !
  !     atomic self-consistent local-density program
  !     atomic rydberg units are used : e^2=2, m=1/2, hbar=1
  !     psi(r) = rR(r), where R(r) is the radial part of the wfct
  !     rho(r) = psi(r)^2 => rho(r) = (true charge density)*(4\pi r^2)
  !                       The same applies to the core charge
  !---------------------------------------------------------------
  !
  USE mp_global,         ONLY : mp_startup, mp_global_end
  USE environment,       ONLY : environment_start
  USE ld1inc,            ONLY : iswitch, write_coulomb, grid
  USE radial_grids,      ONLY : deallocate_radial_grid
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  CHARACTER(LEN=9) :: code = 'LD1'
  !
  !   write initialization information
  !
  CALL mp_startup()
  CALL environment_start(code)
  !
  !    read input, possible pseudopotential and set the main variables
  !
  CALL ld1_readin(input_file_)
  CALL ld1_setup()
  !
  !   four possible working mode:
  !
  IF( iswitch == 1 ) THEN
    !
    !   all-electron calculation
    !
    CALL my_all_electron(.true., 1) ! also compute log-deriv
    IF( write_coulomb ) CALL write_ae_pseudo ( )
    !
  ELSEIF( iswitch == 2 ) THEN
     !
     !   pseudopotential test
     !
     CALL run_test()
     CALL ld1_writeout()
     !
  ELSEIF( iswitch == 3 ) THEN
     !
     !  pseudopotential generation and test
     !
     CALL all_electron(.FALSE., 1) ! do not compute log-deriv
     CALL my_gener_pseudo()
     !if(.not. lgipaw_reconstruction) 
     CALL my_run_test()
     CALL ld1_writeout()
     !
  ELSEIF( iswitch == 4 ) THEN
     !
     ! LDA-1/2 correction to the input pseudopotential 
     !
     CALL run_lda_half( )
     CALL ld1_writeout( )
     !
  ELSE 
     CALL errore('ld1', 'iswitch not implemented',1)
  ENDIF 
  CALL deallocate_radial_grid( grid )

  CALL mp_global_end()

END PROGRAM ld1

