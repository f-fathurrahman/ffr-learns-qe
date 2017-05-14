! ffr: 5 Dec 2015
!
! An example of how to initialize G-vectors

!=-----------------------------------------------------------------------------
PROGRAM main
!=-----------------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE io_global, ONLY: stdout
  USE mp_global, ONLY : mp_startup, mp_global_end
  IMPLICIT NONE

  CALL mp_startup()

  CALL setup_structure_Si8()

  CALL setup_fft( 30_DP, 120_DP )

  CALL setup_symmetry()

  CALL setup_kpoints_noshift( 3, 3, 3 )

  CALL setup_gvect( .FALSE. )

  CALL mp_global_end()

END PROGRAM

