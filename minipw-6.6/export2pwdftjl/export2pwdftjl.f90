! Using json-fortran module to export various data from PWSCF
! 
! Intended to be used for data that cannot be accessed easily from
! Julia, especially, the ones that defined using derived data types.
! Other varibles such as scalars and arrays that are exposed directly
! in the .so file should be accessed directly using cglobal and
! other C-interface functions of Julia.
!
! Current limitations of json-fortran:
!
! - Multidimensional arrays are not supported.
!   Work-around: using reshape to one-dimensional array and
!   also its shape
!
! - Complex scalars and arrays are not supported
!   Work-around: write real and imag parts separately


INCLUDE 'prepare_all.f90'

PROGRAM main

  CALL prepare_all()
  
  CALL export_atoms()
  
  CALL export_pspot_upf()

  CALL export_atom_grid()

  CALL export_vlocal()

  CALL export_vars_uspp_mod()

  CALL export_paw_radial_integrator()

  CALL export_vars_scf_mod()

END PROGRAM main


include 'export_atoms.f90'
include 'export_pspot_upf.f90'
include 'export_atom_grid.f90'
include 'export_vlocal.f90'
include 'export_paw_radial_integrator.f90'
include 'export_vars_uspp_mod.f90'
include 'export_vars_scf_mod.f90'
