Nonlinear core correction
=========================

Important variables:
- upf(nt)%nlcc
- upf(nt)%rho_atc
- rho_core, rhog_core (module scf)

```fortran
REAL(DP), ALLOCATABLE :: rho_core(:) !! the core charge in real space
COMPLEX(DP), ALLOCATABLE :: rhog_core(:) !! the core charge in reciprocal space
```

Important routines:

- set_rhoc
- drhoc

set_rhoc is called in hinit0.
set_rhoc is called after initializing local pseudopotentials.

rho_core and rhog_core affect Hartree and XC terms in v_of_rho.

energy, force, and stress terms should be modified.