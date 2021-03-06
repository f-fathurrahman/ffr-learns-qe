# 3 April 2017


I tried to learn how to apply non-local pseudopotential to wavefunction.

In `h_psi`, the relevant calls seem to be `calbec()` and `add_vuspsi()`,
for non-real space algorithm.

```fortran
CALL calbec( n, vkb, psi, becp, m )
CALL add_vuspsi( lda, n, m, hpsi )
```

Array `vkb` is initialized by calling `init_us_2()`

`calbec()` is defined in `Modules/becmod.f90`.

```fortran
SUBROUTINE calbec_bec_type( npw, beta, psi, betapsi, nbnd )
SUBROUTINE calbec_k( npw, beta, psi, betapsi, nbnd )
```


# 25 December 2016

Updating to qe-6.0. (Now it use the name `qe` instead of `espresso`)

Uploaded to github.

Changed the the name to `ffr-learns-qe`.

Array `igk` (formerly in module `wvfct` file: `pwcom.f90`) is not used
anymore.

The subroutine `realspace_grids_init()` is not used anymore.


# 19 May 2016

Updating to espresso-5.4.0

* No flib directory now. Its functionality is replace by LAXlib

* New GVECW module to contain the variables ecutwfc, ecfixed, etc
  So, any statements that contain USE wvfct, ONLY : ecutwfc will result
  in error.


# 21 Feb 2016


Updating to espresso-5.3.0.

* FFT subroutines are isolated to separate directory: FFTXlib
  Include and library paths must be changed accordingly.

* All calls to flush_unit are replaced by FLUSH( stdout )


Project: implementing steepest descent (SD) method to solve Kohn-Sham
equations.

We need to evaluate total energy, given wavefunction coefficients
as the input. Most likely I need to learn the functions and subroutine needed to
calculate the total energy.

Currently, I don't know yet where the Poisson solver need to be
updated.
Probably the it will be after doing one SD update.

How about the orthogonality problem?

****

Alternatively, the DFT++ formalism may be used.
One initial important task that is needed is to implement various
operators.
