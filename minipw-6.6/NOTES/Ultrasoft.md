Variables:
- okvan
- qrad (combined qfuncl?)
- qq_at, qq_nt
- dvan
- ddeq

Many of these variables as allocated in allocate_nlpot

initialization: init_us_1:

Subroutines:
- compute_qrad: (in init_us_1.f90)
- qvan2:
  The interpolation table for the radial Fourier transform is stored
  in qrad
- addusdens
- s_psi
- qvan2

Two grids: dense and smooth.



`fft_interpolate` is called in subroutine `set_vrs` to interpolate potential
in dense grid to smooth grid. This smooth potential will be used in local potential
operator of Hamiltonian, to be multiplied with wave function in real space.

