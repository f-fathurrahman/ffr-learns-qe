Variables:
- okvan


Subroutines:
- addusdens
- s_psi
- qvan2

Two grids: dense and smooth.



`fft_interpolate` is called in subroutine `set_vrs` to interpolate potential
in dense grid to smooth grid. This smooth potential will be used in local potential
operator of Hamiltonian, to be multiplied with wave function in real space.

