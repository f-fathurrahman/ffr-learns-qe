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
  qvan2 calculate Q(G) for given indices m,n
  They are calculated on the fly? Not stored?

- addusdens
- s_psi
- qvan2

Two grids: dense and smooth.

`fft_interpolate` is called in subroutine `set_vrs` to interpolate potential
in dense grid to smooth grid. This smooth potential will be used in local potential
operator of Hamiltonian, to be multiplied with wave function in real space.

Formula in `qvan2`:
$$
Q_{mn}(\mathbf{G}) = \sum_{lm} (-\mathrm{i})^l
a^{lm}_{ij} Y_{lm}(\widehat{\mathbf{G}}) \, Q_{ij}^{l}(G)
$$
Formula for calculating $a^{lm}_{ij}$ are coefficients of expansions of the product of two real spherical harmonics into real spherical harmonic
$$
Y_{l_i m_i}(\hat{r}) \, Y_{l_j m_j}(\hat{r}) =
\sum_{lm}  a^{lm}_{l_i m_i l_j m_j} Y_{lm}(\hat{r})
$$
In code:

```fortran
Y_limi(r) * Y_ljmj(r) = \sum_LM  ap(LM,limi,ljmj)  Y_LM(r)
```

