!==============================================================================
!
! Routines:
!
! (1) absh() Originally By MLT Last Modified 7/1/2008 (JRD)
!
! Calculates absorption spectra in haydock code.
!
! input: mmts type (has the Haydock coefficients, a_n b_n)
! nspin number of spins
! eta energy resolution
!
! output: file "absorption_eh.dat"
!
! Calculate the dielectric function (real and imaginary parts),
! from coefficients mmts%an, mmts%bn
!
! omega = frequency, given in eV
! eta = energy broadening, given in eV
!
!==============================================================================
!The following macro puts any point/array in the [-0.5, 0.5) range:
!The following macro puts any point/array in the [0, 1) range:
!Integer division of a/b rounded up*/
!Rounds a up to the smallest multiple of b*/
! disable Fortran OMP pragmas if not -DOMP*/
! note: C standard does not permit $ in identifiers, however this seems acceptable
! as an extension, for all versions of cpp I tried. --DAS
! truncate spaces in string
!#!define TRUNC(s) trim(adjustl(s))
! Oracle compiler has a length limit of 132 characters and won`t support these macros
! No checking for faster performance, if not in debug mode
! Use this instead of the intrinsic 'deallocate' for pointers
! Use this instead of the intrinsic 'deallocate' for arrays
!the TOSTRING macro converts a macro into a string
! deprecated identifiers
! Created Sept 2011 by DAS.
! Define characteristics of various compilers, via compiler symbols (e.g. -DGNU)
! to be used directly from the arch.mk files, and then defining what we need to do
! for that compiler via the symbols for various properties (e.g. NOSIZEOF).
! Ideally, to support a new compiler, one need only change this file, adding a
! new block to define what -DNEWCOMPILER would mean.
! NOTE: of course, Makefile-level issues still need to be handled in common-rules.mk
! very ancient version may require NOSIZEOF
! FHJ: Support for Open64 will be removed shortly in favor of OpenUH
! open64 is very similar to path, it is an open-sourced version of it
! omp_lib.f90 needed to do OpenMP, see common-rules.mk.
! cce 7.4.4 and before support sizeof for intrinsic types, but need NOSIZEOF_TYPE
! cce 8.0.0 and later do not allow sizeof for multidimensional arrays, requiring us
! to turn sizeof off everywhere. Why would Cray do this?
! It is considered a bug in OPEN64 that sizeof will not work in our code.
! on some platforms there is a different return value for sizeof if build is 64-bit
! Intrinsic module for OpenMP. Almost all compilers that support OpenMP provide
! a "omp_lib.mod" module, though the OpenMP standard allow them to only ship a
! "omp_lib.h" Fortran header.
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
module absh_m
  use global_m
  implicit none
  private
  public :: &
    absh
contains
subroutine absh(mmts,nspin,nspinor,eta)
  type (mmtsinfo), intent(in) :: mmts
  integer, intent(in) :: nspin, nspinor
  real(DP), intent(in) :: eta
  integer :: ih,iw,nwstep
  real(DP) :: w0,wfin,omega,pref,eps1,eps2
  complex(DPC) :: sum,zo
  pref = 16.d0 * PI_D**2 /(mmts%vol * dble(nspin) * dble(nspinor))
 
  w0 = 0.d0
  wfin = 20.d0
  nwstep = 2000
  call open_file(10,file='absorption_eh.dat',form='formatted',status='replace')
  write(10,*) "# Column 1: omega"
  write(10,*) "# Column 2: eps2(omega)"
  write(10,*) "# Column 3: eps1(omega)"
  do iw=0,nwstep
    eps2 = 0.d0
    eps1 = 0.d0
    omega = w0 + (wfin - w0) * dble(iw) / dble(nwstep)
!------------------------------------
! Absorption contribution
    zo = cmplx(omega,eta,kind=DPC)/RYD
! Assume that all coefficients above mmts%nmax are constant...
!
! sum = sqrt( ((zo - mmts%an(mmts%nmax) )/2.d0)**2 - mmts%bn(mmts%nmax) )
!
! or zero... (it doesn`t matter if mmts%nmax is sufficiently high)
    sum = cmplx(0.d0,0.d0,kind=DPC)
    do ih= mmts%nmax, 2, -1
      sum = mmts%bn(ih-1)/( zo - mmts%an(ih) - sum)
    enddo
    sum = cmplx(mmts%norm / PI_D,0.d0,kind=DPC)/( zo - mmts%an(1) - sum)
    eps2 = -pref*aimag(sum)
    eps1 = -pref*dble(sum) + 1.d0
!-------------------------------------
! Emission contribution
    zo = cmplx(-omega,eta,kind=DPC)/RYD
! Assume that all coefficients above mmts%nmax are constant...
!
! sum = sqrt( ((zo - mmts%an(mmts%nmax) )/2.d0)**2 -
! > mmts%bn(mmts%nmax) )
!
! or zero...
    sum = cmplx(0.d0,0.d0,kind=DPC)
    do ih= mmts%nmax, 2, -1
      sum = mmts%bn(ih-1)/( zo - mmts%an(ih) - sum)
    enddo
    sum = cmplx(mmts%norm / PI_D,0.d0,kind=DPC)/( zo - mmts%an(1) - sum)
    eps2 = pref*aimag(sum) + eps2
    eps1 = -pref*dble(sum) + eps1
    write(10,100) omega,eps2,eps1
  enddo
100 format(3f16.9)
  call close_file(10)
 
  return
end subroutine absh
end module absh_m
