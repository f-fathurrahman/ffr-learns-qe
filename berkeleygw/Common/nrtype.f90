!>================================================================================
!!
!! Modules:
!!
!! (1) nrtype_m Originally By ? Last Modified 08/2019 (FHJ)
!!
!! Global constants, parameters, and routines that convert numeric types.
!!
!!================================================================================
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
module nrtype_m
  use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
  implicit none
  public ! only public parameters here
! Below are the version number for each BerkeleyGW file. These numbers should
! be changed whenever the structure of the file is altered and there`s either
! incompatibility with the previous version or a new feature. The version
! should be -1 if that file is not versioned yet, which corresponds to
! the formats used in the Berkeley 1.0.x family.
  integer, parameter :: VER_WFN_FORT = -1
  integer, parameter :: VER_WFN_HDF5 = 1
  integer, parameter :: VER_WFN = VER_WFN_FORT
  integer, parameter :: VER_EPS_FORT = -1
  integer, parameter :: VER_EPS_HDF5 = 3
  integer, parameter :: VER_EPS = VER_EPS_FORT
  integer, parameter :: VER_BSE_FORT = 1
  integer, parameter :: VER_BSE_HDF5 = 2
  integer, parameter :: VER_BSE = VER_BSE_FORT
! Screening types, typically stored in a variable named `iscreen`
  integer, parameter :: SCREEN_SEMICOND = 0 !< semiconductor screening
  integer, parameter :: SCREEN_GRAPHENE = 1 !< graphene screening
  integer, parameter :: SCREEN_METAL = 2 !< metal screening
! Truncation flags, typically stored in a variable named `itrunc` or `icutv`
  integer, parameter :: TRUNC_NONE = 0 !< no truncation (3D)
  integer, parameter :: TRUNC_SPHERICAL = 2 !< 0D spherical truncation
  integer, parameter :: TRUNC_WIRE = 4 !< 1D cell wire truncation
  integer, parameter :: TRUNC_BOX = 5 !< 0D cell box truncation
  integer, parameter :: TRUNC_SLAB = 6 !< 2D slab truncation
  integer, parameter :: TRUNC_SUPERCELL = 7 !< supercell truncation (3D), experimental
!> Maximum number of bands supported by the *inread* routines. This sets the
!! size of arrays such as "occupations". These arrays should all be allocated
!! dynamically in the future.
  integer, parameter :: MAX_BANDS = 1000000 ! "occupations" array => 7MB
!> Maximum number of {k,q}-points supported by the *inread* routines.
!! The actual number of k-points/q-points in the WFN/bsemat/epsmat files
!! can be larger.
  integer, parameter :: MAX_KPTS = 100000 ! "kpt_read" array => 0.8 MB
!> parameters for real-space resolution in cell-truncation schemes
  integer, parameter :: n_in_box = 2
  integer, parameter :: n_in_wire = 4
!> parameter for construction of Wigner-Seitz cell
  integer, parameter :: ncell = 3
!> number of Monte-Carlo integration points
  integer, parameter :: nmc_coarse = 250000
  integer, parameter :: nmc_fine = 2500000
  integer, parameter :: nmc = nmc_fine
!> type definitions following the convention of Numerical Recipes
!! do not ever use single-precision!!
! integer, parameter :: SP = kind(1.0)
  integer, parameter :: DP = kind(1.0d0)
! integer, parameter :: SPC = kind((1.0,1.0))
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  integer, parameter :: size_of_complex_DPC= DPC * 2
  integer, parameter :: size_of_SCALAR = DP
!> a shift on the grid in order to avoid the singularity for truncation
  real(DP), parameter :: trunc_shift(3) = (/0.5d0, 0.5d0, 0.5d0/)
!> physical constants
!!
!! These are the "2010 CODATA recommended values" taken from
!! "The NIST Reference on Constants, Units, and Uncertainty"
!! http://physics.nist.gov/cuu/
!!
!! The following variables are used throughout the package:
!! 'BOHR', 'bohr' is Bohr radius, in Angstrom
!! 'RYD', 'ryd2eV', 'rydberg' is Rydberg constant times hc, in eV
!! 'HARTREE', 'hartree' is Hartree energy, in eV
!! 'LIGHTSPEED' is inverse alpha (fine-structure constant)
!!
!! These variables are defined in the following files:
!! Common/nrtype.f90
!! Common/wfn_utils.cpp
!! MeanField/EPM/ff2vq.py
!! MeanField/EPM/sysParams.f90
!! MeanField/EPM/vca.py
!! MeanField/ICM/icm.cpp
!! Visual/common.py
!!
  real(DP), parameter :: BOHR = 0.52917721092_dp
  real(DP), parameter :: RYD = 13.60569253_dp
  real(DP), parameter :: LIGHTSPEED = 137.035999074_dp
!> mathematical constants
!! real(SP), parameter :: PI_S = 3.1415926535897932384626433832795_sp
  real(DP), parameter :: PI_D = 3.1415926535897932384626433832795_dp
  real(DP), parameter :: TOL_Small = 1.0d-6
  real(DP), parameter :: TOL_Zero = 1.0d-12
  real(DP), parameter :: TOL_Degeneracy = 1.0d-6
  real(DP), parameter :: INF = 1.0d12
!> Do direct diagonalization for BSE
  integer, parameter :: BSE_ALGO_DIAG = 1
!> Solve BSE with Lanczos alg. by M. Shao and C. Yang.
  integer, parameter :: BSE_ALGO_LANCZOS = 2
!> Solve BSE with Haydock scheme. Only works with TDA.
  integer, parameter :: BSE_ALGO_HAYDOCK = 3
!> Iterative diagonalization of BSE with PRIMME solver
  integer, parameter :: BSE_ALGO_DIAG_PRIMME = 4
contains
  !> Return a pointer that is a complex reinterpretation of the double precision
  !! input array data_in. No data is copied if the input array is contiguous.
  !! `siz` is the total number of elements of the input array.
  function ptr_complex2real(data_in, siz) result(data_out)
    complex(DPC), intent(inout), target :: data_in(*)
    integer, intent(in) :: siz
    real(DP), pointer :: data_out(:)
    call c_f_pointer(c_loc(data_in), data_out, [2*siz])
  end function ptr_complex2real
  !> Return a pointer that is a double precision reinterpretation of the complex
  !! input array data_in. No data is copied if the input array is contiguous.
  !! `siz` is the total number of elements of the input array.
  function ptr_real2complex(data_in, siz) result(data_out)
    real(DP), intent(inout), target :: data_in(*)
    integer, intent(in) :: siz
    complex(DPC), pointer :: data_out(:)
    call c_f_pointer(c_loc(data_in), data_out, [siz/2])
  end function ptr_real2complex
end module nrtype_m
