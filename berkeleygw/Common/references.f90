!=============================================================================
!
! Modules:
!
! references_m Originally By FHJ
!
! Reference manager for BerkeleyGW. File references.f90p automatically parses
! references.txt for references and generates references.f90 (via
! mako_preprocess.py).
!
!=============================================================================
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
module references_m
  use, intrinsic :: iso_c_binding
  use message_m
  use nrtype_m
  use os_m
  use peinfo_m
  use push_pop_m
  implicit none
  private
  public :: require_reference, show_references
  integer, parameter, public :: &
    REF_DESLIPPE2012 = 1, &
    REF_HYBERTSEN1986 = 2, &
    REF_ROHLFING2000 = 3, &
    REF_ISMAILBEIGI2016 = 4, &
    REF_SAMSONIDZE2011 = 5, &
    REF_DESLIPPE2013 = 6, &
    REF_LIU2015 = 7, &
    REF_QIU2015 = 8, &
    REF_SHAO2016 = 9, &
    REF_JORNADA2017 = 10, &
    REF_SHAO2018 = 11, &
    REF_DELBEN2019SUBSPACE = 12, &
    REF_DELBEN2019LARGESCALE = 13, &
    REF_JORNADA2020 = 14, &
    REF_LI2019 = 15, &
    REF_BARKER2018 = 16, &
    REF_WU2020 = 17, &
    REF_WU2019 = 18, &
    REF_ZHAO2021 = 19
  integer, private :: i_
  logical, private :: need_reference(19) = [(.false., i_=1,19)]
contains
subroutine require_reference(reference)
  integer, intent(in) :: reference
 
  if (reference<1 .or. reference>19) then
    write(6,'(/1x,a,i0,a/)') &
      'WARNING: Internal reference error, reference ', reference, ' not found.'
  else
    need_reference(reference) = .true.
  endif
 
end subroutine require_reference
subroutine show_references(iunit)
  integer, intent(in), optional :: iunit
  integer :: iunit_
  integer :: iref
  iunit_ = 6
  if (present(iunit)) iunit_ = iunit
 
  if (peinf%inode/=0) then
   
    return
  endif
  if (any(need_reference)) then
    write(iunit_,'(/a)') "********************************************************************************"
    write(iunit_,'(a)') "*                                                                              *"
    write(iunit_,'(a)') "* Your calculation employed methods and algorithms published in peer-reviewed  *"
    write(iunit_,'(a)') "* papers. Please, cite the following references to give proper credit to their *"
    write(iunit_,'(a)') "*          authors and help support future development of BerkeleyGW.          *"
    write(iunit_,'(a)') "*                                                                              *"
    write(iunit_,'(a/)') "********************************************************************************"
  endif
  if (need_reference(1)) then
    write(iunit_,'(1x,a)') "J. Deslippe, G. Samsonidze, D. A. Strubbe, M. Jain, M. L. Cohen, and S. G."
    write(iunit_,'(1x,a)') "Louie, BerkeleyGW: A Massively Parallel Computer Package for the Calculation"
    write(iunit_,'(1x,a)') "of the Quasiparticle and Optical Properties of Materials and Nanostructures,"
    write(iunit_,'(1x,a)') "Computer Physics Communications 183, 1269 (2012)."
    write(iunit_,'()')
  endif
  if (need_reference(2)) then
    write(iunit_,'(1x,a)') "M. S. Hybertsen and S. G. Louie, Electron Correlation in Semiconductors and"
    write(iunit_,'(1x,a)') "Insulators: Band Gaps and Quasiparticle Energies, Phys. Rev. B 34, 5390"
    write(iunit_,'(1x,a)') "(1986)."
    write(iunit_,'()')
  endif
  if (need_reference(3)) then
    write(iunit_,'(1x,a)') "M. Rohlfing and S. G. Louie, Electron-Hole Excitations and Optical Spectra"
    write(iunit_,'(1x,a)') "from First Principles, Phys. Rev. B 62, 4927 (2000)."
    write(iunit_,'()')
  endif
  if (need_reference(4)) then
    write(iunit_,'(1x,a)') "S. Ismail-Beigi, Truncation of Periodic Image Interactions for Confined"
    write(iunit_,'(1x,a)') "Systems, Phys. Rev. B 73, 233103 (2006)."
    write(iunit_,'()')
  endif
  if (need_reference(5)) then
    write(iunit_,'(1x,a)') "G. Samsonidze, M. Jain, J. Deslippe, M. L. Cohen, and S. G. Louie, Simple"
    write(iunit_,'(1x,a)') "Approximate Physical Orbitals For GW Quasiparticle Calculations, Phys. Rev."
    write(iunit_,'(1x,a)') "Lett. 107, 186404 (2011)."
    write(iunit_,'()')
  endif
  if (need_reference(6)) then
    write(iunit_,'(1x,a)') "J. Deslippe, G. Samsonidze, M. Jain, M. L. Cohen, and S. G. Louie, Coulomb-"
    write(iunit_,'(1x,a)') "Hole Summations and Energies For GW calculations with Limited Number of Empty"
    write(iunit_,'(1x,a)') "Orbitals: A Modified Static Remainder Approach, Phys. Rev. B 87, 165124"
    write(iunit_,'(1x,a)') "(2013)."
    write(iunit_,'()')
  endif
  if (need_reference(7)) then
    write(iunit_,'(1x,a)') "F. Liu, L. Lin, D. Vigil-Fowler, J. Lischner, A. F. Kemper, S. Sharifzadeh, F."
    write(iunit_,'(1x,a)') "H. da Jornada, J. Deslippe, C. Yang, J. B. Neaton, and S. G. Louie, Numerical"
    write(iunit_,'(1x,a)') "Integration for Ab Initio Many-Electron Self Energy Calculations within the GW"
    write(iunit_,'(1x,a)') "Approximation, Journal of Computational Physics 286, 1 (2015)."
    write(iunit_,'()')
  endif
  if (need_reference(8)) then
    write(iunit_,'(1x,a)') "D. Y. Qiu, T. Cao, and S. G. Louie, Nonanalyticity, Valley Quantum Phases, and"
    write(iunit_,'(1x,a)') "Lightlike Exciton Dispersion in Monolayer Transition Metal Dichalcogenides:"
    write(iunit_,'(1x,a)') "Theory and First-Principles Calculations, Phys. Rev. Lett. 115, (2015)."
    write(iunit_,'()')
  endif
  if (need_reference(9)) then
    write(iunit_,'(1x,a)') "M. Shao, F. H. da Jornada, C. Yang, J. Deslippe, and S. G. Louie, Structure"
    write(iunit_,'(1x,a)') "Preserving Parallel Algorithms for Solving the Bethe-Salpeter Eigenvalue"
    write(iunit_,'(1x,a)') "Problem, Linear Algebra and Its Applications 488, 148 (2016)."
    write(iunit_,'()')
  endif
  if (need_reference(10)) then
    write(iunit_,'(1x,a)') "F. H. da Jornada, D. Y. Qiu, and S. G. Louie, Nonuniform Sampling Schemes of"
    write(iunit_,'(1x,a)') "the Brillouin Zone for Many-Electron Perturbation-Theory Calculations in"
    write(iunit_,'(1x,a)') "Reduced Dimensionality, Phys. Rev. B 95, 035109 (2017)."
    write(iunit_,'()')
  endif
  if (need_reference(11)) then
    write(iunit_,'(1x,a)') "M. Shao, F. H. da Jornada, L. Lin, C. Yang, J. Deslippe, and S. G. Louie, A"
    write(iunit_,'(1x,a)') "Structure Preserving Lanczos Algorithm for Computing the Optical Absorption"
    write(iunit_,'(1x,a)') "Spectrum, SIAM J. Matrix Anal. & Appl. 39, 683 (2018)."
    write(iunit_,'()')
  endif
  if (need_reference(12)) then
    write(iunit_,'(1x,a)') "M. Del Ben, F. H. da Jornada, G. Antonius, T. Rangel, S. G. Louie, J."
    write(iunit_,'(1x,a)') "Deslippe, and A. Canning, Static Subspace Approximation for the Evaluation of"
    write(iunit_,'(1x,a)') "G0W0 Quasiparticle Energies within a Sum-over-Bands Approach, Phys. Rev. B 99,"
    write(iunit_,'(1x,a)') "125128 (2019)."
    write(iunit_,'()')
  endif
  if (need_reference(13)) then
    write(iunit_,'(1x,a)') "M. Del Ben, F. H. da Jornada, A. Canning, N. Wichmann, K. Raman, R. Sasanka,"
    write(iunit_,'(1x,a)') "C. Yang, S. G. Louie, and J. Deslippe, Large-Scale GW Calculations on Pre-"
    write(iunit_,'(1x,a)') "Exascale HPC Systems, Computer Physics Communications 235, 187 (2019)."
    write(iunit_,'()')
  endif
  if (need_reference(14)) then
    write(iunit_,'(1x,a)') "F. H. da Jornada, L. Xian, A. Rubio, and S. G. Louie, Universal Slow Plasmons"
    write(iunit_,'(1x,a)') "and Giant Field Enhancement in Atomically Thin Quasi-Two-Dimensional Metals,"
    write(iunit_,'(1x,a)') "Nat. Commun. 11, 1013 (2020)."
    write(iunit_,'()')
  endif
  if (need_reference(15)) then
    write(iunit_,'(1x,a)') "Z. Li, G. Antonius, M. Wu, F. H. da Jornada, and S. G. Louie, Electron-Phonon"
    write(iunit_,'(1x,a)') "Coupling from Ab Initio Linear-Response Theory within the GW Method:"
    write(iunit_,'(1x,a)') "Correlation-Enhanced Interactions and Superconductivity in Ba1-xKxBiO3, Phys."
    write(iunit_,'(1x,a)') "Rev. Lett. 122, 186402 (2019)."
    write(iunit_,'()')
  endif
  if (need_reference(16)) then
    write(iunit_,'(1x,a)') "B. A. Barker, and S. G. Louie, Electronic and Optical Properties of Solids"
    write(iunit_,'(1x,a)') "with Strong Spin-Orbit Coupling, U. C. Berkeley PhD dissertation (2018)."
    write(iunit_,'()')
  endif
  if (need_reference(17)) then
    write(iunit_,'(1x,a)') "M. Wu, and S. G. Louie, Spin-Orbit Coupling, Broken Time-Reversal Symmetry, "
    write(iunit_,'(1x,a)') "and Polarizability Self-Consistency in GW and GW-BSE Theory with Applications"
    write(iunit_,'(1x,a)') "to Two-Dimensional Materials, U. C. Berkeley PhD dissertation (2020)."
    write(iunit_,'()')
  endif
  if (need_reference(18)) then
    write(iunit_,'(1x,a)') "M. Wu, Z. Li, T. Cao and S. G. Louie, Physical Origin of Giant Excitonic"
    write(iunit_,'(1x,a)') "and Magneto-Optical Responses in Two-Dimensional Ferromagnetic Insulators,"
    write(iunit_,'(1x,a)') "Nat. Commun. 10, 2371 (2019)."
    write(iunit_,'()')
  endif
  if (need_reference(19)) then
    write(iunit_,'(1x,a)') "F. Zhao, and S. G. Louie, Topological Electronic Properties and Optical"
    write(iunit_,'(1x,a)') "Properties of Graphene Nanoribbons and Carbon Conjugated Systems,"
    write(iunit_,'(1x,a)') "U. C. Berkeley PhD dissertation (2021)."
    write(iunit_,'()')
  endif
  if (any(need_reference)) then
    write(iunit_,'(a/)') "********************************************************************************"
  endif
 
end subroutine show_references
end module references_m
