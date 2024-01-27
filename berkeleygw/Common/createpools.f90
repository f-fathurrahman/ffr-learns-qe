!===============================================================================
!
! Routines:
!
! 1. createpools() Originally By gsm Last Modified 06/11/2012 (DVF)
!
! Create pools of valence bands/diagonal matrix elements for Epsilon/Sigma. Number of pools
! is chosen to minimize memory in calculation. This is also used in BSE, which has multiple
! options for what is pooled based on the parameters of the problem.
!
!===============================================================================
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
module createpools_m
  use global_m
  implicit none
  private
  public :: createpools
contains
  !> Create pools of valence bands (for epsilon) or diagonal matrix elements (for sigma). For BSE,
  !! the description starts at `Finally,` below. You will probably still need to reference the
  !! information above the `Finally,`. The description before `Finally,` is for Sigma and Epsilon.
  !!
  !! Each pool is defined as a group of *valence* bands or diagonal matrix elements.
  !! The behaviour of the code is to minimize the memory requirement.
  !! For epsilon, in a naive algorithm, if we distributed all (v,c) pairs round-robin,
  !! each process could potentially get as many as (Nk*)Nc*Nv/Nproc
  !! bands, which can be pretty large. The usage of pools tackles this
  !! memory issue. The same holds for sigma, except Nv -> Ndiag, Nc -> Nband, where Ndiag
  !! is the number of diagonal bands for which you are computing sigma (offdiagonal calculations
  !! do not use this routine), and Nband is the number of bands in the CH summation (`summation bands`).
  !! Obviously, (v,c) pairs goes to (diagonal matrix elements,summation bands) pairs in the above description.
  !! The current algorithm works in two levels (this is for epsilon; for sigma make the replacements
  !! valence bands -> diagonal matrix elements, conduction bands -> summation bands, Nv -> Ndiag):
  !! (1) First, we divide the valence bands into groups (the "pools")
  !! (2) Then, the conduction bands are spread across the MPI processes.
  !! If the pools don`t have all the same size (i.e., if Nv is not
  !! divisible by Nproc), then some MPI processes will be idle.
  !! Variable description: formal name first, then actual name in epsilon and sigma in parentheses (epsilon first).
  !! The variable name npq means number of pooled quantity (valence bands or diagonal matrix elements),
  !! while nsq means number of spread quantity (conduction bands or summation bands), in the sense above.
  !! If there is better technical language for describing this, feel free to update this description.
  !! \param npq num. of pooled quantity (valence bands or diagonal matrix elements)
  !! \param nsq num. of spread quantity (conduction or summation bands)
  !! \param npes num. of MPI processes
  !! \param npoolsout num. of (valence band or diagonal matrix element) pools created
  !! \param npqownmaxout max num. of pool quantity (valence bands or diagonal matrix elements)/MPI process
  !! \param nsqwnmaxout max num. of spread quantity (conduction or summation bands)/MPI process
  !!
  !! Finally, this routine is also used in BSE, where the pooled quantity and spread quantity depends on the
  !! parameters of the problem. See BSE/distrib_kernel.f90 . The choices are 1) pooled quantity and spread
  !! quantity = num. valence bands, 2) pooled quantity and spread quantity = num. cond. bands, and
  !! 3) pooled quantity and spread quantity = num. k-points . The above description changes accordingly.
  subroutine createpools(npq,nsq,npes,npoolsout,npqownmaxout,nsqownmaxout)
    integer, intent(in) :: npq,nsq,npes
    integer, intent(out) :: npoolsout,npqownmaxout,nsqownmaxout
    integer :: nmemmin
    integer :: npools,nsqownmax,npqownmax,nmem
    integer :: npes_per_npools
   
    nmemmin = npq + nsq + 1
    npoolsout = 0
    npqownmaxout = 0
    nsqownmaxout = 0
    ! FHJ: we can`t have more pools than the number of pool quantities or MPI processes
    do npools = 1, min(npq,npes)
      ! FHJ: npqownmax = max num. of pooled quantity per proc. for this number of pools
      npqownmax = (npq + npools - 1) / npools
      ! FHJ: nsqownmax = max num. of spread quantity per proc. for this number of pools
      ! For a given pool, we have (npes/npools) MPI processes over which the
      ! pooled quantity is parallelized.
      npes_per_npools = npes/npools
      nsqownmax = (nsq + npes_per_npools - 1) / npes_per_npools
      ! FHJ: Note that both npqownmax and nsqownmax are max. number of the pooled/spread quantity, and for
      ! this reason we round the divisions up. We also round (npes/npools)
      ! down, so we are conservative (and pessimistic) about the resources.
      ! FHJ: Upper bound on the memory required for each proc to store the WFNs
      nmem = npqownmax + nsqownmax
      !ntime = npqownmax * nsqownmax
      if (nmem .lt. nmemmin) then
        nmemmin = nmem
        npoolsout = npools
        npqownmaxout = npqownmax
        nsqownmaxout = nsqownmax
      endif
    enddo
   
    return
  end subroutine createpools
end module createpools_m
