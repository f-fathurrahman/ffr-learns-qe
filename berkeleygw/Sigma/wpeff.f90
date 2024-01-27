!=============================================================================
!
! Routines:
!
! (1) wpeff() Originally By SIB Last Modified 5/1/2008 (JRD)
!
! This routine computes effective plasma frequencies (squared).
!
! Specifically, for each G` (igp) it computes for all G (ig),
! the quantity in formula (31) on page 5396 of
! Hybertsen and Louie, PRB vol 34, #8 (1986) given by
!
! Omega(G,G`)^2 = wp^2 * [rho(G-G`)/rho(0)] * (q+G).(q+G`)*vc(q+G)/(8pi)
!
! (vc(q+G) is the Coulomb interaction), and places the numbers
! into wpmtx(ig). Units are eV^2.
!
! Computes effective plasma freq matrix-
! wp(g,g`) = wp**2 (rho(g-g`)/rho(0)) (q+g).(q+g`) vc(q+g)/8pi
!
! isrtrq index array, (rq+g)**2
! nelec,wpsq rho(g=0), feg plasma frequency^2
! rho rho(g)
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
module wpeff_m
  use global_m
  use misc_m
  implicit none
  private
  public :: &
    wpeff
contains
subroutine wpeff(crys, gvec, wpg, sig, neps, isrtrq, igp, ncouls, wpmtx, nspin, qk,vcoul,coulfact)
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (wpgen), intent(in) :: wpg
  type (siginfo), intent(in) :: sig
  integer, intent(in) :: neps
  integer, intent(in) :: isrtrq(neps), igp, ncouls
  real(DP), intent(out) :: wpmtx(neps)
  integer, intent(in) :: nspin
  real(DP), intent(in) :: qk(3), vcoul(ncouls), coulfact
  real(DP) :: rho_g_minus_gp
  integer :: ig, igadd, igpadd, kadd, gg(3)
  real(DP) :: qg(3), qgp(3)
  logical, save :: warned=.false.
!-----------------------------
! SIB: stuff that is remembered between calls in order to speed
! things up
! 'save' is deprecated, replace with an object... DAS
  !logical, save :: first_call = .true.
  !logical, save :: q_is_not_zero
  !real(DP), save :: qk_old(3), fact, qnorm2
  !real(DP), allocatable, save :: precalc(:,:)
  real(DP) :: fact, qnorm2
  logical :: q_is_not_zero
  real(DP), allocatable :: precalc(:,:)
  ! no push_sub, called too frequently
! WARNING: if Fermi level is adjusted, the charge density here does not change,
! which will make the plasma frequency inconsistent. -- DAS
  if((abs(sig%efermi_input) > TOL_Zero .or. .not. sig%rfermi) &
    .and. (peinf%inode == 0) .and..not.warned) then
    write(0,'(a)') "WARNING: GPP plasma frequency does not respond to manual adjustment of the Fermi level."
    warned = .true.
  endif
!--------------- Begin Calculation ------------------------------------------------------
!------------------------
! SIB: Calculate stuff that is the same between calls for a given qk.
! The precalculated table contains:
!
! precalc(:,ig) = wp**2/rho(0)*(vc(q+G)/8pi)*(q+G)
  !if (sum(abs(qk-qk_old)) > 1.0d-12) then
    allocate(precalc (3,ncouls))
    precalc = 0.0d0
! Common factor
    fact=(wpg%wpsq(1)+(nspin-1)*wpg%wpsq(2))/(wpg%nelec(1)+(nspin-1)*wpg%nelec(2))
! Square length of qk
    qnorm2=dot_product(qk,matmul(crys%bdot,qk))
! Set a flag if q is not "zero"
    if (qnorm2 > 1.0d-12) then
      q_is_not_zero = .true.
    else
      q_is_not_zero = .false.
    endif
!---------- Loop over g and tabulate -------------------------------
    do ig=1,ncouls
! Get q+g and |q+g|^2
      igadd=isrtrq(ig)
! If g<>0 or q<>0, we just calculate the formula
      if (igadd.ne.1 .or. q_is_not_zero) then
        qg(:)=gvec%components(:,igadd)+qk(:)
        precalc(:,ig) = fact*vcoul(ig)*qg(:)/coulfact
      else
! If g=q=0, we have to avoid dividing by zero;
! we handle this special case separately below
        precalc(:,ig) = 0.0d0
      endif
    enddo ! g loop
  !endif
!----------------- Here starts the main calculation! ----------------------------
! Get q+gp
  igpadd=isrtrq(igp)
  qgp(:)=gvec%components(:,igpadd)+qk(:)
!!----------- Loop over g and calculate Omega^2(g,gp) -------------------------
  do ig=1,ncouls
    igadd=isrtrq(ig)
    wpmtx(ig) = 0.0d0
! Compute address of g-gp, and if it is a vector for which
! we have the density, get rho(g-gp); if out of bounds,
! skip this g
    gg(:) = gvec%components(:, igadd) - gvec%components(:, igpadd)
    call findvector(kadd, gg, gvec)
    if(kadd.eq.0) cycle
    rho_g_minus_gp = sum(wpg%rho(kadd,1:nspin))
! Using precalculated stuff, assemble together wpmtx
! if g<>0 or q<>0, we can just do what the formula says
    if (igadd.ne.1 .or. q_is_not_zero) then
      wpmtx(ig) = dot_product(qgp(:),matmul(crys%bdot,precalc(:,ig)))*rho_g_minus_gp
! The special case q=g=0
    else
! If gp=0, then wpmtx=wp**2 times a possible
! (1-cos(0))=0 if truncating
      if (igpadd.eq.1) then
        wpmtx(ig) = fact*rho_g_minus_gp
! if (sig%icutv/=TRUNC_NONE) wpmtx(ig) = 0.0d0
! When q=g=0, gp<>0, the result is just set to zero
! (it is infinite for true coulomb interaction and zero
! if truncated).
      else
        wpmtx(ig) = 0.0d0
      endif
    endif
  enddo !! g loop
  if(allocated(precalc))then;deallocate(precalc);endif
  return
end subroutine wpeff
end module wpeff_m
