!===============================================================================
!
! Routines:
!
! (1) mtxel_m() Originally By MLT Last Modified 6/5/2008 JRD
!
! input: crys, wfnc, wfnvq, gvec, eqp, types
! ik label of k-point in FBZ
!
! output: s0 matrix element of the momentum operator at point ik
!
! Calculates the momentum operator between two sets of wavefunctions
! < ic,k | P dot 2 (G+k+q) exp(-i q.r) | iv,k+q > / | P |
! Division by ( E_c^LDA - E_v^LDA ) is done only if divide_energy = .true.
! Each set has its own isort vector and the number of bands is nband
! The momentum operator is divided by electron mass m = 0.5 (in Ry atomic units)
! q is an optional small shift to k in reciprocal space
! P is the polarization vector
!
! (2) mtxel_v() Originally By MLT Last Modified: 6/5/2008 (JRD)
!
! input: wfnc, wfnvq, gvec types
! qshift length of the q-shift vector
!
! output: s0 velocity matrix elements at a k-point
!
! Calculates the velocity operator between two sets of wavefunctions
! < ic,k | exp(-i q.r) | iv,k+q > / q
! Note that this form is also correct for intraband transitions. --FHJ
! Each set has its own isort vector and the number of bands is nband
! q is a small but finite shift to k in reciprocal space
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
module mtxel_optical_m
  use global_m
  implicit none
  private
  public :: mtxel_m, mtxel_v
contains
subroutine mtxel_m(crys,wfnc,wfnvq,gvec,eqp,pol,s0_dim1,s0_dim2,s0,ik,divide_energy,kpt)
  type (crystal), intent(in) :: crys
  type (wavefunction), intent(in) :: wfnc, wfnvq
  type (gspace), intent(in) :: gvec
  type (eqpinfo), intent(in) :: eqp
  real(DP), intent(in) :: pol(3) !< light polarization for transition matrix elements
  integer, intent(in) :: s0_dim1, s0_dim2
  real(DP), intent(out) :: s0(:,:,:) !< (s0_dim1, s0_dim2, wfnc%nspin)
  integer, intent(in) :: ik
  logical, intent(in) :: divide_energy
  real(DP),intent(in),optional::kpt(3)
  real(DP) :: lpol, de
  real(DP) :: kpg(3), kpt_(3)
  integer :: ig, igq, ic, iv, isc, isp
  integer, allocatable :: isorti(:)
  real(DP) :: fac
  real(DP) :: sum
!---------------------------------
! Initialize isorti array
 
  s0 = 0.0d0
  allocate(isorti (gvec%ng))
  isorti(:)=0
  do ig=1, gvec%ng
    isorti(wfnvq%isort(ig)) = ig
  enddo
  kpt_(:) = 0d0
  if (present(kpt)) kpt_ = kpt
!----------------------------------
! Check if the polarization vector is properly defined
  lpol=sqrt(DOT_PRODUCT(pol,MATMUL(crys%bdot,pol)))
  if (abs(lpol).lt.TOL_Zero) then
    write(0,*) lpol, pol(:)
    call die("zero length polarization vector")
  endif
!----------------------------------
! Calculate s0(ic,iv) = < ic,k | P dot 2 (G+k+q) exp(-i q.r) | iv,k+q > / | P |
! / ( E_c^LDA - E_v^LDA )
! Here, q = 0 and (P dot 2 (G+k)) is replaced with (P dot 2 G)
! because < ic,k | P dot 2 k | iv, k > = P dot 2 k < ic,k | iv,k > = 0
! (only true for interband transitions. --DAS)
! (but we need k+G for nonlinear optics --DYQ)
  ! FHJ: This does not work with the PGI compiler. Don`t know if this is
  ! compiler bug yet.
  !!disabled PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) &
  !!disabled PRIVATE(isc,ic,iv,sum,ig,igq,kpg,fac,isp,de)
  do isc=1,wfnc%nspin
    do ic=1,s0_dim1
      do iv=1,s0_dim2
        sum=0.0d0
        !DIR$ LOOP COUNT MAX=2, MIN=1, AVG=1
        do isp=1,wfnc%nspinor
          !!OMP DO SIMD
          do ig=1, wfnc%ng
            igq = isorti(wfnc%isort(ig))
            kpg(1:3) = gvec%components(1:3,wfnc%isort(ig)) + kpt_(1:3)
            fac = DOT_PRODUCT(pol,MATMUL(crys%bdot,kpg))
            if (igq<=wfnvq%ng .and. igq>0) then
              sum = sum + (wfnc%cg(ig,ic,isc*isp)) * wfnvq%cg(igq,iv,isc*isp) * fac
            endif
          enddo
          !!OMP END DO SIMD
        enddo
        s0(ic,iv,isc) = 2.d0 * sum / lpol
        if(divide_energy) then
          de = (eqp%eclda(ic,ik,isc)-eqp%evlda(iv,ik,isc))
          if (abs(de) .lt. TOL_Degeneracy) then
            s0(ic,iv,isc) = 0.0d0
          else
            s0(ic,iv,isc) = s0(ic,iv,isc) / de
          end if
        endif
      enddo
    enddo
  enddo
  !!disabled END PARALLEL DO
  if(allocated(isorti))then;deallocate(isorti);endif
 
  return
end subroutine mtxel_m
!===============================================================================
subroutine mtxel_v(wfnc,wfnvq,gvec,qshift,s0_dim1,s0_dim2,s0)
  type (wavefunction), intent(in) :: wfnc, wfnvq
  type (gspace), intent(in) :: gvec
  real(DP), intent(in) :: qshift
  integer, intent(in) :: s0_dim1, s0_dim2
  real(DP), intent(out) :: s0(:,:,:) !< (s0_dim1, s0_dim2, wfnc%nspin)
  integer :: ig, igq, ic, iv, isc, isp
  integer, allocatable :: isorti(:)
  real(DP) :: sum
!--------------------------------
! Initialize isorti array
 
  s0 = 0.0d0
  allocate(isorti (gvec%ng))
  isorti(:)=0
  do ig=1, gvec%ng
    isorti(wfnvq%isort(ig)) = ig
  enddo
!--------------------------------
! Calculate s0(ic,iv) = < ic,k | exp(-i q.r) | iv,k+q > / q
  do isc=1,wfnc%nspin
    do ic=1,s0_dim1
      do iv=1,s0_dim2
        sum=0.0d0
        do ig=1, wfnc%ng
          igq=isorti(wfnc%isort(ig))
          if (igq.gt.wfnvq%ng) exit
          do isp=1,wfnc%nspinor
            sum = sum + (wfnc%cg(ig,ic,isc*isp)) * wfnvq%cg(igq,iv,isc*isp)
          enddo
        enddo
        s0(ic,iv,isc) = sum / qshift
      enddo ! iv
    enddo ! ic
  enddo ! isc
  if(allocated(isorti))then;deallocate(isorti);endif
 
  return
end subroutine mtxel_v
end module mtxel_optical_m
