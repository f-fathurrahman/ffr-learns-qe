!============================================================================
!
! Routines:
!
! (1) mtxel_occ() Originally By ? Last Modified 8/15/2015 (FHJ)
!
! Adapted from subroutine mtxel. Subroutine computes required matrix
! elements that involve only occupied states, of the form <nk|exp{i g.r}|nk>
! They are required for COHSEX and exact CH calculations.
!
! input n,m band indices
! input ncoul number of matrix elements required
! input isrtrq index array for g-vectors in
! <nk|exp{i g.r}|nk>
! output aqs matrix elements required
!
!============================================================================
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
module mtxel_occ_m
  use fftw_m
  use global_m
  use misc_m
  implicit none
  private
  public :: &
    mtxel_occ
contains
subroutine mtxel_occ(n,m,gvec,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: n, m
  type (gspace), intent(in) :: gvec
  type (wfnkstates), intent(in) :: wfnk
  type (kpoints), intent(inout) :: kp
  integer, intent(in) :: ncoul
  integer, intent(in) :: isrtrq(gvec%ng)
  real(DP), intent(out) :: aqs(ncoul)
  integer, intent(in) :: ispin
!-------------------------
! If we are using FFT to calculate matrix elements...
! We use FFT to compute <u_nk|e^(iG.r)|u_nk> elements where
! u_nk is the periodic part of the wave function.
! The calculation is done in real space, and integration over
! the grid is replaced by the sum over the grid points p:
!
! <u_nk|e^(iG.r)|u_nk> =
! Volume/Np * sum_p { conj(u_nk(p))*e^(iG.p)*u_nk(p) }
!
! Since u_nk(p) = Volume^-0.5 * sum_G { cnk(G)*e^(iG.p) },
! and FFT is defined as FFT(cnk,+,p) = sum_G { cnk(G)*e^{+iG.p} },
! we must compute
!
! <u_nk|e^(iG.r)|u_nk>
! = 1/Np * sum_p { conj(FFT(cnk,+,p))*e^(iG.p)*FFT(cnk,+,p) }
! = 1/Np * FFT(conj(FFT(cnk,+,:)).*FFT(cnk,+,:),+,G)
!
! where .* is a point by point multiplication on the grid
  complex(DPC), dimension(:,:,:), allocatable :: fftbox1,fftbox2
  integer :: jsp
  integer, dimension(3) :: Nfft
  real(DP) :: scale
  real(DP), dimension(:), allocatable :: tmparray
 
! Compute size of FFT box we need and scale factor
  call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)
! Allocate FFT boxes
  allocate(fftbox1 (Nfft(1),Nfft(2),Nfft(3)))
  allocate(fftbox2 (Nfft(1),Nfft(2),Nfft(3)))
! Put the data for band n into FFT box 1 and do the FFT,zk(:,1)
  allocate(tmparray (ncoul))
  do jsp = ispin,ispin*kp%nspinor
    call put_into_fftbox(wfnk%nkpt,wfnk%zk((n-1)*wfnk%nkpt+1:,jsp),gvec%components,wfnk%isrtk,fftbox1,Nfft)
    call do_FFT(fftbox1,Nfft,1)
! We need the complex conjugate of the |nk> band actually
    call conjg_fftbox(fftbox1,Nfft)
! Now we get the matrix elements:
! Get n wave function and put it into box 2,
! do FFT,
! multiply by box1 contents,
! do FFT again,
! and extract the resulting matrix elements
    call put_into_fftbox(wfnk%nkpt,wfnk%zk((m-1)*wfnk%nkpt+1:,jsp),gvec%components,wfnk%isrtk,fftbox2,Nfft)
    call do_FFT(fftbox2,Nfft,1)
    call multiply_fftboxes(fftbox1,fftbox2,Nfft)
    call do_FFT(fftbox2,Nfft,1)
    call get_from_fftbox(ncoul,tmparray,gvec%components,isrtrq,fftbox2,Nfft,scale)
    if (kp%nspinor.eq.1 .or. jsp.eq. 1) then
      aqs(:) = tmparray(:)
    else
      aqs(:) = aqs(:) + tmparray(:)
    endif
  enddo
  if(allocated(tmparray))then;deallocate(tmparray);endif
! We are done, so deallocate FFT boxes
  if(allocated(fftbox1))then;deallocate(fftbox1);endif
  if(allocated(fftbox2))then;deallocate(fftbox2);endif
 
  return
end subroutine mtxel_occ
end module mtxel_occ_m
