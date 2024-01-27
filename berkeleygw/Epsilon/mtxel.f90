!===========================================================================
!
! Routines()
!
! (1) mtxel() Originally By (?) Last Modified 5/6/2012 (FHJ)
!
! Compute matrix elements (gme) for valence state iv with all
! conduction bands and for all G-vectors.
!
! <c,k,ispin|exp{-i(q+G).r}|v,k+q,ispin> = [M_{vc}(k,q,G)]^*
!
! On exit,
! pol%eden(band,spin) = 1/(e_val-e_cond) = energy denominators
! pol%gme(band,g-vector,spin) = plane wave matrix elements
! pol%isrtx orders the |G(i)|^2 i=1,pol%nmtx
! vwfn%isort orders |qk+g|^2 (in vwfn type)
!
! energies are apparently assumed in Rydbergs.
!
!===========================================================================
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
module mtxel_m
  use global_m
  use fftw_m
  use misc_m
  use lin_denominator_m
  use timing_m, only: timing => epsilon_timing
  implicit none
  private
  public :: mtxel_init_FFT_cond, mtxel_free_FFT_cond, mtxel
contains
!> Precalculates the FFTs for all conduction bands
subroutine mtxel_init_FFT_cond(gvec,pol,cwfn,kp)
  type (gspace), intent(in) :: gvec
  type (polarizability), intent(inout) :: pol
  type (conduction_wfns), intent(inout) :: cwfn
  type (kpoints), intent(inout) :: kp
  real(DP) :: scale
  integer, dimension(3) :: Nfft
  integer :: j, offset_g, jsp, ict
 
  if (kp%nspinor>1) then
      call die('Epsilon on Steroids only supports one spin at the moment.')
  endif
  call timing%start(timing%opt_fft)
  call timing%start(timing%opt_fft_fft)
  call setup_FFT_sizes(pol%FFTgrid,Nfft,scale)
  allocate(cwfn%wfn_fft (Nfft(1),Nfft(2),Nfft(3),peinf%ncownactual))
  jsp = 1
! JRD XXX We should we be doing a many_fft call here or what is the point of allocating a bigger array
  do j=1,peinf%ncownactual
    offset_g = (j-1)*cwfn%ngc
    call put_into_fftbox(cwfn%ngc,cwfn%zc(offset_g+1:,jsp),gvec%components,cwfn%isort,cwfn%wfn_fft(:,:,:,j),Nfft)
    call do_FFT(cwfn%wfn_fft(:,:,:,j),Nfft,1)
  enddo
  call timing%stop(timing%opt_fft)
  call timing%stop(timing%opt_fft_fft)
 
  return
end subroutine mtxel_init_FFT_cond
!> Frees ffts_cond buffer
subroutine mtxel_free_FFT_cond(cwfn)
  type (conduction_wfns), intent(inout) :: cwfn
 
  if(associated(cwfn%wfn_fft))then;deallocate(cwfn%wfn_fft);nullify(cwfn%wfn_fft);endif
 
  return
end subroutine mtxel_free_FFT_cond
subroutine mtxel(iv,gvec,vwfn,cwfn,pol,ispin,irk,kp,kpq,rank_mtxel,kfact)
  integer, intent(in) :: iv
  type (gspace), intent(in) :: gvec
  type (valence_wfns), intent(in) :: vwfn
  type (conduction_wfns), intent(in) :: cwfn
  type (polarizability), intent(inout) :: pol
  type (kpoints), intent(inout) :: kp,kpq
  integer, intent(in) :: ispin,irk
  integer, intent(in) :: rank_mtxel
  real(DP), intent(in) :: kfact
  integer :: ic_loc, ic, ic_FE, offset_g, jsp, ig, iv_loc
  integer :: freq_idx
  integer :: ia, ib, ict
  integer, dimension(3) :: Nfft
  real(DP) :: scale, eden
  complex(DPC), dimension(:,:,:), allocatable :: fftbox1,fftbox2
  real(DP), dimension(:), allocatable :: tmparray
  logical :: keep_transition
 
  call timing%start(timing%mtxel_fft)
  ! FHJ: notation:
  ! iv_loc -> local valence band index
  ! iv -> global valence band index (iv==1 := lowest-energy band)
  ! ic_loc -> local conduction band index
  ! ic -> global conduction band (ic==nv+1 := lowest-energy unocc. band)
  ! ic_FE -> global conduction band, but starting at FE
  ! (ic_FE==1 := first cond band). Rarely used.
  iv_loc = peinf%indexv(iv)
  if(pol%nfreq_group .gt. 1 .and. pol%gcomm .eq. 0) then
    freq_idx = 1
  else
    freq_idx = rank_mtxel+1
  endif
!--------------------------
! Use FFTs to calculate matrix elements
! Compute size of FFT box we need
  call setup_FFT_sizes(pol%FFTgrid,Nfft,scale)
! Allocate FFT boxes
  allocate(fftbox2 (Nfft(1),Nfft(2),Nfft(3)))
! Put the data for valence band iv into FFT box 1 and do the FFT
  if (pol%os_opt_ffts/=2) then
    allocate(fftbox1 (Nfft(1),Nfft(2),Nfft(3)))
  endif
! JRD XXX This needs to be threaded
!disabled PARALLEL DO collapse(2)
  do ic_loc = 1, peinf%ncownactual
    do ig = 1, pol%nmtx
      pol%gme(ig, ic_loc, iv_loc, ispin, irk, freq_idx) = 0.0d0
    enddo
  enddo
  allocate(tmparray (pol%nmtx))
  do jsp=ispin,ispin*kp%nspinor
    if (pol%os_opt_ffts/=2) then
      call put_into_fftbox(vwfn%ngv,vwfn%zv(:,jsp),gvec%components,vwfn%isort,fftbox1,Nfft)
      call do_FFT(fftbox1,Nfft,1)
      ! We need the complex conjugate of u_{vk+q)(r) for the cross correlation
      call conjg_fftbox(fftbox1,Nfft)
    endif
! Now we loop over the conduction states and get the matrix elements:
! 1. Get conduction wave function and put it into box 2,
! 2. do FFT, get u_{ck}(r)
! 3. multiply by box1 contents, get F(r) = [u_{vk+q)(r)]^* u_{ck}(r)
! 4. do FFT again, and extract the resulting matrix elements and put the into pol
! We conjugate the final result since we really want <ck|e^{-i(q+G).r}|vk+q>
! but we have calculated <vk+q|e^{i(q+G).r}|ck>.
    do ic_loc = 1, peinf%ncownactual
      ic_FE = peinf%invindexc(ic_loc)
      ic = vwfn%nband + ic_FE
      offset_g = (ic_loc-1)*cwfn%ngc
      if (pol%os_opt_ffts==2) then
        ! FHJ: optimization level 2 precomputed all the FFTs
        call timing%start(timing%fft_put)
!disabled PARALLEL DO collapse(2) PRIVATE(ia,ib,ict)
        do ia = 1, Nfft(3)
        do ib = 1, Nfft(2)
          do ict = 1, Nfft(1)
            fftbox2(ict,ib,ia) = vwfn%wfn_fft(ict,ib,ia,iv_loc) * cwfn%wfn_fft(ict,ib,ia,ic_loc)
          enddo
        enddo
        enddo
!disabled END PARALLEL DO
        call timing%stop(timing%fft_put)
     else
        if (pol%os_opt_ffts==1) then
          call timing%start(timing%fft_put)
          ! FHJ: Optimization level 1 precalculated at least these cond. FFTs
!disabled PARALLEL DO collapse(2) PRIVATE (ia,ib,ict)
          do ia = 1, Nfft(3)
          do ib = 1, Nfft(2)
            do ict = 1, Nfft(1)
              fftbox2(ict,ib,ia) = fftbox1(ict,ib,ia) * cwfn%wfn_fft(ict,ib,ia,ic_loc)
            enddo
          enddo
          enddo
!disabled END PARALLEL DO
          call timing%stop(timing%fft_put)
        else
          call put_into_fftbox(cwfn%ngc,cwfn%zc(offset_g+1:,jsp),gvec%components,cwfn%isort,fftbox2,Nfft)
          call do_FFT(fftbox2,Nfft,1)
          call multiply_fftboxes(fftbox1,fftbox2,Nfft)
        endif
      endif
      call do_FFT(fftbox2,Nfft,1)
      call get_from_fftbox(pol%nmtx,tmparray,gvec%components,pol%isrtx,fftbox2,Nfft,scale)
      keep_transition = .true.
      call timing%start(timing%fft_mltply)
      if (keep_transition) then
!disabled PARALLEL DO
        do ia = 1, pol%nmtx
          pol%gme(ia, ic_loc, iv_loc, ispin, irk, freq_idx) = &
            pol%gme(ia, ic_loc, iv_loc, ispin, irk, freq_idx) + (tmparray(ia))*kfact
        enddo
!disabled END PARALLEL DO
      endif
      ! Get energy denominator (static), or transition energy (FF)
      call get_eden()
      if (kp%nspinor.eq.1 .or. jsp.eq.2) then
        if (pol%freq_dep .eq. 0) then
!disabled PARALLEL DO
          do ia = 1, pol%nmtx
            pol%gme(ia, ic_loc, iv_loc, ispin, irk, 1) = &
              pol%gme(ia, ic_loc, iv_loc, ispin, irk, 1) * &
              sqrt(-1D0*eden)
          enddo
!disabled END PARALLEL DO
        endif
      endif
      call timing%stop(timing%fft_mltply)
    enddo
  enddo
  if(allocated(tmparray))then;deallocate(tmparray);endif
! We are done, so deallocate FFT boxes
  if (pol%os_opt_ffts.ne.2) then
    if(allocated(fftbox1))then;deallocate(fftbox1);endif
  endif
  if(allocated(fftbox2))then;deallocate(fftbox2);endif
! End FFT Case
!---------------------------
  call timing%stop(timing%mtxel_fft)
 
  return
contains
  !> FHJ: Get energy denominator (static), or transition energy (FF)
  subroutine get_eden()
    type(cvpair_info) :: lin_edenTemp
    real(DP) :: eval, econd, occ_v, occ_c, occ_diff
    real(DP) :: vk(2), vkq(2)
   
    call timing%start(timing%mtxel_denom)
    ! FHJ: See convention for conduction/valence band indices above.
    eval = vwfn%ev(iv,ispin)
    econd = cwfn%ec(ic, ispin)
    eden = 0d0
    ! guess occupations based on efermi; eventually this should be replaced by use of kp%occ
    if(eval*ryd > pol%efermi + TOL_Degeneracy) then
      occ_v = 0d0
    else if (eval*ryd < pol%efermi - TOL_Degeneracy) then
      occ_v = 1d0
    else
      occ_v = 0.5 ! within TOL_Degeneracy of the Fermi level, use FD(E_F) = 1/2
    endif
    if(econd*ryd > pol%efermi + TOL_Degeneracy) then
      occ_c = 0d0
    else if (econd*ryd < pol%efermi - TOL_Degeneracy) then
      occ_c = 1d0
    else
      occ_c = 0.5 ! within TOL_Degeneracy of the Fermi level, use FD(E_F) = 1/2
    endif
    occ_diff = occ_v - occ_c
    ! FHJ: Note that eden means different things depending on pol%freq_dep
    ! static: eden := 1/(transition energy).
    ! FF: eden := (transition energy). (I know...)
    ! In the static case, we lump sqrt(eden) into the matrix elements gme.
    ! In the FF case, we have to put it by hand for each frequency we evaluate chi.
    if (pol%freq_dep==0) then
      if(eval - econd < TOL_Degeneracy .and. occ_diff > TOL_Zero) then
        ! avoid dividing by zero or making eden > 0
        eden = occ_diff / (eval - econd)
      else
        eden = 0d0 ! in this case, occ_diff = 0 too
      endif
    elseif (pol%freq_dep==2 .or. pol%freq_dep==3) then
      ! FHJ: In chi_summation, we explicitly neglect transitions if eden<TOL_ZERO.
      ! That`s the way we keep track of forbidden transitions.
      if(eval - econd < TOL_Degeneracy .and. occ_diff > TOL_Zero) then
        eden = (eval - econd) / occ_diff
      else
        eden = 0.0d0
      endif
      pol%edenDyn(iv_loc, ic_loc, ispin, irk, freq_idx) = eden
    endif !pol%freq_dep
    call timing%stop(timing%mtxel_denom)
   
  end subroutine get_eden
end subroutine mtxel
end module mtxel_m
