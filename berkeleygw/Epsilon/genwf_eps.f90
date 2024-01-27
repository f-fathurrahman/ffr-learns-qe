!==========================================================================
!
! Module:
!
! genwf_eps Originally By FHJ Last Modified 04/2012 (FHJ)
!
!
!==========================================================================
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
module genwf_eps_m
  use global_m
  use input_utils_m
  use fftw_m
  use genwf_mpi_m
  use gmap_m
  use susymmetries_m
  use sort_m
  use timing_m, only: timing => epsilon_timing
  use genwf_mpi_m
  implicit none
  private
  !> communicator object for the WFN FFTs
  type wfn_FFT_comm_t
    logical :: done !< Have we received all the buffers?
    integer, pointer :: req_recvv(:), req_recvc(:) !< Array of MPI_REQUEST
    integer, pointer :: req_sendv(:), req_sendc(:) !< Array of MPI_REQUEST
    integer :: recv_cntv, recv_cntc !< Number of requests
    integer :: send_cntv, send_cntc !< Number of requests
    integer :: nv !< Number of valence bands
  end type wfn_FFT_comm_t
  public :: genwf_gen, genwf_FFT, free_wfns, genwf_lvl2, &
    genwf_FFT_Isend, genwf_FFT_Wait, wfn_FFT_comm_t, &
    get_wfn_fftgrid, get_eps_fftgrid
contains
  !> FHJ: Figure out what is the min/max gvec components give an isort array
  !! Note: you should manually initialize box_min and box_max to zero!
  subroutine get_gvecs_bounds(gvec, ng, isort, box_min, box_max)
    type(gspace), intent(in) :: gvec
    integer, intent(in) :: ng
    integer, intent(in) :: isort(:)
    integer, intent(inout) :: box_min(3), box_max(3)
    integer :: ig
   
    do ig=1,ng
      box_min(1:3) = min(box_min(1:3), gvec%components(1:3, isort(ig)))
      box_max(1:3) = max(box_max(1:3), gvec%components(1:3, isort(ig)))
    enddo
   
    return
  end subroutine get_gvecs_bounds
  !> FHJ: Figure out the minimum fftbox that holds all the WFNs
  subroutine get_wfn_fftgrid(pol, gvec, kp, intwfn)
    type(polarizability), intent(inout) :: pol
    type(gspace), intent(in) :: gvec
    type(kpoints), intent(in) :: kp
    type(int_wavefunction), intent(in) :: intwfn
    integer :: ik, wfn_box_min(3), wfn_box_max(3)
   
    wfn_box_min(:) = 0; wfn_box_max(:) = 0
    do ik=1,kp%nrk
      call get_gvecs_bounds(gvec, intwfn%ng(ik), intwfn%isort(:,ik), wfn_box_min, wfn_box_max)
    enddo
    pol%WFN_FFTgrid(1:3) = wfn_box_max(1:3) - wfn_box_min(1:3) + 1
    if (peinf%verb_debug .and. peinf%inode==0) then
      write(6,*) 'WFN min. FFT grid:',pol%WFN_FFTgrid
    endif
   
    return
  end subroutine get_wfn_fftgrid
  !> FHJ: Figure out the minimum fftbox that allows us to convolve the WFNs within
  !! an energy window of nmtx G vectors.
  subroutine get_eps_fftgrid(pol, gvec)
    type(polarizability), intent(inout) :: pol
    type(gspace), intent(in) :: gvec
    integer :: eps_box_min(3), eps_box_max(3)
   
    eps_box_min(:) = 0; eps_box_max(:) = 0
    call get_gvecs_bounds(gvec, pol%nmtx, pol%isrtx, eps_box_min, eps_box_max)
    ! FHJ: Note: the amount of padding is actually N_sig + N_window - 1
    pol%FFTgrid(1:3) = pol%WFN_FFTgrid(1:3) + (eps_box_max(1:3) - eps_box_min(1:3))
    pol%FFTgrid(1:3) = min(pol%FFTgrid(1:3), gvec%FFTgrid(1:3))
    if (peinf%verb_debug .and. peinf%inode==0) then
      write(6,'(1x,a,3(1x,i0))') 'Original FFT grid:', gvec%FFTgrid
      write(6,'(1x,a,3(1x,i0))') 'Minimal  FFT grid:', pol%FFTgrid
    endif
   
    return
  end subroutine get_eps_fftgrid
  !> FHJ: to be used internally with genwf_FFT_Isend
  subroutine do_my_FFTs(this,gvec,Nfft,wfn_fft,intwfn,my_bands,my_cnt,ng,tmp_wfn,isort,ind,ph,is_val)
    type (wfn_FFT_comm_t), intent(inout), target :: this !< communicator object for the WFN FFTs
    type(gspace), intent(in) :: gvec
    complex(DPC), intent(inout) :: wfn_fft(:,:,:,:)
    integer, intent(in) :: Nfft(3)
    type(int_wavefunction), intent(in) :: intwfn
    integer, intent(in) :: my_bands(:) !< my bands
    integer, intent(in) :: my_cnt !< number of FFTs to do = sizeof(my_bands)
    integer, intent(in) :: ng !< number of gvectors for the k-pt in question
    real(DP), intent(inout) :: tmp_wfn(:)!< buffer to reorder WFN using ind and ph
    integer, intent(in) :: isort(:) !< wfn isort
    integer, intent(in) :: ind(:)
    real(DP), intent(in) :: ph(:)
    logical, intent(in) :: is_val !< .true. to take the conjg_fftbox
    integer :: is, ig, fft_size
    integer :: ib_list, ib_local, ib, iproc, offset
    integer, pointer :: invindex(:), send_cnt, req_send(:)
    logical, pointer :: does_it_own(:,:)
   
    is = 1
    if ( is_val ) then
      invindex => peinf%invindexv
      does_it_own => peinf%does_it_ownv
      req_send => this%req_sendv
      send_cnt => this%send_cntv
      offset = 0
    else
      invindex => peinf%invindexc
      does_it_own => peinf%does_it_ownc
      req_send => this%req_sendc
      send_cnt => this%send_cntc
      offset = this%nv
    endif
    fft_size = product(Nfft(1:3))
    do ib_list = 1, my_cnt
      ib_local = my_bands(ib_list)
      do ig=1,ng
        tmp_wfn(ig) = intwfn%cg(ind(ig), ib_local, is)*ph(ig)
      enddo
      call timing%start(timing%opt_fft_fft)
      call put_into_fftbox(ng, tmp_wfn, gvec%components, &
        isort, wfn_fft(:,:,:,ib_local), Nfft)
      call do_FFT( wfn_fft(:,:,:,ib_local), Nfft, 1)
      if ( is_val ) call conjg_fftbox( wfn_fft(:,:,:,ib_local), Nfft )
      call timing%stop(timing%opt_fft_fft)
    enddo
   
    return
  end subroutine do_my_FFTs
  !> Generates all the real-space wavefunctions. Used only if pol%os_opt_fft==2
  !! This version avoids communication, and it`s under development!
  !!TODO`s:
  !! (1) we are just supporting one spin and one kpt/qpt.
  !! (2) support serial code
  subroutine genwf_FFT_Isend(this,crys,gvec,syms,kp,kpq,vwfn,pol,cwfn,intwfnv,intwfnvq,intwfnc)
    type (wfn_FFT_comm_t), intent(inout) :: this !< communicator object for the WFN FFTs
    type (crystal), intent(in) :: crys
    type (gspace), intent(in) :: gvec
    type (symmetry), intent(in) :: syms
    type (kpoints), target, intent(in) :: kp
    type (kpoints), target, intent(in) :: kpq
    type (valence_wfns), intent(inout) :: vwfn
    type (polarizability), intent(inout) :: pol
    type (conduction_wfns), intent(inout) :: cwfn
    type (int_wavefunction), intent(inout) :: intwfnv
    type (int_wavefunction), intent(inout) :: intwfnvq
    type (int_wavefunction), intent(inout) :: intwfnc
    integer :: npes_per_pool, nc_groups, nproc_max, nc
    integer, allocatable :: nv_bands(:), v_owners(:)
    integer :: Nfft(3), fft_size
    real(DP) :: scale
    integer, allocatable :: my_vbands(:)
    integer :: my_vcnt, iv, iv_local
    integer :: ipool, isubrank, inode, ioffset
    integer, allocatable :: grp_global_ranks(:), ntot_bands(:), &
      grp_local_owners(:)
    integer, allocatable :: my_cbands(:)
    integer :: my_ccnt, ic, ic_local
    integer :: my_grp_rank, grp_nprocs
    integer :: min_bands, iproc, iworker
    integer :: ik, is, ib
    !sort stuff
    real(DP), allocatable :: tmp_wfn(:)
    real(DP), allocatable :: ph(:)
    real(DP), allocatable :: ekin(:)
    integer, allocatable :: ind(:), isorti(:)
    integer, allocatable, target :: wfn_isort(:)
    integer :: ng0, ig
   
    call logit('generating all real-space wavefunctions')
    call timing%start(timing%opt_fft)
    if(pol%nq>1.or.kp%nrk>1.or.kpq%nrk>0.or.pol%need_WFNq) &
      call die('FFT opt. level 2 only works for 1 qpt and 1 kpt, and without WFNq',&
      only_root_writes=.true.)
    call logit('done generating real-space wavefunctions')
    call timing%stop(timing%opt_fft)
   
    return
  end subroutine genwf_FFT_Isend
  ! FHJ: call me after genwf_FFT_Isend, but just before you actually need the data
  subroutine genwf_FFT_Wait(this)
    type (wfn_FFT_comm_t), intent(inout) :: this !< communicator object for the WFN FFTs
   
   
    return
  end subroutine genwf_FFT_Wait
  subroutine genwf_lvl2(kp,kpq,vwfn,pol,cwfn)
    type (kpoints), target, intent(in) :: kp
    type (kpoints), target, intent(in) :: kpq
    type (valence_wfns), intent(inout) :: vwfn
    type (polarizability), intent(in) :: pol
    type (conduction_wfns), intent(inout) :: cwfn
    type(kpoints), pointer :: kp_point
   
    if(pol%need_WFNq) then ! FIXME I think this is wrong if pol%nq1>0
      kp_point => kpq
    else
      kp_point => kp
    endif
    allocate(vwfn%ev (vwfn%nband+pol%ncrit, kp%nspin))
    allocate(cwfn%ec (cwfn%nband,kp%nspin))
    vwfn%ev(1:vwfn%nband+pol%ncrit,1:kp%nspin) = &
      kp_point%el(1:vwfn%nband+pol%ncrit, vwfn%idx_kp, 1:kp%nspin)
    cwfn%ec(1:cwfn%nband,1:kp%nspin) = &
      kp%el(1:cwfn%nband, cwfn%idx_kp, 1:kp%nspin)
   
    return
  end subroutine genwf_lvl2
  !> Generates all the real-space wavefunctions. Used only if pol%os_opt_fft==2
  !!TODO`s:
  !! (1) we are just supporting one spin and one kpt/qpt.
  !! (2) communication can be reduced if we distribute the WFNs in a smarter way
  !! (3) support serial code
  subroutine genwf_FFT(crys,gvec,syms,kp,kpq,vwfn,pol,cwfn,intwfnv,intwfnvq,intwfnc)
    type (crystal), intent(in) :: crys
    type (gspace), intent(in) :: gvec
    type (symmetry), intent(in) :: syms
    type (kpoints), target, intent(in) :: kp
    type (kpoints), target, intent(in) :: kpq
    type (valence_wfns), intent(inout) :: vwfn
    type (polarizability), intent(inout) :: pol
    type (conduction_wfns), intent(inout) :: cwfn
    type (int_wavefunction), intent(inout) :: intwfnv
    type (int_wavefunction), intent(inout) :: intwfnvq
    type (int_wavefunction), intent(inout) :: intwfnc
    integer :: ik, is, ipe
    integer :: own_max
    integer, allocatable :: band_owners(:)
    integer :: ib_loc, ib, ng0
    integer :: receiver, recv_cnt, send_cnt
    integer :: local_band_idx
    integer, allocatable :: req_send(:), req_recv(:)
    integer, allocatable :: my_bands(:), isort0(:)
    real(DP), allocatable :: bufs_wfn(:,:)
    complex(DPC), allocatable :: work_ffts(:,:,:,:)
    integer :: Nfft(3), fft_size
    real(DP) :: scale
    real(DP), allocatable :: tmp_wfn(:)
    real(DP), allocatable :: ph(:)
    real(DP), allocatable :: ekin(:)
    integer, allocatable :: ind(:), isorti(:)
    integer, allocatable, target :: wfn_isort(:)
    integer :: ig
   
    call logit('generating all real-space wavefunctions')
    call timing%start(timing%opt_fft)
    if(pol%nq>1.or.kp%nrk>1.or.kpq%nrk>0.or.pol%need_WFNq) &
      call die('FFT opt. level 2 only works for 1 qpt and 1 kpt, and without WFNq',&
      only_root_writes=.true.)
   
    return
  contains
    logical function should_send(own_arr)
      logical, intent(in) :: own_arr(:,:)
      integer :: sender, ib_
     
      ib_ = ib
      if (ib_ > vwfn%nband) ib_ = ib_ - vwfn%nband
      should_send = .false.
     
      return
    end function should_send
  end subroutine genwf_FFT
  !> Generic routine that generates wavefunctions
  subroutine genwf_gen(syms,gvec,crys,kp,kpq,irk,rk,qq,vwfn,pol,cwfn,use_wfnq,intwfnv,intwfnvq,intwfnc,iv)
    type (symmetry), intent(in) :: syms
    type (gspace), intent(in) :: gvec
    type (crystal), intent(in) :: crys
    type (kpoints), target, intent(in) :: kp
    type (kpoints), target, intent(in) :: kpq
    integer, intent(in) :: irk
    real(DP), intent(in) :: rk(3)
    real(DP), intent(in) :: qq(3)
    type (valence_wfns), intent(inout) :: vwfn
    type (polarizability), intent(in) :: pol
    type (conduction_wfns), intent(inout) :: cwfn
    logical, intent(in) :: use_wfnq
    type (int_wavefunction), intent(in) :: intwfnv
    type (int_wavefunction), intent(in) :: intwfnvq
    type (int_wavefunction), intent(in) :: intwfnc
    integer, intent(in) :: iv
    integer :: ivr
   
    call logit('calling genwf')
    call timing%start(timing%genwf)
    if (iv .le. peinf%nvownactual) then
      call genwf_mpi(syms,gvec,crys,kp,kpq,irk,rk,qq,vwfn,pol,cwfn,use_wfnq,intwfnv,intwfnvq,intwfnc,iv)
    endif
    call timing%stop(timing%genwf)
    call logit('done genwf')
   
    return
  end subroutine genwf_gen
  !> Deallocate all "intermediate" wavefunctions
  subroutine free_wfns(pol, intwfnv, intwfnvq, intwfnc, free_all)
    type(polarizability), intent(in) :: pol
    type(int_wavefunction), intent(inout) :: intwfnv, intwfnvq, intwfnc
    logical, intent(in) :: free_all !< if .false., then only %ng and %isort will be preserved
   
    if (free_all) then
      if(associated(intwfnv%ng))then;deallocate(intwfnv%ng);nullify(intwfnv%ng);endif
      if(associated(intwfnv%isort))then;deallocate(intwfnv%isort);nullify(intwfnv%isort);endif
    endif
    if(associated(intwfnv%cg))then;deallocate(intwfnv%cg);nullify(intwfnv%cg);endif
    if(associated(intwfnv%qk))then;deallocate(intwfnv%qk);nullify(intwfnv%qk);endif
    if (pol%need_WFNq) then
      if (free_all) then
        if(associated(intwfnvq%ng))then;deallocate(intwfnvq%ng);nullify(intwfnvq%ng);endif
        if(associated(intwfnvq%isort))then;deallocate(intwfnvq%isort);nullify(intwfnvq%isort);endif
      endif
      if(associated(intwfnvq%cg))then;deallocate(intwfnvq%cg);nullify(intwfnvq%cg);endif
      if(associated(intwfnvq%qk))then;deallocate(intwfnvq%qk);nullify(intwfnvq%qk);endif
    endif
    if (free_all) then
      if(associated(intwfnc%ng))then;deallocate(intwfnc%ng);nullify(intwfnc%ng);endif
      if(associated(intwfnc%isort))then;deallocate(intwfnc%isort);nullify(intwfnc%isort);endif
    endif
    if(associated(intwfnc%cg))then;deallocate(intwfnc%cg);nullify(intwfnc%cg);endif
    if(associated(intwfnc%cbi))then;deallocate(intwfnc%cbi);nullify(intwfnc%cbi);endif
    if(associated(intwfnc%qk))then;deallocate(intwfnc%qk);nullify(intwfnc%qk);endif
   
    return
  end subroutine free_wfns
end module genwf_eps_m
