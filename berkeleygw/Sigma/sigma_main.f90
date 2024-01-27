!===============================================================================
!
! Routines:
!
! (1) sigma Originally By MSH Last Modified 10/5/2009 (gsm)
!
! This is the main routine for the Sigma code. Please see the documentation
! in the README for more information on this code.
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
program sigma
  use acc_mtxel_kernels_m, only : acc_mtxel_sig, &
                                  allocate_acc_mtxel_sig, &
                                  deallocate_acc_mtxel_sig
  use algos_sigma_m
  use ch_converge_m
  use check_screening_m
  use checkbz_m
  use checkgriduniformity_m
  use epscopy_m
  use epsread_hdf5_m
  use fftw_m
  use fixwings_m
  use fullbz_m
  use genwf_mpi_m
  use global_m
  use gmap_m
  use input_m
  use input_outer_m
  use input_utils_m
  use io_utils_m
  use irrbz_m
  use misc_m
  use mtxel_cor_m
  use mtxel_m
  use mtxel_occ_m
  use mtxel_vxc_m
  use references_m
  use shiftenergy_dyn_m
  use shiftenergy_m
  use sort_m
  use subgrp_m
  use timing_m, only: common_timing, timing => sigma_timing
  use vcoul_generator_m
  use wfn_rho_vxc_io_m
  use write_program_header_m
  use write_result_dyn_hp_m
  use write_result_dyn_m
  use write_result_hp_m
  use write_result_m
  implicit none
!---------------------
! Derived Types
  type (crystal) :: crys
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (kpoints) :: kp
  type (siginfo) :: sig
  type (wpgen) :: wpg
  type (wfnkstates) :: wfnk,wfnkoff
  type (wfnkqstates) :: wfnkq
  type (epsmpiinfo) :: epsmpi
  type (wfnkqmpiinfo) :: wfnkqmpi
  type (wfnkmpiinfo) :: wfnkmpi
  type (twork_scell) :: work_scell
  type(progress_info) :: prog_info !< a user-friendly progress report
!------ Electron Phonon (EP) ----------
! ZL: allocate same variables to read wave functions
! while without overwriting the useful information
! phonq/phonv: electron phonon q/v (v: nu, phonon mode band index)
! for perturbed wave functions: dqv wfn(k)
  type (crystal) :: ep_crys
  type (symmetry) :: ep_syms
  type (gspace) :: ep_gvec
  type (kpoints) :: ep_kp
  ! type (siginfo) :: ep_sig ! sig will be re-used because somehow inread
                                 ! does not allow to be called twice
  type (wpgen) :: ep_wpg
  type (wfnkstates) :: ep_dq_wfnk
  type (wfnkqstates) :: ep_dq_wfnkq
  type (wfnkqmpiinfo) :: ep_dq_wfnkqmpi
  type (wfnkmpiinfo) :: ep_dq_wfnkmpi
  ! extra variables
  type (wfnkstates) :: wfnk_phonq, wfnk_phonq_off ! wfn(k+phonq), regular wfn at a different point
                                                  ! wfnk_phonq plyas the role of wfnk
                                                  ! wfnk_phonq_off, for dVXC mtxel
  type (wfnkmpiinfo) :: wfnk_phonq_mpi ! distributed
  type (wfnkstates) :: ep_dq_wfnk_tmp ! temporary variable
  type (wfnkmpiinfo) :: ep_dq_wfnkmpi_tmp ! temporary variable
  ! define variables associated with the second term in formalism
  type (wfnkqstates) :: wfnkq_phonq ! inner wfn(k-q+phonq)
  type (wfnkqstates) :: ep_dmq_wfnkq_phonq ! inner d_{-phonq,v} wfn(k-q+phonq), need time-reversal
                                       ! "mq" represents "minus q": "-q"
!---------------------
! k-points for the sum over BZ (the rq-points)
  integer :: nm, nrq, iout, iparallel, igp, igp_loc, igp_owner, nq0
  integer, allocatable :: neq(:),indrq(:),itnrq(:),kg0(:,:)
  real(DP), allocatable :: rq(:,:)
  type(grid) :: gr
!---------------------
! (k-q) kpoints involved in sigma summations
! ===> see data type cwfnkq and vwfnkq
  real(DP) :: rkq(3)
  real(DP) :: rkq_ep(3), rk_phonq(3), r_mkq_mphonq(3) ! used for EP
  integer :: ib_trs_ep, is_trs_ep ! looper for bands and spins, in the time-reversal case
!---------------------
! Dielectric matrices
! ZL: for EP under the static screening approximation, use the same eps
  integer :: ngq, neps, nmtx, ncoul, ncoulch, ncoulb, ncouls, ngpown_q
  integer :: ngqt,nmtxt,nfreqgpp
  integer, allocatable :: isrtq(:),isrtqi(:)
  integer, allocatable :: isrtrq(:)
  integer, pointer :: isrtrqi(:)
  real(DP), pointer :: eps(:,:)
  real(DP), allocatable :: vcoul(:), ekin(:)
  complex(DPC), pointer :: epsR(:,:,:),epsA(:,:,:)
!---------------------
! Matrix elements for sigma
  real(DP), allocatable :: aqs(:,:), aqsaug(:,:,:,:), alda(:,:), alda2(:,:), &
    ax(:,:), asx(:,:,:), ach(:,:,:), asig(:,:), ach_n1(:,:,:), achcor_n1(:,:,:), akih(:,:)
  ! ZL: alda is for VXC, and akih(:,:) is for KIH: Kinetic+Ionic+Hartree
  real(DP), pointer :: aqsch(:), aqsaugchd(:,:,:), aqsaugcho(:,:,:), acht_n1(:), achtcor_n1(:)
  real(DP), allocatable :: enew(:,:),efsto(:,:),zrenorm(:,:)
  complex(DPC), allocatable :: achcor(:,:)
  complex(DPC), allocatable :: asig_imag(:,:)
  complex(DPC), allocatable :: asxDyn(:,:,:), achDyn(:,:,:), achDyn_cor(:,:,:), &
    achDyn_corb(:,:,:), ach2Dyn(:,:,:), asigDyn(:,:), achD_n1(:,:,:)
  complex (DPC), pointer :: achtD_n1(:)
  complex(DPC), allocatable :: efstoDyn(:,:), enewDyn(:,:), enewDyn_nosr(:,:)
  integer, allocatable :: neqp1(:,:), neqp1_nosr(:,:)
!----------------------
! eps distrib variables
  real(DP), allocatable :: epstemp(:)
  character :: tmpstr*120
  character :: tmpfn*16
  character*20 :: fnc,fnk,fne
  character*16 :: routnam(100)
  integer :: routsrt(59),nullvec(3)
  logical :: xflag,imagvxcflag,imagxflag,found,q0flag,bExactlyZero, imagkihflag ! ZL: last one for KIH
  logical :: eqp1_warns(4) ! (GPP extrap, FF extrap, FF multiple solns, FF no soln)
  integer :: ig,i,j,k,itran,ikn,ika,ioff,error
  integer :: in,im,iw,ib,jb,idum,kg(3),jj,ii,ispin,jsp,g1,g2
  integer :: ncount,ndum,nbandi,tag,dest,source,nfold,ifold
  integer :: iwlda,irq,irq_,irq_min,n1,ierr
  integer :: s2,iunit_c,iunit_k,iunit_eps,ndv_ikn,iunit
  integer, allocatable :: ind(:), indinv(:)
  real(DP) :: fact,coulfact,weight,tempval,occ
  real(DP) :: qshift(3),oneoverq,qlen,q0len,vq(3),qk(3)
  real(DP) :: tsec(2),diffmin,diff,e_lk,avgcut,subcut,freq0
  complex(DPC), pointer :: achtDyn(:),achtDyn_cor(:),asxtDyn(:),ach2tDyn(:),achtDyn_corb(:)
  real(DP) :: achtcor,axt,epshead,asigt_imag
  real(DP), pointer :: asxt(:), acht(:)
  real(DP), allocatable :: ph(:)
  logical :: skip_checkbz, is_subq
!------------- ZL: Electron phonon variables --------------------------------
  logical :: ep_read ! controls reading dWFN of EP
  logical :: ep_debug ! controls debug of EP
  integer :: ik_phonq_idx ! similar role of ikn
  real(DP) :: k_phonq_coord(3) ! similar role of qk(3)
  logical :: check_norms_save
  integer :: ioff_check
  real(DP) :: Eo_save
! sigma-subspace variables --------------------------------------------------
  integer :: ipe_wing, my_pos_loc_wing, iproc_dum
  integer :: nfreq_fixwings, ifreq
  complex(DPC), pointer :: epsR_corrections(:,:,:)
  integer, dimension(3) :: Nfft
  real(DP) :: dummy_scale
  ep_read = .false.
  ep_debug = .false.
!--------------- Begin Program -------------------------------------------------
  call peinfo_init()
!----------------------
! Initialize random numbers
  peinf%jobtypeeval = 1
!------------------------
! Initialize timer
  call timing%init()
  call common_timing%init()
  call timing%start(timing%total)
!------------------------
! Initialize files
  call open_file(55,file='sigma.inp',form='formatted',status='old')
  if(peinf%inode == 0) then
    call open_file(8,file='sigma_hp.log',form='formatted',status='replace')
    call open_file(30,file='eqp0.dat',form='formatted',status='replace')
    call open_file(31,file='eqp1.dat',form='formatted',status='replace')
  endif
  call write_program_header('Sigma', .false.)
!------- Read crys data and wavefunctions from WFN_inner ----------------------------
! JRD: Included in input is the inread routine which reads the
! job parameters from sigma.inp and initializes the XC potential
  call timing%start(timing%input)
  call input(crys,gvec,syms,kp,wpg,sig,wfnk,iunit_c,iunit_k,fnc,fnk,wfnkqmpi,wfnkmpi,wfnk_phonq,wfnk_phonq_mpi,ep_read_in=.false.)
  call output_algos()
!------------------------
! Initialize GPU
  if (sigma_gpp_algo == OPENACC_ALGO .or. mtxel_algo == OPENACC_ALGO) then
    call initialize_gpu(OPENACC_ALGO)
  end if
  if (sigma_gpp_algo == OMP_TARGET_ALGO .or. mtxel_algo == OMP_TARGET_ALGO) then
    call initialize_gpu(OMP_TARGET_ALGO)
  end if
  call timing%stop(timing%input)
  call timing%start(timing%input_outer)
  ! ZL: for EP, we do not use WFN_outer. Everything should be in WFN_inner and dWFN
  if (sig%elph) then
  else
    call input_outer(crys,gvec,syms,kp,sig,wfnk,iunit_k,fnk,wfnkmpi)
  endif
  call timing%stop(timing%input_outer)
  if(associated(sig%kpt))then;deallocate(sig%kpt);nullify(sig%kpt);endif
  if(allocated(kp%ifmin))then;deallocate(kp%ifmin);endif
  if(allocated(kp%ifmax))then;deallocate(kp%ifmax);endif
!-------------------
! Initialize Various Parameters from inread
  eqp1_warns(:)=.false.
! ZL: KIH not supported with Real flavor
  if (sig%use_kihdat) then
    call die("KIH usage for arbitrary functional starting point MUST BE complied with COMPLEX version!")
  endif
! imaginary parts of diagonal vxc or exchange matrix elements
  imagvxcflag = .false.
  imagxflag = .false.
  imagkihflag = .true. ! ZL: always true for kih
! fraction of bare exchange
  if (abs(sig%xfrac).GT.TOL_Small) then
    xflag=.true.
  else
    xflag=.false.
  endif
!---------------------
! Open Various Files
  if (peinf%inode .eq. 0) then
    if (.not.(sig%freq_dep .eq. 0 .and. sig%exact_ch .eq. 1) .and. .not. (sig%freq_dep == -1)) then
      call open_file(127,file='ch_converge.dat',form='formatted',status='replace')
    endif
    if (sig%iwritecoul .eq. 1) then
      call open_file(19,file='vcoul',form='formatted',status='replace')
    endif
   ! This if for the hybrid functional calculations (one shot) otherwise just open x.dat
    if (sig%coul_mod_flag .and. (.not. sig%use_vxc2dat)) then
      call open_file(121,file='vxc2.dat',form='formatted',status='replace')
    else if ((.not. sig%use_xdat) .and. xflag .and. (.not. sig%coul_mod_flag)) then
      call open_file(119,file='x.dat',form='formatted',status='replace')
    endif
    if (.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX) then
      if(.not.sig%use_kihdat) then ! ZL: only meaningful to generate vxc.dat if not running with KIH
        call open_file(120,file='vxc.dat',form='formatted',status='replace')
      endif
    endif
  endif
!---------------------
! Write header of sigma_hp.log file
  if (peinf%inode.eq.0) then
    write(8,601) sig%freq_dep
    write(8,602) sig%bmin,sig%bmax
    write(8,603) sig%loff,sig%toff
    write(8,604) sig%fdf
    if(sig%fdf /= -2) write(8,605) sig%dw
    write(8,606) syms%ntran
    do itran=1,syms%ntran
      write(8,607) itran, ((syms%mtrx(i,j,itran),i=1,3),j=1,3)
    enddo
    write(8,*)
  endif
601 format(/,1x,"frequency_dependence",i4)
602 format(/,1x,"band_index",2i6)
603 format(1x,"sigma_matrix",i6,i4)
604 format(/,1x,"finite_difference_form",i4)
605 format(1x,"finite_difference_spacing",f10.6)
606 format(/,1x,"symmetries",/,1x,"ntran  =",i3)
607 format(1x,"mtrx",i2.2,1x,"=",9i3)
!---------------------
! JRD: Initialize the Full Frequency output files
  if (peinf%inode.eq.0 .and. (sig%freq_dep.eq.2 .or. (sig%fdf.eq.-3 .and. sig%freq_dep.eq.1))) then
    call open_file(8000,file='spectrum.dat',form='formatted',status='replace')
  endif
!---------------------
! Determine nq and neps
! JRD: This performs significantly better with hdf5
  call timing%start(timing%read_neps)
  call epscopy_init(sig, neps)
  if (sig%freq_dep/=-1) then
    ! FHJ: sig%nq and sig%qpt already defined if this is a HF calculation
    sig%nq = sig%nq0 + sig%nq1
    allocate(sig%qpt (3,sig%nq))
  endif
  if (sig%nq0==0) call die('There is no q->0 point in your calculation!', only_root_writes=.true.)
  epsmpi%nb = 1
  epsmpi%ngpown = NUMROC(neps, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
  epsmpi%ngpown_max = NUMROC(neps, epsmpi%nb, 0, 0, peinf%npes_pool)
  ! define distribution for subspace matrices (MDB)
  if(sig%do_sigma_subspace) then
    ! this has to 1 (no change)
    sig%epssub%nb_sub = 1
    !
    sig%epssub%Nbas_own_max = NUMROC(sig%neig_sub_max, sig%epssub%nb_sub, 0, 0, peinf%npes_pool)
    sig%epssub%Nbas_own = MIN(sig%epssub%Nbas_own_max*peinf%pool_rank + 1, sig%neig_sub_max)
    sig%epssub%Nbas_own = MIN(sig%epssub%Nbas_own_max*(peinf%pool_rank+1), sig%neig_sub_max) - sig%epssub%Nbas_own + 1
    ! note the block size here is the same as the epsmpi (not a good idea, this will be removed in the future)
    sig%epssub%ngpown_sub_max = NUMROC(neps, epsmpi%nb, 0, 0, peinf%npes_pool)
    sig%epssub%ngpown_sub = MIN(sig%epssub%ngpown_sub_max*peinf%pool_rank + 1, neps)
    sig%epssub%ngpown_sub = MIN(sig%epssub%ngpown_sub_max*(peinf%pool_rank+1), neps) - sig%epssub%ngpown_sub + 1
    sig%epssub%neps = neps
    !
    !XXX write(*,*) peinf%inode, peinf%pool_rank, sig%epssub%Nbas_own, sig%epssub%Nbas_own_max,&
    !XXX sig%epssub%ngpown_sub, sig%epssub%ngpown_sub_max
  end if
!----------------------------
! Allocate arrays
! ZL: allocate wfnkq contents
! wfn(k-q) where q is the internal regular q (or p in the equation)
  allocate(wfnkq%isrtkq (gvec%ng))
  allocate(wfnkq%ekq (sig%ntband,kp%nspin))
  allocate(alda (sig%ndiag+sig%noffdiag,sig%nspin))
  ! ZL: add for array KIH
  allocate(akih (sig%ndiag+sig%noffdiag,sig%nspin))
! ZL: 'a' means ARRAY, for bare exchange
  allocate(ax (sig%ndiag+sig%noffdiag,sig%nspin))
! achcor for static remainder
  allocate(achcor (sig%ndiag+sig%noffdiag,sig%nspin))
  allocate(asig_imag (sig%ndiag+sig%noffdiag,sig%nspin))
  allocate(achcor_n1 (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))
  if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
      ! static remainder for 1 band, allocated for all bands
    allocate(achtcor_n1 (sig%ntband))
  endif
  if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
    if (sig%fdf .eq. -3) then ! fdf: finite difference form, way to eval quasi energy, -3 many freqs
      ! ZL: we do not consider unusual cases for EP now
      allocate(asx (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
      allocate(ach (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
      if (sig%elph) call die('Electron phonon not set up for fdf = -3 yet')
    else
      ! at most use 3 freqs
      allocate(asx (3,sig%ndiag+sig%noffdiag,sig%nspin))
      allocate(ach (3,sig%ndiag+sig%noffdiag,sig%nspin))
    endif
    ! total self energy for final results (maybe)
    allocate(asig (sig%ndiag+sig%noffdiag,sig%nspin))
    ! t means temperory
    allocate(acht_n1 (sig%ntband))
    ! COH resolved for all bands, convergence usage
    allocate(ach_n1 (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))
    ! enew eqp1
    allocate(enew (sig%ndiag,sig%nspin))
    ! efsto eqp0
    allocate(efsto (sig%ndiag,sig%nspin))
    ! Z renorm factor
    allocate(zrenorm (sig%ndiag,sig%nspin))
    if (sig%fdf.eq.-3) then
      nfreqgpp=sig%nfreqeval
      allocate(asxt (sig%nfreqeval))
      allocate(acht (sig%nfreqeval))
      ! ZL: if this is the case for EP, it should have happened already
      if(sig%elph) call die('Electron phonon not set up for fdf = -3 yet')
    else
      nfreqgpp=3
      allocate(asxt (3))
      allocate(acht (3))
    endif
  endif
  ! full frequency calculation, labeled by Dyn
  if (sig%freq_dep.eq.2) then
    allocate(asxDyn (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
    allocate(achDyn (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
    allocate(achDyn_cor (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
    allocate(achDyn_corb (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
    allocate(ach2Dyn (sig%nfreqeval,sig%ndiag+sig%noffdiag,sig%nspin))
    allocate(asigDyn (sig%ndiag+sig%noffdiag,sig%nspin))
    allocate(achtD_n1 (sig%ntband))
    allocate(achD_n1 (sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin))
    allocate(efstoDyn (sig%ndiag,sig%nspin))
    allocate(enewDyn (sig%ndiag,sig%nspin))
    allocate(enewDyn_nosr (sig%ndiag,sig%nspin))
    allocate(neqp1 (sig%ndiag,sig%nspin))
    allocate(neqp1_nosr (sig%ndiag,sig%nspin))
    allocate(asxtDyn (sig%nfreqeval))
    allocate(achtDyn (sig%nfreqeval))
    allocate(achtDyn_cor (sig%nfreqeval))
    allocate(achtDyn_corb (sig%nfreqeval))
    allocate(ach2tDyn (sig%nfreqeval))
    if(sig%elph) call die('Electron phonon not set up for full frequency calculation')
  endif
  ! ZL: isrtrq is of the length of all gvec pool
  allocate(isrtrq (gvec%ng))
  ! GPP is freq_dep = 1, the only case we consider for EP for now
  if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then
    allocate(isrtrqi (gvec%ng))
    if (sig%elph) call die('Electron phonon not set up for sig%freq_dep.eq.0.or.sig%exact_ch.eq.1')
  endif
  ! ZL: ekin is of the length gvec%ng total G
  allocate(ekin (gvec%ng))
  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
    ! ZL: this includes GPP
    ! isrtq is of the length of all gvec pool
    allocate(isrtq (gvec%ng))
    allocate(isrtqi (gvec%ng))
    allocate(ind (gvec%ng))
    allocate(indinv (gvec%ng))
    allocate(ph (gvec%ng))
    if (peinf%inode .eq. 0) then
      allocate(epsmpi%isrtq (gvec%ng,sig%nq))
      allocate(epsmpi%isrtqi (gvec%ng,sig%nq))
    endif
    allocate(epsmpi%qk (3,sig%nq))
    allocate(epsmpi%nmtx (sig%nq))
    allocate(epsmpi%inv_igp_index (epsmpi%ngpown_max))
  endif
  call timing%stop(timing%read_neps)
!----------------------------
! Read eps^-1 from eps0mat/epsmat
!
! JRD: The matrices are read in from eps0mat/epsmat files and writen
! to temporary INT_EPS files on unit iunit_eps. The q->0 matrix is not
! symmetrized. The wavevector q is also read from subroutine epscopy.
  call timing%start(timing%epscopy)
  epshead = 0.0d0
  if (sig%freq_dep/=-1) then
    call logit('Calling epscopy')
    ! FHJ: distribute columns of epsinv with standard ScaLAPACK 1d block-column
    ! layout. For convenience, the ScaLAPACK auxiliary functions are defined
    ! in BerkeleyGW even if you compile the code without ScaLAPACK.
    ! TODO: get rid of array inv_igp_index, it can be computed on the fly.
    epsmpi%inv_igp_index = 0
    do igp_loc = 1, epsmpi%ngpown
      igp = INDXL2G(igp_loc, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
      epsmpi%inv_igp_index(igp_loc) = igp
    enddo
    call epscopy(crys,gvec,sig,neps,epsmpi,epshead,iunit_eps,fne)
  endif
  call timing%stop(timing%epscopy)
  if (sig%subsample) then
    allocate(sig%subweights (sig%nq0))
    call open_file(666, file='subweights.dat', form='formatted', status='old')
    read(666,*) nq0
    if (nq0/=sig%nq0) then
      call die('Inconsistency between nq0 in subweights.dat and eps0mat.h5 files.')
    endif
    do irq=1,sig%nq0
      read(666,*) sig%subweights(irq)
    enddo
    call close_file(666)
    if (peinf%inode==0) write(6,'(/,1x,a,f0.6,/)') &
      'Sum of subweights before renormalization: ', sum(sig%subweights(:))
    sig%subweights(:) = sig%subweights(:) / sum(sig%subweights(:))
    ! FHJ: Choose the smallest cutoff that doesn`t include a lattice vector
    ! in a periodic direction. FIXME: there might be a combination, such as
    ! |Gx+Gy|^2, that gives a smaller cutoff than |Gx|^2 or |Gy|^2 alone!
    subcut = INF
    do ii=1,3
      if (sig%qgrid(ii)>1) subcut = min(subcut,crys%bdot(ii,ii)*(1d0-TOL_SMALL))
    enddo
    if (all(sig%qgrid<=1)) call die('Can`t do subsampling for molecules.', only_root_writes=.true.)
    if (peinf%inode==0) write(6,'(1x,a,f0.9)') &
      'Cutoff for subsampled q-points: ', subcut
  endif
!----------------------------
! Generate full Brillouin zone from irreducible wedge q -> gr%f
  call timing%start(timing%fullbz)
  ! gr%nr: number in reduced zone
  ! gr%nf: number in full zone
  ! gr%r: points in reduced zone
  ! gr%f: points in full zone
  gr%nr = sig%nq
  allocate(gr%r (3, sig%nq))
  ! reduced q points
  gr%r(1:3,1:sig%nq) = sig%qpt(1:3,1:sig%nq)
  ! wigner_seitz = .false. for Epsilon, Sigma, use usual "box" BZ
  ! = .true. for BSE, use wigner seitz cell
  call fullbz(crys,syms,gr,syms%ntran,skip_checkbz,wigner_seitz=.false.,paranoid=.true.,nfix=sig%nq0)
  qshift(:)=0.0d0
  if (sig%freq_dep.eq.-1) then
    ! for Hartree-Fock, there is no epsmat/eps0mat file
    tmpfn="sigma.inp"
  else
    if (sig%igamma.ne.0) then
      tmpfn='eps0mat'
    else
      tmpfn='epsmat'
    endif
  endif
  if (.not. skip_checkbz) then
    !FHJ: TODO: ignore change checkbz to support nq0
    !call checkbz(gr%nf,gr%f,sig%qgrid,qshift,crys%bdot,tmpfn,'q',.false.,sig%freplacebz,sig%fwritebz,nfix=sig%nq0)
    call checkbz(gr%nf,gr%f,sig%qgrid,qshift,crys%bdot,tmpfn,'q',.false.,sig%freplacebz,sig%fwritebz)
  endif
  call timing%stop(timing%fullbz)
  if (peinf%inode==0) then
    write(6,'(1x,a,i0)') 'Number of k-points in WFN_inner: ', kp%nrk
    write(6,'(1x,a,i0)') 'Number of k-points in the full BZ of WFN_inner: ', gr%nf
  endif
! ZL: above are all preparing procedures for sigma and EP perturbed self-energy calculations
!-------- Start computation of sigma operator ----------------------------------
! ZL: here is general start
  fact = 1D0/(dble(gr%nf-sig%nq0+1)*crys%celvol)
  coulfact = 8D0*PI_D/(dble(gr%nf-sig%nq0+1)*crys%celvol)
!----------------------------
! Initialize distribution of epsilon
  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
    allocate(eps (neps,epsmpi%ngpown))
    allocate(epstemp (neps))
  endif
  if (sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
    if(sig%do_sigma_subspace) then
      !MDB only static epsinv will be calculated in the original PW basis
      nfreq_fixwings = 1
      allocate(epsR (neps,epsmpi%ngpown,1))
    else
      nfreq_fixwings = sig%nFreq
      allocate(epsR (neps,epsmpi%ngpown,sig%nFreq))
      if (sig%need_advanced) then
        allocate(epsA (neps,epsmpi%ngpown,sig%nFreq))
      endif
    end if
  endif
  if (kp%nspinor == 2) then
    call require_reference(REF_Barker2018)
  endif
!---------- Check grid for uniformity
! icutv: truncation type
  if(peinf%inode == 0) call checkgriduniformity(sig%qgrid, crys, sig%icutv)
!-------- Loop over kpoints rkn-------------------------------------------------
! nkn: number of k-points on which to calculate Sigma (from sigma.inp)
  do ika=1,sig%nkn
!----------------------------
! Read wavefunctions for rkn (if sig%nkn.gt.1)
! (else it is already in wfnk)
! ZL: indkn: mapping of k-points from sigma.inp to those in kp%rk from WFN files
! get the kpoint index in the WFN file
!
    ! ==================================================================
    ! ZL: begin preparation of wfnk
    ! in ika for outer
    ! in ikn for inner
    ! ZL: the reason for this lengthy generation of wfnk is that wfnkmpi is
    ! distributed in input.f90
    ! ika is ik for outer in sigma.inp, ikn is ik for inner in kp%rk(1:3, 1:kp%nrk)
    ikn = sig%indkn(ika)
    call timing%start(timing%wf_comm)
    if(sig%nkn.gt.1) then
      ! ZL: ntband !< number of bands in dynamical sigma summation
      ! nvband !< number of bands in bare exchange
      nbandi=sig%ntband
      ! ZL: in input.f90:
      ! wfnkmpi stores in the order of sigma input, outer wfn
      ! wfnkqmpi stores wavefunction read in, inner wfn
      wfnk%nkpt=wfnkmpi%nkptotal(ika)
      ! ndiag_max=sig%ndiag/npools
      ! ndv = ndiag * ngk
      ! ZL: def in input.f90
      ! wfnkqmpi%nkptotal(irk) = kp%ngk(irk)
      ! wfnkqmpi%isort(1:kp%ngk(irk),irk) = isort(1:kp%ngk(irk))
      wfnk%ndv=peinf%ndiag_max*wfnk%nkpt
      ! ZL: only first ngk numbers are effective
      ! this is done in read_wavefunctions() in input.f90, with findvector
      wfnk%isrtk(1:wfnk%nkpt)=wfnkmpi%isort(1:wfnk%nkpt,ika)
      ! ZL: in input.f90: allocate(wfnkmpi%qk (3,sig%nkn)), nkn: kpoints
      ! in outer wfn
      ! ZL: in read_wavefunctions():
      ! do irk=1,kp%nrk
      ! qk(:)=kp%rk(:,irk)
      ! then others copy this qk
      ! therefore qk(1:3) stores coordinates of k-points
      qk(1:3)=wfnkmpi%qk(1:3,ika)
      ! ZL: in ../Common/typedefs.f90
      ! el(:,:,:) !< band energies (band, kpoint, spin)
      ! elda(:,:,:) !< band energies before eqp correction
      wfnk%ek(1:sig%ntband,1:sig%nspin) = wfnkmpi%el(1:sig%ntband,1:sig%nspin,ika)
      wfnk%elda(1:sig%ntband,1:sig%nspin) = wfnkmpi%elda(1:sig%ntband,1:sig%nspin,ika)
      allocate(wfnk%zk (wfnk%ndv,sig%nspin*kp%nspinor))
      wfnk%zk=0.0d0
      ! ZL: here k loops over spin index (unbelievably!)
      do k=1,sig%nspin*kp%nspinor
        wfnk%zk(1:wfnk%ndv,k)=wfnkmpi%cg(1:wfnk%ndv,k,ika)
        ! wfnk done with collection from wfnkmpi (outer), which is distributed
        ! ZL: outer wfn is stored in wfnk when in calculations
      enddo ! k=1,sig%nspin*kp%spinor
    endif ! sig%nkn.gt.1
    call timing%stop(timing%wf_comm)
    if(peinf%inode.eq.0) then
      write(6,*)
      call print_dealing_with(ika, sig%nkn, kp%rk(:,ikn), 'k')
    endif
    ! ZL: end preparation of wfnk
    ! in ika for outer
    ! in ikn for inner
    ! ==================================================================
    call timing%stop(timing%wf_comm)
    ! =====================================================================
    ! =================== done wfnk_phonq from wfnk_phonq_mpi=====================
    ! =====================================================================
    ! ZL: for DEBUG use ONLY!
    ! wfnk_phonq passes test
! if(sig%elph) ep_debug = .true.
! if(ep_debug) then
! ikn = ik_phonq_idx
! wfnk = wfnk_phonq
! endif
!----------------------------
! Initialize Matrix Elements
! ZL: here is the start of the calculation of sigma matrix elements
    alda=0.0d0
    ! alda = vxc
    ! ZL: initialization for kih
    akih=0.0d0
    ax=0.0d0
    achcor(:,:)=(0.0d0,0.0d0)
    asig_imag(:,:)=(0.0d0,0.0d0)
    if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
      achtcor_n1(:) = 0.0d0
    endif
    achcor_n1(:,:,:) = 0.0d0
    if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
      ! array size: allocate(asx (3,sig%ndiag+sig%noffdiag,sig%nspin))
      asx(:,:,:)=0.0d0
      ach(:,:,:)=0.0d0
      asig(:,:)=0.0d0
      ach_n1(:,:,:)=0.0d0
    endif
    if (sig%freq_dep.eq.2) then
      asxDyn(:,:,:)=(0.0d0,0.0d0)
      achDyn(:,:,:)=(0.0d0,0.0d0)
      achDyn_cor(:,:,:)=(0.0d0,0.0d0)
      achDyn_corb(:,:,:)=(0.0d0,0.0d0)
      ach2Dyn(:,:,:)=(0.0d0,0.0d0)
      asigDyn(:,:)=(0.0d0,0.0d0)
      achD_n1(:,:,:)=(0.0d0,0.0d0)
    endif
!----------------------------
! Read matrix elements of Vxc from file vxc.dat
! or compute them on the fly from Vxc potential
    call timing%start(timing%vxc)
    if(sig%use_vxcdat .and. .not. sig%is_EXX) then
      if(peinf%inode == 0) write(6,*) 'Reading vxc.dat'
      call open_file(120,file='vxc.dat',form='formatted',status='old')
      qk(:)=INF
      ierr=0
      do while (ierr.eq.0)
        call read_matrix_elements_type(120, ierr, qk, sig, alda)
        if (all(abs(kp%rk(1:3,ikn)-qk(1:3)) .lt. TOL_Small)) exit
      enddo
      call close_file(120)
! Check k-point
      if(any(abs(kp%rk(1:3,ikn)-qk(1:3)) .ge. TOL_Small)) then
        call die('cannot find k-point in vxc.dat', only_root_writes = .true.)
      endif
! ZL: IMPORTANT: here is where the units of diag and offdiag become different
! Divide by ryd for diag
! this will be undone by shift_energy routines later
      do s2=1,sig%nspin
        do in=1,sig%ndiag
        alda(in,s2) = alda(in,s2)/ryd ! ZL: if directly read from vxc.dat
        enddo
      enddo
    elseif (.not.sig%sigma_correction .and. .not. sig%is_EXX .and. .not.sig%use_kihdat) then ! not using vxc.dat
      ! ZL: here we calculate mtxel_vxc
      call logit('Calling mtxel_vxc')
      call mtxel_vxc(kp,gvec,sig,wfnk,wfnkoff,alda,1) ! ZL: 1 for VXC
    elseif (sig%use_kihdat) then
      ! ZL: conflict check between kihdat and vxcdat has been done in inread.f90
      if(peinf%inode == 0) write(6,*) 'Reading kih.dat'
      call open_file(876,file='kih.dat',form='formatted',status='old')
      qk(:)=INF
      ierr=0
      do while (ierr.eq.0)
        call read_matrix_elements_type(876, ierr, qk, sig, akih)
        if (all(abs(kp%rk(1:3,ikn)-qk(1:3)) .lt. TOL_Small)) exit
      enddo
      call close_file(876)
! Check k-point
      if(any(abs(kp%rk(1:3,ikn)-qk(1:3)) .ge. TOL_Small)) then
        call die('cannot find k-point in vxc.dat', only_root_writes = .true.)
      endif
! ZL: IMPORTANT: again the units of diag and offdiag become different
      do s2=1,sig%nspin
        do in=1,sig%ndiag
        akih(in,s2) = akih(in,s2)/ryd ! ZL: if directly read from kih.dat
        enddo
      enddo
    endif
    call timing%stop(timing%vxc)
!----------------------------
! Read ax from existing data
! ZL: generally it is directly calculated
    if(sig%use_xdat .and. xflag .and. (.not. sig%coul_mod_flag)) then
      if(sig%elph) then
        call die('EP calculation cannot read x.dat')
      endif
      if(peinf%inode == 0) write(6,*) 'Reading x.dat'
      call open_file(119,file='x.dat',form='formatted',status='old')
      qk(:)=INF
      ierr=0
      do while (ierr.eq.0)
        call read_matrix_elements_type(119, ierr, qk, sig, ax)
        if (all(abs(kp%rk(1:3,ikn)-qk(1:3)) .lt. TOL_Small)) exit
      enddo
      ax(:,:) = ax(:,:) * sig%xfrac
      call close_file(119)
! Check k-point
      if(any(abs(kp%rk(1:3,ikn)-qk(1:3)) .ge. TOL_Small)) then
        call die('cannot find k-point in x.dat', only_root_writes = .true.)
      endif
! Divide by ryd for diag
      do s2=1,sig%nspin
        do in=1,sig%ndiag
          ax(in,s2) = ax(in,s2)/ryd
        enddo
      enddo
    endif ! using x.dat
!----------------------------
! Find subgroup which leaves kn invariant
! Indices of group operations in subgroup stored in array indsub
! stored in structure syms
! ZL: for phonon purpose, there are k, q, phonq points, may ignore this for now
! we now disable this for EP
    ! ZL: note that this has already been dealt with in input.f90
    ! ZL: here get qk(:), and qk(:) is the coord of ikn, inner wfn kpoint
    qk(:) = kp%rk(:,ikn)
    call timing%start(timing%subgrp)
    if (sig%qgridsym) then
      call subgrp(qk,syms)
    else
      syms%ntranq=1
      syms%indsub(1)=1
      syms%kgzero(1:3,1)=0
    endif
    call timing%stop(timing%subgrp)
!----------------------------
! Reduce qpoints with respect to group of kn
! Keep track of number of q-points equivalent
! to give qpoint in irr-bz neq(irq)
!
! Keep track of q-point in the set given in epsmat
! RQ = R(Q) + G0 with transformation R and umklapp G0
! In order to unfold the inverse dielectric matrix from Q to RQ
!
! The q-points are the points for which we have epsilon
! (should be the unshifted irreducible grid)
    ! ZL: grid%nr !< number in reduced zone
    ! grid%nf !< number in full zone
    allocate(indrq (gr%nf)) ! full zone
    allocate(neq (gr%nf)) ! full zone
    allocate(itnrq (gr%nf)) ! full zone
    allocate(rq (3,gr%nf)) ! full zone
    allocate(kg0 (3,gr%nf)) ! full zone
    call timing%start(timing%irrbz)
    call irrbz(syms,gr%nf,gr%f,nrq,neq,indrq,rq,sig%nq,sig%qpt,itnrq,kg0,nfix=sig%nq0)
    call timing%stop(timing%irrbz)
   if ( mtxel_algo /= CPU_ALGO ) then
      call setup_FFT_sizes(gvec%FFTgrid,Nfft,dummy_scale)
      acc_mtxel_sig%mtxel_band_block_size = sig%acc_mtxel_band_block_size
      acc_mtxel_sig%mtxel_band_block_size = MIN(acc_mtxel_sig%mtxel_band_block_size, peinf%ntband_node)
      acc_mtxel_sig%mtxel_band_block_size = MAX(1, acc_mtxel_sig%mtxel_band_block_size)
      call allocate_acc_mtxel_sig ( Nfft, kp%ngkmax, kp%ngkmax, gvec%ng, kp%ngkmax, mtxel_algo )
      if (peinf%inode == 0) then
        write(6,'(1x,A)')
        write(6,'(1x,A)') 'Using GPU support via OpenACC/OMP-Target implementation'
        write(6,'(1x,A,I6)') 'Band block size in MTXEL:', acc_mtxel_sig%mtxel_band_block_size
        if ( sigma_gpp_algo == OMP_TARGET_ALGO .and. (sig%freq_dep==1 .or. sig%freq_dep==3) ) then
          write(6,'(1x,A)') 'Using GPU support via OMP-Target for Sigma-GPP'
          write(6,'(1x,A,I6)') 'Band block size:', sig%gpp_band_block_size
          write(6,'(1x,A,I6)') 'G-vec block size:', sig%gpp_ig_block_size
        end if
      end if
   end if
!!---------- Loop over k-points in irr bz with respect to kn (rq) --------------
    irq_min = 1
! FHJ: For subsampling calculations, we calculate the self energy for the first
! q-point (irq_==0, irq==1) using the input epsilon cutoff (ecuts), but only a
! correction to the self energy for the other q-points, which is calculated
! using a smaller cuttof. The way this is implemented is by looping twice over
! the first q-point: first (irq_==0, irq==1) using a weight of 1, and then
! (irq_==1, irq==1) using a weight of subweight(1)-1. All other subsampled
! q-points are calculated with the weight subweight(irq).
    if (sig%subsample) irq_min = 0
    call progress_init(prog_info, 'calculating Sigma', 'block', &
      (nrq-irq_min+1)*(peinf%ndiag_max+peinf%noffdiag_max)*sig%nspin)
    !======================================================================
    ! ZL: here is the start loop over q points
    do irq_ = irq_min, nrq
      irq = irq_
      is_subq = .false.
      ! FHJ: is this a q-point used in the subsampling of the voronoi cell for q==0?
      if (sig%subsample.and.irq<=sig%nq0) then
        is_subq = .true.
        if (irq==0) irq=1
      endif
      if (peinf%verb_debug .and. peinf%inode.eq.0) then
        write(6,60) irq,nrq
        write(6,*)
      endif
60 format(/,3x,'qpoint',i5,' out of ',i5)
! Compute energies: |q+g|**2
! ZL: stored as ekin(ig)
      if (is_subq) then
        call kinetic_energies(gvec, crys%bdot, ekin)
      else
        call kinetic_energies(gvec, crys%bdot, ekin, qvec = rq(1:3, irq))
      endif
! Sort ekin in ascending order for this q
! The indices are placed in array isrtrq
      !ZL: Note that ekin of |q+g|**2 is used
      call sortrx(gvec%ng, ekin, isrtrq, gvec = gvec%components)
      if ((sig%freq_dep.eq.0.or.sig%exact_ch.eq.1).and.irq_==irq_min) then
        isrtrqi=0
        do j=1,gvec%ng
          if (isrtrq(j).ge.1.and.isrtrq(j).le.gvec%ng) &
            isrtrqi(isrtrq(j))=j
        enddo
      endif
      !=========================================================
      ! ZL: calculate cutoffs, unchanged
! Compute cutoff in sums over G,G`
! definition: real(DP) :: ecutb !< energy cutoff of bare coulomb interaction in Ry
! definition: real(DP) :: ecuts !< energy cutoff of screened coulomb interaction in Ry
      ncoulb = gcutoff(gvec%ng, ekin, isrtrq, sig%ecutb)
      ncouls = gcutoff(gvec%ng, ekin, isrtrq, sig%ecuts)
      if (is_subq.and.irq_>0) then !FHJ: see comment before loop over irq_
        ncouls = min(ncouls, gcutoff(gvec%ng, ekin, isrtrq, subcut))
      endif
      ! ncoul is the maximum possible value
      ncoul = max(ncouls,ncoulb)
      if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then
        if(irq_==irq_min) ncoulch = ncoul
      else
        ncoulch = 0
      endif
      ! ZL: important
      ! this condition is assumed later in the code (e.g. wpeff array sizes), so we must check it
      if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
        if (ncouls.gt.neps) then
          write(tmpstr,'(a,i6,a,i6)')&
            "screened Coulomb cutoff is bigger than epsilon cutoff"//CHAR(10)&
            //' ncouls = ',ncouls,' neps = ',neps
          ! NB: CHAR(10) is the carriage return.
          call die(tmpstr, only_root_writes = .true.)
        endif
      endif
      ! ZL: done cutoffs, unchanged
      !=========================================================
!----------------------------
! Allocate arrays for q-point rq(:,irq)
      ! ZL: V(G)
      allocate(vcoul (ncoul))
      ! ZL: aqs - matrix elements M(N_G, N_band) now for a given k and q point (in the loops)
      allocate(aqs (ncoul,peinf%ntband_max))
      if (sig%noffdiag.gt.0) then
        allocate(aqsaug (ncoul,peinf%ntband_max,sig%ndiag,sig%nspin))
      endif
      ! ZL: allocating for EP
      ! ZL: for GPP, we dont consider this if-statement
      if ((sig%freq_dep.eq.0.or.sig%exact_ch.eq.1).and.irq_==irq_min) then
        ! ZL: EP implemented for GPP only for now
        allocate(aqsch (ncoulch))
        if (nrq.gt.1) then
          allocate(aqsaugchd (ncoulch,peinf%ndiag_max,sig%nspin))
          if (sig%noffdiag.gt.0) then
            allocate(aqsaugcho (ncoulch,peinf%noffdiag_max,sig%nspin))
          end if
        endif
      endif
      ! ZL: nm is index for q
      nm = indrq(irq)
!!!------- Calculate Vcoul -----------------------------------------------------
      call timing%start(timing%vcoul)
      vq=rq(:,irq)
      qlen = sqrt(DOT_PRODUCT(vq,MATMUL(crys%bdot,vq)))
      if(sig%freq_dep /= -1) call check_screening_trunc(sig%icutv,sig%iscreen,sig%q0vec,crys%bdot)
      iparallel=1
      avgcut = sig%avgcut
      ! ZL: MC subsampling
      ! FHJ: Disable MC averages when we perform the subsampling of the voronoi cell
      if (is_subq) avgcut = 0d0
      if (.not. sig%coul_mod_flag) then
        call vcoul_generator(sig%icutv,sig%truncval,gvec,crys%bdot,crys%celvol, &
          gr%nf-sig%nq0+1,ncoul,isrtrq,sig%iscreen,vq,sig%q0vec,vcoul, &
          sig%iwritecoul,iparallel,avgcut,oneoverq,sig%qgrid,epshead, &
          work_scell,sig%averagew,sig%wcoul0)
      endif
      fact = 1D0/(dble(gr%nf-sig%nq0+1)*crys%celvol)
      coulfact = 8D0*PI_D/(dble(gr%nf-sig%nq0+1)*crys%celvol)
      ! FHJ: we won`t actually use sig%wcoul0 when sig%subsample==.true.
      ! Rescale vcoul according to the subsample weights
      if (is_subq) then
        weight = sig%subweights(nm)
        if (irq_==0) then ! FHJ: see comment before loop over irq_
          weight = 1d0
        elseif (irq_==1) then
          weight = weight - 1d0
        endif
        fact = fact*weight
        coulfact = coulfact*weight
        if (irq/=nm) then
          if (peinf%inode==0) write(0,'(a,i0,a,i0)') 'ERROR: irq=',irq,' nm=',nm
          call die('Inconsistent indices for a sub-sampled q-point', only_root_writes=.true.)
        endif
      endif
      ! ZL: now we have V(G)
      do ig = 1, ncoul
        vcoul(ig)=fact*vcoul(ig)
      enddo
      if (ika.eq.1.and.irq_==irq_min) then
        sig%wcoul0 = sig%wcoul0 * fact
      endif
      call timing%stop(timing%vcoul)
!!!------- Read inverse dielectric matrix for q-point rq(:,irq) ----------------
! ZL: for EP, epsinv is unchanged, which determines wtilde and Omega,
! and they are for the same q point (or p point in the formalism)
      call timing%start(timing%epsread)
      if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
        q0len = sqrt(DOT_PRODUCT(sig%q0vec,MATMUL(crys%bdot,sig%q0vec)))
!----------------------------
! Processor 0 read eps^-1 for this q and broadcast to others
!
! Note: the total number of g-vectors used during
! computation of the inverse dielectric matrix
! may be different than in present calculation
! although the sets must coincide for small g
        ngq = gvec%ng
        nmtx = epsmpi%nmtx(nm)
        if (peinf%inode .eq. 0) then
          ! ZL: isrtq now represents the order of q-points in eps
          isrtq(:) = epsmpi%isrtq(:,nm)
          isrtqi(:) = epsmpi%isrtqi(:,nm)
        endif
        qk(:) = epsmpi%qk(:,nm)
! Find g=0 in main gvec list and eps gvector list
        ! write explicitly to avoid possible warning about array temporary
        nullvec(1:3) = 0
        call findvector(iout,nullvec,gvec)
        iout = isrtqi(iout)
        if (peinf%verb_debug .and. peinf%inode==0) write(6,*) 'Reading Eps Back'
! JRD: XXX This is sort of a waste of memory... Can we use pointers for this sort of thing?
        if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
          eps(:,:)=epsmpi%eps(:,:,nm)
        else ! sig%freq_dep.eq.2 .or. 3
          epsR(:,:,:)=epsmpi%epsR(:,:,:,nm)
          if (sig%need_advanced) then
            epsA(:,:,:)=epsmpi%epsA(:,:,:,nm)
          endif
        endif
! CHP: By first distributing and then doing the wing fix, we can let
! all the processors work together, thus, saving some time.
! ZL: TODO q0->0 behavior not determined yet, to be fixed or double checked
!-------------------------------------------------------------------------------
! Fix wing divergence for semiconductors and graphene
! This should really be done for all "|q+G| < avgcut" - but for now,
! it is done if "|q| < avgcut and G=0"
        ! make sure corrections for wings are zero at the beginning
        if (sig%do_sigma_subspace) then
           sig%epssub%eps_wings_correction_cols(:,:) = (0.0d0, 0.0d0)
           sig%epssub%eps_wings_correction_rows(:,:) = (0.0d0, 0.0d0)
           ! save actual nm
           sig%epssub%actual_nm = nm
        end if
        if (.not.sig%subsample) then
        ngpown_q = NUMROC(nmtx, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
          if (sig%do_sigma_subspace) then
            ! initialize wings corrections
            ipe_wing = INDXG2P( iout, epsmpi%nb, iproc_dum, 0, peinf%npes_pool)
            if (ipe_wing == peinf%pool_rank) then
              ! make a copy of the unmodified wing
              my_pos_loc_wing = INDXG2L(iout, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
              sig%epssub%eps_wings_correction_cols(1:nmtx,1:sig%nFreq) = &
              sig%epssub%eps_wings_cols(1:nmtx,1:sig%nFreq,nm)
            end if
            do igp_loc = 1, ngpown_q
               igp = INDXL2G(igp_loc, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
               sig%epssub%eps_wings_correction_rows(igp,1:sig%nFreq) = &
               sig%epssub%eps_wings_rows(igp,1:sig%nFreq,nm)
            end do
            ! MDB calculate correction factors for fixing wings (here we try
            ! to avoid copying stuff from other modules even if less efficient)
            allocate(epsR_corrections (neps,epsmpi%ngpown,1))
            epsR_corrections = (1.0,0.0)
          end if
! fix wings here
        do igp_loc = 1, ngpown_q
          igp = INDXL2G(igp_loc, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
          if (nm .eq. 1) then
            q0flag=.true.
            if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
              call fixwings(vcoul(1),sig%wcoul0,eps(:,igp_loc),sig%icutv, &
                sig%iscreen,igp,nmtx,iout,q0len,oneoverq,fact,q0flag,sig%averagew,crys%bdot)
            else ! sig%freq_dep.eq.2 .or. 3
! JRD XXX bad locality
              call fixwings_dyn(vcoul(1),epsR(1:nmtx,igp_loc,:), &
                sig%icutv,sig%iscreen,igp,nfreq_fixwings,nmtx,iout,q0len,oneoverq,fact,q0flag,crys%bdot)
              if (sig%need_advanced) then
                call fixwings_dyn(vcoul(1),epsA(1:nmtx,igp_loc,:), &
                  sig%icutv,sig%iscreen,igp,nfreq_fixwings,nmtx,iout,q0len,oneoverq,fact,q0flag,crys%bdot)
              endif
              if (sig%do_sigma_subspace) then
                call fixwings_dyn(vcoul(1),epsR_corrections(1:nmtx,igp_loc,:), &
                  sig%icutv,sig%iscreen,igp,1,nmtx,iout,q0len,oneoverq,fact,q0flag,crys%bdot)
              end if
            endif
          else if (qlen**2 .lt. sig%avgcut) then
            q0flag=.false.
            if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
              call fixwings(vcoul(1),sig%wcoul0,eps(:,igp_loc),sig%icutv, &
                sig%iscreen,igp,nmtx,iout,qlen,oneoverq,fact,q0flag,sig%averagew,crys%bdot)
            else ! sig%freq_dep.eq.2 .or. 3
! JRD XXX bad locality
              call fixwings_dyn(vcoul(1),epsR(1:nmtx,igp_loc,:), &
                sig%icutv,sig%iscreen,igp,nfreq_fixwings,nmtx,iout,qlen,oneoverq,fact,q0flag,crys%bdot)
              if (sig%need_advanced) then
                call fixwings_dyn(vcoul(1),epsA(1:nmtx,igp_loc,:), &
                  sig%icutv,sig%iscreen,igp,nfreq_fixwings,nmtx,iout,qlen,oneoverq,fact,q0flag,crys%bdot)
              endif
              if (sig%do_sigma_subspace) then
                 call fixwings_dyn(vcoul(1),epsR_corrections(1:nmtx,igp_loc,:), &
                   sig%icutv,sig%iscreen,igp,1,nmtx,iout,qlen,oneoverq,fact,q0flag,crys%bdot)
              end if
            endif
          endif
        enddo ! igp
          if (sig%do_sigma_subspace) then
            ! calculte corrections
            ipe_wing = INDXG2P( iout, epsmpi%nb, iproc_dum, 0, peinf%npes_pool)
            if (ipe_wing == peinf%pool_rank) then
              my_pos_loc_wing = INDXG2L(iout, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
              do ifreq = 1, sig%nFreq
                sig%epssub%eps_wings_correction_cols(1:nmtx,ifreq) = &
                sig%epssub%eps_wings_correction_cols(1:nmtx,ifreq) * &
                (epsR_corrections(1:nmtx,my_pos_loc_wing,1) - 1.0d0)
              end do
              ! broadcast to your fellow
            else
              ! receive the wing / col
            end if
            ! wing row
            do igp_loc = 1, ngpown_q
              igp = INDXL2G(igp_loc, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
              sig%epssub%eps_wings_correction_rows(igp,1:sig%nFreq) = &
              sig%epssub%eps_wings_correction_rows(igp,1:sig%nFreq)*(epsR_corrections(iout,igp_loc,1) - 1.0d0)
              if(igp == iout) then
                ! here compensate for the 1.0 on diagonal
                sig%epssub%eps_wings_correction_rows(igp,1:sig%nFreq) = &
                sig%epssub%eps_wings_correction_rows(igp,1:sig%nFreq) + epsR_corrections(iout,igp_loc,1)
              end if
            end do
            ! sum up
            if(associated(epsR_corrections))then;deallocate(epsR_corrections);nullify(epsR_corrections);endif
          end if ! subspace
        endif !.not.sig%subsample
        call logit('Read eps from memory')
      endif ! sig%freq_dep
      if (sig%freq_dep.eq.-1) then
        ngq=0
        nmtx=0
        qk(:)=sig%qpt(:,nm)
      endif
      call timing%stop(timing%epsread)
      if (peinf%verb_debug .and. peinf%inode==0) then
        write(6,'(3(a,i10))') 'nmtx =',nmtx,' ncouls =',ncouls
        write(6,*)
      endif
      if(nmtx.gt.neps) then
        call die('nmtx.gt.neps')
      endif
! Check q-vector
      if(any(abs(sig%qpt(1:3, nm) - qk(1:3)) .gt. TOL_SMALL)) then
        write(0,*) peinf%inode,sig%qpt(:,nm),qk(:)
        call die('q-vector check wrong')
      endif
      if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
        if(ncouls .gt. nmtx) ncouls = nmtx
      endif
      if (peinf%verb_debug .and. peinf%inode.eq.0) then
        if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
          write(6,100) (rq(i,irq),i=1,3),ncouls,dble(eps(1,1))
        endif
        if (sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
          write(6,100) (rq(i,irq),i=1,3),ncouls,dble(epsR(1,1,1))
        endif
      endif
100 format(3x,'q=',3f8.5,2x,'n=',i6,2x,'head of epsilon inverse =',f12.6,/)
! Map g-vectors required for eps**(-1)(r(q)) to those
! for the known matrix eps**(-1)(q) and calculate phases
      itran = itnrq(irq)
      kg(:) = kg0(:,irq)
      call timing%start(timing%gmap)
      call logit('Calling gmap')
      if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
        ind(:)=0
        indinv(:)=0
        ph(:)=0.0d0
        call gmap(gvec,syms,ncouls,itran,kg,isrtrq,isrtqi,ind,ph, &
         sig%die_outside_sphere) ! TROUBLE
        do ig = 1, gvec%ng
          if (ind(ig) .gt. 0 .and. ind(ig) .le. gvec%ng) then
            indinv(ind(ig))=ig
          endif
        enddo
      endif
      call timing%stop(timing%gmap)
!!!------- Done reading inverse dielectric matrix for q-point rq(:,irq) --------
!=================================================================
!ZL: generate wfnkq from wfnkqmpi (WFN read in) for rkq = rk - rq
! the format of wfnkq is wfnkqinfo
!--------------------
! Generate needed wavefunctions for rkq = rkn - rq
! stored in derived type wfnkq
      rkq(1:3) = kp%rk(1:3, ikn) - rq(1:3, irq)
      ! FHJ: when we perform subsampling, we use the actual q-vector for the
      ! Coulomb potential, but we pretend we have q=0 for the matrix elements.
      if (is_subq) rkq(1:3) = kp%rk(1:3, ikn)
      call timing%start(timing%genwf)
      call logit('Calling genwf')
      if (.not.is_subq.or.irq_==0) then ! FHJ: subsampling uses the same WFNs and matrix elements
        call genwf_mpi(rkq, syms, gvec, crys, kp, sig, wfnkq, wfnkqmpi)
      endif
      call timing%stop(timing%genwf)
!ZL: done generating wfnkq
!=================================================================
!=================================================================
!!-------- Loop over spins for which Sigma is computed -------------------------
! ZL: Now get into core part of sigma
      do ispin=1,sig%nspin
!!-------- Loop over bands for which diag Sigma is computed --------------------
! Bands are relabelled according to sig%diag(1:sig%ndiag)
! Loop is parallelized according to peinf%index_diag(1:peinf%ndiag_max)
! ZL: Loop over diagonal bands, in parallel
        do in=1,peinf%ndiag_max
          call progress_step(prog_info)
          if (peinf%verb_debug .and. peinf%inode.eq.0) then
            if (peinf%npools.eq.1) then
              write(6,999) peinf%index_diag(in)+(ispin-1)*sig%ndiag, &
                peinf%ndiag_max*peinf%npools*sig%nspin
            else
              write(6,997) peinf%index_diag(in)+(ispin-1)*sig%ndiag, &
                peinf%index_diag(in)+(ispin-1)*sig%ndiag+peinf%npools-1, &
                peinf%ndiag_max*peinf%npools*sig%nspin
            endif
          endif
999 format(1x,"Computing Sigma diag",i4,1x,"of",i4)
997 format(1x,"Computing Sigma diag",i4,1x,"to",i4,1x,"of",i4)
          write(tmpstr,*) 'Working on band ', sig%diag(peinf%index_diag(in)), ' 1st pool'
          call logit(tmpstr)
!---------------------
! Compute planewave matrix elements of g <n,k|exp{i(q+g).r}|n1,k-q>
! Note: wfnk%zk array keeps only the bands specified in sig%diag(:)
! Must keep track of the right band label
          call logit('Calling mtxel')
          call timing%start(timing%mtxel)
          ! ZL: role of isrtrq
          ! index array for g-vectors in <nk|exp(i(q+g).r)|n1k-q>
          ! sorted with |q+g|**2
          ! in ../Common/fftw_inc.f90
          ! isrtrq is only used to find where the ncoul cutoff is
          ! when getting aqs data from already-calculated FFT
          !
          ! ZL: wfnk, wfnkq will use their own isrt in putting themselves in FFT
          ! box, and that will be cutoff-ed as well
          !
          ! ZL: aqs = <nk|exp(i(q+g).r)|n1k-q>
          ! aqs(nG,num of n1 distributed on each process)
          ! aqs is distributed over bands peinf%ntband_max
          call mtxel(in,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
          ! ZL: for EP: gvec passed to mtxel provides ng, FFTgrid, components,
          ! which are the same for WFN_inner and dWFN
          ! ZL: in kp and ep_kp, ngk are not the same
          call timing%stop(timing%mtxel)
          if (sig%noffdiag.gt.0.and.peinf%flag_diag(in)) then
            do j=1,peinf%ntband_node
              do i=1,ncoul
                ! ZL: aqsaug accumulates all aqs
                ! this is used for off-diag, because even in off-diag
                ! all matrix elements have already been calculated
                aqsaug(i,j,peinf%index_diag(in),ispin)=aqs(i,j)
                ! ZL: for EP with off-diag in band index
              enddo
            enddo
          endif
          if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then
            ! ZL: EP does not include this situation for now
            if (irq_==irq_min) then
              call logit('Calling mtxel_occ')
              call timing%start(timing%mtxel_ch)
! JRD ALL PROCS DO THIS NOW.
              call mtxel_occ(in,in,gvec,wfnk,ncoulch,isrtrq,aqsch,ispin,kp)
              call timing%stop(timing%mtxel_ch)
              if (nrq.gt.1) aqsaugchd(:,in,ispin)=aqsch(:)
            else
              aqsch(:)=aqsaugchd(:,in,ispin)
            endif
          endif
!---------------------
! Compute diag SX and CH
          if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
            call logit('Calling mtxel_cor for diagonal matrix elements')
            call timing%start(timing%mtxel_cor_tot)
            ! ZL: actually there are two M - aqsn and aqsm, following
            ! < psi_n ( k ) | Sigma_cor ( E ) | psi_m ( k ) >
            ! aqsn is the regular M directly calculated e^{i(q+G).r}
            ! aqsm is calculated in the same way, so it need to be
            ! complex-conjugated when being used
            ! so be careful when first calculating aqsm
            !
            ! ZL: in mtxel_cor(), wfnk and wfnkq only provide band energies
            ! - wfnk%ek is used as the energy argument E in Sigma(E), in
            ! general with +/- dE in GPP, for interpolation
            ! - wfnkq%eqk is used, as explicitly shown in equations, band
            ! enregy at e(k-q)
            ! ZL: when dealing with EP, just need to put the correct energy
            ! (i.e. the corresponding wfn variables) into mtxel_cor
            ! Note: we may need to average e(k+q) and e(k) for off-diag
            ! nature of EP, for E.
            ! ZL: NOTE: TODO: if regular GW is done already, maybe we can use
            ! quasiparticle energies instead of LDA energies
            !
            ! ZL Note:
            ! ind, indinv: for eps q+G gvec sort
            ! isrtrq, isrtrqi: sorting of |q+G|^2, only used to determin cutoff
            ! aqsn, aqsm: two M in such way M_n * M_m^*
            ! epsR, epsA, and *D*: are related to full frequency
            ! aqsch: for sig%freq_dep==0 or sig%exact_ch==1, not considered for EP for now
            ! achtcor: for exact_ch, not considered for EP now
            !
            ! Useful variables with out intent:
            ! asigt_imag, acht_n1, asxt, acht
            ! We also allocate for EP of the following as place holder:
            ! achtcor, achtcor_n1
            !
            ! aqsch: not even place holder for EP. Not even allocated for GPP
            !
            call mtxel_cor(peinf%index_diag(in), &
             peinf%index_diag(in),ispin,ncouls,neps,gvec,eps,ph, &
             ind,indinv,isrtrqi,isrtrq,vcoul,crys,sig,wpg,wfnk,wfnkq,ncoulch, &
             aqs,aqs,aqsch,asigt_imag,acht_n1,asxt,acht,achtcor,achtcor_n1, &
             kp%nspin,rq(:,irq),coulfact, &
             epsmpi%inv_igp_index,epsmpi%ngpown, &
             epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,1, &
             .false.) ! ZL: set ep_mtxel to be .false. because this is standard
                       ! diagonal quasiparticle energy calculation
            call timing%stop(timing%mtxel_cor_tot)
          else
            achtcor = 0.0d0
            if (sig%freq_dep/=0 .and. sig%exact_ch==1) achtcor_n1 = 0.0d0
            asigt_imag = 0.0d0
          endif
          if (peinf%flag_diag(in)) then
            ! FHJ: Store temporary arrays in other variables, and put in
            ! degeneracy factor due to star of q-points.
            jb = peinf%index_diag(in)
            achcor(jb,ispin) = achcor(jb,ispin) + neq(irq)*achtcor
            asig_imag(jb,ispin) = asig_imag(jb,ispin) + neq(irq)*asigt_imag
            if (sig%freq_dep==2) then
              asxDyn(:,jb,ispin) = asxDyn(:,jb,ispin) + neq(irq)*asxtDyn(:)
              achDyn(:,jb,ispin) = achDyn(:,jb,ispin) + neq(irq)*achtDyn(:)
              achDyn_cor(:,jb,ispin) = achDyn_cor(:,jb,ispin) + neq(irq)*achtDyn_cor(:)
              achDyn_corb(:,jb,ispin) = achDyn_corb(:,jb,ispin) + neq(irq)*achtDyn_corb(:)
              ach2Dyn(:,jb,ispin) = ach2Dyn(:,jb,ispin) + neq(irq)*ach2tDyn(:)
              achD_n1(:,jb,ispin) = achD_n1(:,jb,ispin) + neq(irq)*achtD_n1(:)
            elseif (sig%freq_dep/=-1) then
              asx(:,jb,ispin) = asx(:,jb,ispin) + neq(irq)*asxt(:)
              ach(:,jb,ispin) = ach(:,jb,ispin) + neq(irq)*acht(:)
              ach_n1(:,jb,ispin) = ach_n1(:,jb,ispin) + neq(irq)*acht_n1(:)
            endif
            if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
              achcor_n1(:,jb,ispin) = achcor_n1(:,jb,ispin) + neq(irq)*achtcor_n1(:)
            endif
          endif
!---------------------
! Compute diag bare exchange (SKIP this computation if you already know it)
          call timing%start(timing%bare_x)
          if ((.not. sig%use_xdat .and. xflag .and. (.not. sig%coul_mod_flag)) .or. (sig%coul_mod_flag &
            .and. (.not. sig%use_vxc2dat))) then
            call logit('Computing bare X')
            axt=0.0d0
            ! ZL: EP, bare perturbed exchange part
! XXX THREAD?
            do n1=1,peinf%nvband_node
              tempval = wfnkq%ekq(peinf%indext(n1),ispin) - sig%efermi
              if (tempval < sig%tol) then ! ZL: sig%tol is positive, any negative values are smaller than sig%tol
                if(abs(tempval) < sig%tol) then
                  occ=0.5 ! Fermi-Dirac distribution = 1/2 at Fermi level
                else
                  occ=1D0
                endif
                do ig=1,ncoulb
                  ! ZL: here it calculates bare X term: GV
                  axt = axt + abs(aqs(ig,n1))**2 * occ * vcoul(ig)
                  ! ZL: for EP, involve more terms
                enddo
              endif
                ! sig%ncrit = 0 and tempval > sig%tol should never happen!
            enddo
            ! if(peinf%inode.eq.0) write(6,*) "axt, axt_1, axt_2", axt, axt_ep_one, axt_ep_two
            if (peinf%flag_diag(in)) then
              ib = peinf%index_diag(in)
              ax(ib,ispin) = ax(ib,ispin) - neq(irq)*axt*sig%xfrac
            endif
          endif ! not use x.dat
          call timing%stop(timing%bare_x)
        enddo ! in (loop over bands for which we need diag Sigma)
        if (ispin.eq.sig%nspin) then
          if (.not.is_subq.or.irq_==sig%nq0) then
            if(associated(wfnkq%zkq))then;deallocate(wfnkq%zkq);nullify(wfnkq%zkq);endif
            ! ZL: EP
          endif
          if(allocated(aqs))then;deallocate(aqs);endif
          if ((sig%freq_dep.eq.0.or.sig%exact_ch.eq.1).and.irq_==nrq) then
            if (nrq.gt.1) then
              if(associated(aqsaugchd))then;deallocate(aqsaugchd);nullify(aqsaugchd);endif
            end if
          endif
        endif
!!-------- End diag band loop --------------------------------------------------
! (gsm) begin distributing aqsaug matrix elements for offdiag calculation
! $$$ inefficient communication, this should be rewritten $$$
        call timing%start(timing%mtxel_comm)
        call timing%stop(timing%mtxel_comm)
! (gsm) end distributing aqsaug matrix elements for offdiag calculation
!!-------- Loop over bands for which offdiag Sigma is computed -----------------
! Bands are relabelled according to sig%off*(1:sig%noffdiag)
! Loop is parallelized according to peinf%index_offdiag(1:peinf%noffdiag_max)
        do ioff=1,peinf%noffdiag_max
          call progress_step(prog_info)
          if (peinf%verb_debug .and. peinf%inode==0) then
            if (peinf%npools.eq.1) then
              write(6,998) peinf%index_offdiag(ioff)+(ispin-1)*sig%noffdiag, &
                peinf%noffdiag_max*peinf%npools*sig%nspin
            else
              write(6,996) peinf%index_offdiag(ioff)+(ispin-1)*sig%noffdiag, &
                peinf%index_offdiag(ioff)+(ispin-1)*sig%noffdiag+peinf%npools-1, &
                peinf%noffdiag_max*peinf%npools*sig%nspin
            endif
          endif
998 format(1x,"Computing Sigma offdiag",i4,1x,"of",i4)
996 format(1x,"Computing Sigma offdiag",i4,1x,"to",i4,1x,"of",i4)
          write(tmpstr,'(a,2i6,a)') 'Working on bands ', sig%off1(peinf%index_offdiag(ioff)), &
            sig%off2(peinf%index_offdiag(ioff)), ' 1st pool'
          call logit(tmpstr)
! <n|Sigma|m> = 0 if n and m belong to different irreducible representations
! Even without assigning representations, we can tell they are different if
! the size of the degenerate subspace is different for n and m.
! This saves time and helps enforce symmetry. -- DAS
          bExactlyZero = .false.
          if(.not. sig%wfn_outer_present .and. sig%offdiagsym .and. &
             kp%degeneracy(sig%off1(peinf%index_offdiag(ioff))-sig%ncore_excl, ikn, ispin) /= &
             kp%degeneracy(sig%off2(peinf%index_offdiag(ioff))-sig%ncore_excl, ikn, ispin)) then
            ! ZL: comment out the following line due to messy output
            !if (peinf%inode.eq.0) write(6,'(a)') 'Zero by symmetry -- not computing.'
            ! the matrix elements are zeroed at the beginning, and at the end of the loop
            ! so we can just leave them as they are
            ! JRD - We cannot cycle here because other pools may have legitimate work to do
            ! and we need to be around for the communication. Particularly inside the subroutine
            ! mtxel_cor. Thus we can`t really save time here. We could still zero out elements
            ! at bottom of loop if want. I leave it up to DAS.
            !cycle
            bExactlyZero=.true.
          endif
          ! ZL: for EP, we always calculate everything because phonq may break symmetry
          if(sig%elph) bExactlyZero = .false.
!---------------------
! Compute planewave matrix elements for exact static CH
! ZL: wfnkoff is only used for this purpose
          if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then
            if (irq_==irq_min) then
              call timing%start(timing%wf_ch_comm)
              wfnkoff%nkpt=wfnk%nkpt
              if (mod(peinf%inode,peinf%npes/peinf%npools).eq.0) then
                allocate(wfnkoff%isrtk (gvec%ng))
                wfnkoff%isrtk=wfnk%isrtk
                allocate(wfnkoff%zk (2*wfnk%nkpt,kp%nspinor))
              endif
! (gsm) begin gathering wavefunctions over pools
! $$$ inefficient communication, this should be rewritten $$$
! BAB note: loop over spinors
              do jj=1,peinf%npools
                dest=(jj-1)*(peinf%npes/peinf%npools)
                do ii=1,2
                  i=sig%offmap(peinf%index_offdiag(ioff),ii)
                  j=(i-1)/peinf%npools+1
                  source=mod(i-1,peinf%npools)*(peinf%npes/peinf%npools)
                  if (peinf%inode.eq.source.and.peinf%inode.eq.dest) then
                    do jsp=1,kp%nspinor
                      do k=1,wfnk%nkpt
                        wfnkoff%zk((ii-1)*wfnk%nkpt+k,jsp) = wfnk%zk((j-1)*wfnk%nkpt+k,ispin*jsp)
                      enddo
                    enddo
                  else
                    do jsp=1,kp%nspinor
                      do k=1,wfnk%nkpt
                        wfnkoff%zk((ii-1)*wfnk%nkpt+k,jsp) = wfnk%zk((j-1)*wfnk%nkpt+k,ispin*jsp)
                      enddo
                    enddo
                  endif
                enddo
              enddo
! (gsm) end gathering wavefunctions over pools
              call timing%stop(timing%wf_ch_comm)
              call logit('Calling mtxel_occ')
              call timing%start(timing%mtxel_ch)
! JRD Everyone does this now YYYY
              if (mod(peinf%inode,peinf%npes/peinf%npools).eq.0) then
                call mtxel_occ(1,2,gvec,wfnkoff,ncoulch,isrtrq,aqsch,1,kp)
              endif
              call timing%stop(timing%mtxel_ch)
              if (mod(peinf%inode,peinf%npes/peinf%npools).eq.0) then
                if(associated(wfnkoff%isrtk))then;deallocate(wfnkoff%isrtk);nullify(wfnkoff%isrtk);endif
                if(associated(wfnkoff%zk))then;deallocate(wfnkoff%zk);nullify(wfnkoff%zk);endif
              endif
              if (nrq.gt.1) aqsaugcho(:,ioff,ispin)=aqsch(:)
            else
              aqsch(:)=aqsaugcho(:,ioff,ispin)
            endif
          endif ! sig%freq_dep.eq.0.or.sig%exact_ch.eq.1
!---------------------
! Compute offdiag SX and CH
! ZL: EP related start from here
          if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
            call logit('Calling mtxel_cor for offdiagonal matrix elements')
            call timing%start(timing%mtxel_cor_tot)
            call mtxel_cor(sig%offmap(peinf%index_offdiag(ioff),1), &
              sig%offmap(peinf%index_offdiag(ioff),3), &
              ispin,ncouls,neps,gvec,eps,ph,ind,indinv,isrtrqi, &
              isrtrq,vcoul,crys,sig,wpg,wfnk,wfnkq,ncoulch, &
              aqsaug(:,:,sig%offmap(peinf%index_offdiag(ioff),1),ispin), &
              aqsaug(:,:,sig%offmap(peinf%index_offdiag(ioff),2),ispin), &
              aqsch,asigt_imag,acht_n1,asxt,acht,achtcor,achtcor_n1,kp%nspin, &
              rq(:,irq),coulfact, &
              epsmpi%inv_igp_index,epsmpi%ngpown, &
              epsR,epsA,achtD_n1,asxtDyn,achtDyn,achtDyn_cor,achtDyn_corb,ach2tDyn,2, &
              .false.) ! ZL: set ep_mtxel to be .false. because this is just
                        ! standard off-diagonal quasiparticle Sigma
                        ! calculation
             ! ZL: The following piece calculates off-diag (in band index) for
             ! perturbed Sigma operator matrix element
            call timing%stop(timing%mtxel_cor_tot)
          endif
          if (bExactlyZero) cycle
          if (peinf%flag_offdiag(ioff)) then
            ! FHJ: Store temporary arrays in other variables, and put in
            ! degeneracy factor due to star of q-points.
            jb = peinf%index_offdiag(ioff) + sig%ndiag
            achcor(jb,ispin) = achcor(jb,ispin) + ryd*neq(irq)*achtcor
            if (sig%freq_dep==2) then
              asxDyn(:,jb,ispin) = asxDyn(:,jb,ispin) + ryd*neq(irq)*asxtDyn(:)
              achDyn(:,jb,ispin) = achDyn(:,jb,ispin) + ryd*neq(irq)*achtDyn(:)
              achDyn_cor(:,jb,ispin) = achDyn_cor(:,jb,ispin) + ryd*neq(irq)*achtDyn_cor(:)
              achDyn_corb(:,jb,ispin) = achDyn_corb(:,jb,ispin) + ryd*neq(irq)*achtDyn_corb(:)
              ach2Dyn(:,jb,ispin) = ach2Dyn(:,jb,ispin) + ryd*neq(irq)*ach2tDyn(:)
              achD_n1(:,jb,ispin) = achD_n1(:,jb,ispin) + neq(irq)*achtD_n1(:)
            elseif (sig%freq_dep/=-1) then
              asig_imag(jb,ispin) = asig_imag(jb,ispin) + ryd*neq(irq)*asigt_imag
              asx(:,jb,ispin) = asx(:,jb,ispin) + ryd*neq(irq)*asxt(:)
              ach(:,jb,ispin) = ach(:,jb,ispin) + ryd*neq(irq)*acht(:)
              ach_n1(:,jb,ispin) = ach_n1(:,jb,ispin) + neq(irq)*acht_n1(:)
            endif
            if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
              achcor_n1(:,jb,ispin) = achcor_n1(:,jb,ispin) + neq(irq)*achtcor_n1(:)
            endif
          endif
!---------------------
! Compute offdiag bare exchange (SKIP this computation if you already know it)
          call timing%start(timing%bare_x)
          if ((.not. sig%use_xdat .and. xflag .and. (.not. sig%coul_mod_flag)) .or. (sig%coul_mod_flag &
            .and. (.not. sig%use_vxc2dat))) then
            axt=0.0d0
! XXX THREAD?
            do n1=1,peinf%nvband_node
              tempval = wfnkq%ekq(peinf%indext(n1),ispin) - sig%efermi
              if (tempval < sig%tol) then
                if(abs(tempval) < sig%tol) then
                  occ=0.5 ! Fermi-Dirac distribution = 1/2 at Fermi level
                else
                  occ=1D0
                endif
                do ig=1,ncoulb
                  axt=axt+aqsaug(ig,n1,sig%offmap(peinf%index_offdiag(ioff),1),ispin) &
                    *(aqsaug(ig,n1,sig%offmap(peinf%index_offdiag(ioff),2),ispin))*occ &
                    *vcoul(ig)
                enddo
              endif
            enddo
            if (peinf%flag_offdiag(ioff)) then
              ib = peinf%index_offdiag(ioff) + sig%ndiag
              ax(ib,ispin) = ax(ib,ispin) - neq(irq)*axt*ryd*sig%xfrac
            endif
          endif ! not using x.dat
          call timing%stop(timing%bare_x)
        enddo ! ioff (loop over bands for which we need offdiag Sigma)
        if (ispin.eq.sig%nspin) then
          if(allocated(vcoul))then;deallocate(vcoul);endif
          if (sig%noffdiag.gt.0) then
            if(allocated(aqsaug))then;deallocate(aqsaug);endif
          end if
          if ((sig%freq_dep.eq.0.or.sig%exact_ch.eq.1).and.irq_==nrq) then
            if(associated(aqsch))then;deallocate(aqsch);nullify(aqsch);endif
            if (nrq.gt.1.and.sig%noffdiag.gt.0) then
              if(associated(aqsaugcho))then;deallocate(aqsaugcho);nullify(aqsaugcho);endif
            end if
          endif
        endif
!!-------- End offdiag band loop -----------------------------------------------
      enddo ! ispin
!!-------- End spin loop -------------------------------------------------------
    enddo ! irq (loop over rq point in BZ summation)
    call progress_free(prog_info)
!!-------- End loop over k-points in irr bz with respect to kn (rq) ------------
    ! deallocate GPU stuff
    if ( mtxel_algo /= CPU_ALGO ) then
      call deallocate_acc_mtxel_sig ( )
    end if
!----------------------------
! Output unsymmetrized values of X,SX,CH
! Symmetrize X,SX,CH matrix elements over degenerate states
! Convert to eV and output symmetrized values of X,SX,CH
    if (peinf%inode.eq.0) then
      ! This beautiful piece of code should never be executed and here only because it somehow
      ! prevents gfortran from an incorrect optimization of this routine that will produce NaNs.
      ! It has some mysterious relation to the igp loop above reading from iunit_eps.
      ! Insanely, the presence of function 'isNaN' appears to be crucial here. --DAS
      if(ikn == -1) then
        call die("BUG: ikn = -1")
        write(0,*) isNaN(fact)
      endif
! DVF : sig%ncore_excl has to be substracted here because wfnk%ek/elda are defined in the
! read_wavefunction subroutine in input.f90 to be referenced to the case with
! no core states. Same goes for all subsequent occurences below.
      do ispin=1,sig%nspin
        write(6,989)ikn,sig%spin_index(ispin)
        if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
          write(6,*)
          ! ZL: either print Vxc or KIH
          if(.not. sig%use_kihdat) then ! standard Vxc output
            write(6,900) &
              "n", "Emf", "Eo", "Vxc", "X", "SX-X", "CH", "Cor", "Sig"
            do i=1,sig%ndiag
              write(6,901) sig%diag(i), wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), &
                wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin), dble(alda(i,ispin))*ryd, &
                dble(ax(i,ispin))*ryd, dble(asx(2,i,ispin))*ryd, &
                dble(ach(2,i,ispin)+achcor(i,ispin))*ryd, &
                dble(asx(2,i,ispin)+ach(2,i,ispin)+achcor(i,ispin))*ryd, &
                dble(ax(i,ispin)+asx(2,i,ispin)+ach(2,i,ispin)+achcor(i,ispin))*ryd
            enddo
          else ! KIH output
            write(6,900) &
              "n", "Emf", "Eo", "KIH", "X", "SX-X", "CH", "Cor", "Sig"
            do i=1,sig%ndiag
              write(6,901) sig%diag(i), wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), &
                wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin), dble(akih(i,ispin))*ryd, &
                dble(ax(i,ispin))*ryd, dble(asx(2,i,ispin))*ryd, &
                dble(ach(2,i,ispin)+achcor(i,ispin))*ryd, &
                dble(asx(2,i,ispin)+ach(2,i,ispin)+achcor(i,ispin))*ryd, &
                dble(ax(i,ispin)+asx(2,i,ispin)+ach(2,i,ispin)+achcor(i,ispin))*ryd
            enddo
          endif
        endif
        ! ZL: this is full frequency, EP does not support yet
        if (sig%freq_dep.eq.2) then
          write(6,*)
          if (sig%freq_dep_method==2) then
            ! ZL: Vxc vs. KIH
            if(.not. sig%use_kihdat) then
              write(6,900) &
                "n", "Emf", "Eo", "X", "Vxc", "Re Res", "Re Int", "Re Sig"
              write(6,900) &
                '', '', '', '', '', "Im Res", "Im Int", "Im Sig"
            else ! KIH
              write(6,900) &
                "n", "Emf", "Eo", "X", "KIH", "Re Res", "Re Int", "Re Sig"
              write(6,900) &
                '', '', '', '', '', "Im Res", "Im Int", "Im Sig"
            endif
          else
            ! ZL: Vxc vs. KIH
            if(.not. sig%use_kihdat) then
              write(6,900) &
                "n", "Emf", "Eo", "X", "Vxc", "Re SX-X", "Re CH", "Re Sig"
              write(6,900) &
                '', '', '', '', '', "Im SX-X", "Im CH", "Im Sig"
            else ! KIH
              write(6,900) &
                "n", "Emf", "Eo", "X", "KIH", "Re SX-X", "Re CH", "Re Sig"
              write(6,900) &
                '', '', '', '', '', "Im SX-X", "Im CH", "Im Sig"
            endif
          endif
          do i=1,sig%ndiag
! JRD: Find iw closest to e_lk
            diffmin = INF
            e_lk = wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin)
            ! FHJ: Figure out starting frequency for freq. grid
            if (sig%freq_grid_shift<2) then
              freq0 = sig%freqevalmin
            else
              freq0 = e_lk - sig%freqevalstep*(sig%nfreqeval-1)/2
            endif
            do iw=1,sig%nfreqeval
              diff = abs(freq0 + (iw-1)*sig%freqevalstep - e_lk)
              if (diff .lt. diffmin) then
                diffmin=diff
                iwlda=iw
              endif
            enddo
            ! ZL: print Vxc or KIH
            if(.not. sig%use_kihdat) then
              write(6,901) sig%diag(i), wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), &
                wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin), dble(ax(i,ispin))*ryd, &
                dble(alda(i,ispin))*ryd, dble(asxDyn(iwlda,i,ispin))*ryd, &
                dble(achDyn(iwlda,i,ispin)+achcor(i,ispin))*ryd, &
                dble(ax(i,ispin)+asxDyn(iwlda,i,ispin)+achDyn(iwlda,i,ispin)+achcor(i,ispin))*ryd, &
                dble(alda(i,ispin))*ryd
              write(6,902) aimag(asxDyn(iwlda,i,ispin))*ryd, aimag(achDyn(iwlda,i,ispin))*ryd, &
                aimag(ax(i,ispin)+asxDyn(iwlda,i,ispin)+achDyn(iwlda,i,ispin))*ryd, &
                aimag(ax(i,ispin)+asxDyn(iwlda,i,ispin)+ach2Dyn(iwlda,i,ispin))*ryd
            else ! KIH
              write(6,901) sig%diag(i), wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), &
                wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin), dble(ax(i,ispin))*ryd, &
                dble(akih(i,ispin))*ryd, dble(asxDyn(iwlda,i,ispin))*ryd, &
                dble(achDyn(iwlda,i,ispin)+achcor(i,ispin))*ryd, &
                dble(ax(i,ispin)+asxDyn(iwlda,i,ispin)+achDyn(iwlda,i,ispin)+achcor(i,ispin))*ryd, &
                dble(akih(i,ispin))*ryd
              ! ZL: turns out that following format 901, the last redundant akih will not be printed anyways
              ! but we keep it the same as the original code
              write(6,902) aimag(asxDyn(iwlda,i,ispin))*ryd, aimag(achDyn(iwlda,i,ispin))*ryd, &
                aimag(ax(i,ispin)+asxDyn(iwlda,i,ispin)+achDyn(iwlda,i,ispin))*ryd, &
                aimag(ax(i,ispin)+asxDyn(iwlda,i,ispin)+ach2Dyn(iwlda,i,ispin))*ryd
            endif
          enddo
        endif
      enddo
    endif ! node 0
989 format(1x,"Unsymmetrized values for ik =",i4,1x,"spin =",i2)
900 format(a6,8(a9))
901 format(i6,8(f9.3))
902 format(6x,4(9x),4(f9.3))
905 format(a6,8(a9),(a11))
906 format(i6,9(f9.3))
!----------------------------
! Symmetrize matrix elements
    if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
      call shiftenergy(sig, wfnk, alda, asx, ach, achcor, ach_n1, achcor_n1, &
        ax, efsto, asig, enew, zrenorm, nfreqgpp, sig%ncore_excl, akih) ! ZL: akih is optional argument
    endif
    if (sig%freq_dep.eq.2) then
      call shiftenergy_dyn(sig, wfnk, alda, asxDyn, achDyn, achDyn_cor, &
        achDyn_corb, ach2Dyn, achcor, achD_n1, achcor_n1, &
        ax, efstoDyn, asigDyn, enewDyn, enewDyn_nosr, &
        neqp1, neqp1_nosr, ikn, kp,sig%ncore_excl, akih) ! ZL: akih is optional argument
    endif
!----------------------------
! Write out matrix elements
    if (peinf%inode.eq.0) then
      write(6,'(a)')
      write(6,'(a)') ' Symmetrized values from band-averaging:'
      ! FHJ: write eqp0.dat and eqp1.dat to units 30 and 31.
      ! ZL: TODO - write separate files for EP
      write(30,'(3(f13.9),i8)') kp%rk(:,ikn), sig%nspin*sig%ndiag
      write(31,'(3(f13.9),i8)') kp%rk(:,ikn), sig%nspin*sig%ndiag
      if (sig%freq_dep==2) then
        do ispin = 1, sig%nspin
          do i = 1, sig%ndiag
            write(30,'(2i8,3f15.9)') ispin, sig%diag(i), &
              wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), efstoDyn(i,ispin)+achcor(i,ispin)
            write(31,'(2i8,3f15.9)') ispin, sig%diag(i), &
              wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), enewDyn(i,ispin)
          enddo
        enddo
      else
        do ispin = 1, sig%nspin
          do i = 1, sig%ndiag
            write(30,'(2i8,3f15.9)') ispin, sig%diag(i), &
              wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), efsto(i,ispin)+dble(achcor(i,ispin))
            write(31,'(2i8,3f15.9)') ispin, sig%diag(i), &
              wfnk%elda(sig%diag(i)-sig%ncore_excl,ispin), enew(i,ispin)+dble(achcor(i,ispin))*zrenorm(i,ispin)
          enddo
        enddo
      endif
      if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
        call write_result(kp, wfnk, sig, ax, asx, ach, &
          achcor, asig, alda, efsto, enew, zrenorm, ikn, sig%ncore_excl, akih)
        call write_result_hp(kp, wfnk, sig, ax, asx, ach, achcor, asig, alda, &
          efsto, enew, zrenorm, ikn, sig%ncore_excl, akih)
        if (sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
          do ispin=1,sig%nspin
            do i=1,sig%ndiag
              if (abs(efsto(i,ispin)+real(achcor(i,ispin),kind=DP) - wfnk%ek(sig%diag(i)-sig%ncore_excl,ispin)).gt.sig%dw) then
                eqp1_warns(1) = .true.
              endif
            enddo
          enddo
        endif
        call ch_converge(kp, sig, cmplx(ach_n1,kind=DPC), cmplx(achcor_n1,kind=DPC), ikn)
      endif
      if (sig%freq_dep==2) then
        call write_result_dyn(kp, wfnk, sig, ax, asxDyn, &
          achDyn, achDyn_corb, achcor, asigDyn, alda, efstoDyn, enewDyn, &
          enewDyn_nosr, neqp1, neqp1_nosr, ikn,sig%ncore_excl, akih)
        call write_result_dyn_hp(kp, wfnk, sig, ax, asxDyn, achDyn, &
          achDyn_corb, achcor, asigDyn, alda, efstoDyn, enewDyn, enewDyn_nosr, &
          neqp1, neqp1_nosr, ikn,sig%ncore_excl, akih)
        if (any(neqp1<0)) eqp1_warns(2) = .true. !FF extrap
        if (any(neqp1>1)) eqp1_warns(3) = .true. !FF mult solns
        if (any(neqp1==0)) eqp1_warns(4) = .true. !FF no solns
        call ch_converge(kp, sig, achD_n1, cmplx(achcor_n1,kind=DPC), ikn)
      endif
    endif
    if(peinf%inode == 0) then
!----------------------------
! If not using vxc.dat, create and write Vxc in it
      call timing%start(timing%vxc)
      if(.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX) then
        if(.not.sig%use_kihdat) then ! ZL: only if not using KIH
          call write_matrix_elements_type(120, kp%rk(:, ikn), sig, cmplx(alda(:,:),kind=DPC))
        endif
      endif
      call timing%stop(timing%vxc)
!----------------------------
! If not using x.dat, create and write X in it
      call timing%start(timing%bare_x)
      if (sig%coul_mod_flag .and. (.not. sig%use_vxc2dat)) then
        call write_matrix_elements_type(121, kp%rk(:, ikn), sig, cmplx(ax(:,:),kind=DPC))
      else if((.not. sig%use_xdat) .and. xflag .and. (.not. sig%coul_mod_flag)) then
        ax(:,:) = ax(:,:) / sig%xfrac
        call write_matrix_elements_type(119, kp%rk(:, ikn), sig, cmplx(ax(:,:),kind=DPC))
      endif
      call timing%stop(timing%bare_x)
    endif
    if(allocated(indrq))then;deallocate(indrq);endif
    if(allocated(neq))then;deallocate(neq);endif
    if(allocated(itnrq))then;deallocate(itnrq);endif
    if(allocated(rq))then;deallocate(rq);endif
    if(allocated(kg0))then;deallocate(kg0);endif
    if(associated(wfnk%zk))then;deallocate(wfnk%zk);nullify(wfnk%zk);endif
    if (peinf%inode==0) then
      !FHJ: Flush a couple of units so that the user can retrieve some data if
      !the calculation crashes
      if (sig%freq_dep==2 .or. (sig%fdf==-3 .and. sig%freq_dep==1)) then
        FLUSH(8000)
      endif
      FLUSH(0)
      FLUSH(6)
      FLUSH(8)
    endif
  enddo ! ika (loop over k-points sig%nkn)
  call dealloc_grid(gr)
  call destroy_qran() ! from vcoul_generator
!--------- End loop over rkn points (for which we calculate sigma) -------------
  if (peinf%inode==0) then
    do iunit=6,8,2
      write(iunit,'()')
      write(iunit,'(a)') repeat('=', 80)
      write(iunit,'()')
      write(iunit,'(1x,a)') '   n = band index.'
      write(iunit,'(1x,a)') ' Emf = "inner" mean-field energy eigenvalue used to construct Sigma(E),'
      write(iunit,'(1x,a)') '       read from WFN_inner.'
      write(iunit,'(1x,a)') '  Eo = "outer" mean-field energy eigenvalue where we center the evaluation'
      write(iunit,'(1x,a)') '       frequency grid {E} of Sigma(E). Defaults to Emf, unless'
      write(iunit,'(1x,a)') '       you use WFN_outer and eqp_outer.dat / scissors_outer.'
      if (sig%freq_dep==2.and.sig%freq_grid_shift/=2) then
        write(iunit,'(1x,a)') '       (Note: your freq. grid {E} does not contain Eo, so we actually report'
        write(iunit,'(1x,a)') '        below the self energy evaluated at the freq. Eo` closest to Eo.)'
      endif
      write(iunit,'(1x,a)') ' Vxc = exchange-correlation pot., calculated from VXC or read from vxc.dat.'
      ! ZL: for KIH
      write(iunit,'(1x,a)') ' KIH = Kinetic energy + Ionic potential + Hartree, read from kih.dat.'
      write(iunit,'(1x,a)') '   X = bare exchange.'
      if (sig%freq_dep==2.and.sig%freq_dep_method==2) then
        write(iunit,'(1x,a)') ' Res = residue contrib. to Sigma(E) at energy E=Eo.'
        write(iunit,'(1x,a)') ' Int = contrib. to Sigma(E) at energy E=Eo due to integral over imag. freqs.'
        write(iunit,'(1x,a)') ' Cor = Res + Int = correlation portion of Sigma(E) at energy E=Eo.'
        write(iunit,'(1x,a)') ' Sig = X + Cor = self energy Sigma(E) at energy E=Eo.'
      elseif(sig%freq_dep==2) then
        write(iunit,'(1x,a)') '  SX = screened exchange contrib. to Sigma(E) at energy E=Eo.'
        write(iunit,'(1x,a)') '  CH = Coulomb hole contrib. to Sigma(E) at energy E=Eo.'
        write(iunit,'(1x,a)') ' Cor = correlation portion of Sigma(E) at energy E=Eo'
        write(iunit,'(1x,a)') ' Sig = SX + CH = X + Cor = self energy, Sigma(E), at energy E=Eo.'
      else
        write(iunit,'(1x,a)') '  SX = screened exchange contrib. to Sigma(E) at energy E=Eo'
        write(iunit,'(1x,a)') '  CH = Coulomb hole contrib. to Sigma(E) at energy E=Eo'
        write(iunit,'(1x,a)') ' Cor = SX-X + CH = correlation portion of Sigma(E) at energy E=Eo.'
        write(iunit,'(1x,a)') ' Sig = X + Cor = self energy, Sigma(E), at energy E=Eo.'
      endif
      write(iunit,'(1x,a)') 'Eqp0 = on-shell QP energy = Emf - Vxc + Sig(Eo)'
      write(iunit,'(1x,a)') '       Eqp0 is *not* the recommended quantity to use for QP properties.'
      write(iunit,'(1x,a)') 'Eqp1 = off-shell solution to the linearized  Dyson`s equation'
      write(iunit,'(1x,a)') '     = Eqp0 + (dSig/dE) / (1 - dSig/dE) * (Eqp0 - Eo),'
      write(iunit,'(1x,a)') '       or a full linear interpolation if more freq. points where computed.'
      write(iunit,'(1x,a)') '       Eqp1 is the recommended quantity to use for QP properties.'
      if (sig%freq_dep==2) then
        write(iunit,'(1x,a)') 'Soln = Whether we found a solution to Dyson`s equation for Eqp1 (see notes).'
      endif
      write(iunit,'(1x,a)') ' Znk = quasiparticle renormalization factor'
      if (sig%exact_ch==1) then
        write(iunit,'()')
        write(iunit,'(1x,a)') 'Notes on the static remainder:'
        write(iunit,'(1x,a)') '- Unprimed values, such as Eqp0, are evaluated WITH the static remainder'
        write(iunit,'(1x,a)') '- Primed values, such as Eqp0`, are evaluated WITHOUT the static remainder'
      endif
      if (sig%freq_dep==2) then
        write(iunit,'()')
        write(iunit,'(1x,a)') 'Notes on the solutions to Dyson`s equation (Soln):'
        write(iunit,'(1x,a)') '- When we solve Dyson`s eqn. for Eqp1, we might find several or no roots.'
        write(iunit,'(1x,a)') '- If Soln=unique, there is only one solution and the answer is unambiguous'
        write(iunit,'(1x,a)') '- If Soln=extrap{-/+}, we found a solution by extrapolating Sigma(E) to values of'
        write(iunit,'(1x,a)') '  E smaller/greater than our frequency grid. PROCEED WITH CAUTION'
        write(iunit,'(1x,a)') '- If Soln=MULT:, there are multiple solutions to Dyson`s equation. We'
        write(iunit,'(1x,a)') '  report the solution closer to Eqp0. PROCEED WITH CAUTION'
        write(iunit,'(1x,a)') '- If Soln=NO_SOLN!, there is no solution and extrap. is not possible. This'
        write(iunit,'(1x,a)') '  can happen if your frequency grid for Sigma is too small, and Sigma(E) is'
        write(iunit,'(1x,a)') '  ill-behaved. We simply copy Eqp0 to Eqp1, so DON`T TRUST THE ANSWER.'
        write(iunit,'(1x,a)') '- The number of sols. in full-freq. calcs. is affected by self-consistently'
        write(iunit,'(1x,a)') '  updating the eigenvalues in WFN_inner. Consider reruning your calculation'
        write(iunit,'(1x,a)') '  with eqp or scissors corrections if you have many unstable solutions.'
      elseif (sig%freq_dep==1.or.sig%freq_dep==3) then
        write(iunit,'()')
        write(iunit,'(1x,a)') 'Notes on the finite_difference_form from sigma.inp file:'
        write(iunit,'(1x,a)') '  none    : -2 => dSig/dE = 0 (skip the expansion)'
        write(iunit,'(1x,a)') '  backward: -1 => dSig/dE = (Sig(Eo) - Sig(Eo-dE)) / dE'
        write(iunit,'(1x,a)') '  central :  0 => dSig/dE = (Sig(Eo+dE) - Sig(Eo-dE)) / (2*dE)'
        write(iunit,'(1x,a)') '  forward :  1 => dSig/dE = (Sig(Eo+dE) - Sig(Eo)) / dE'
        write(iunit,'(1x,a)') '  default :  2 => forward for diagonal and none for off-diagonal'
        write(iunit,'(1x,a)') '  dE is finite_difference_spacing from Sigma.inp file.'
        write(iunit,'(1x,a,i0,a,f0.3,a)') '  We are using the form #', sig%fdf ,' with dE = ', sig%dw, ' eV.'
      endif
      write(iunit,'()')
      write(iunit,'(1x,a)') 'General notes:'
      write(iunit,'(1x,a)') '- All energies are reported here in eV.'
      write(iunit,'(1x,a)') '- Both Emf and Vxc contain the average pot. Vxc0, so Vxc0 doesn`t affect Sigma.'
      write(iunit,'(1x,a)') '- Eqp1 and Eqp0 are Eqs. (36-37) from Hybertsen & Louie PRB 34 5390.'
      write(iunit,'(1x,a)') '- We recommend you use Eqp1 for QP properties of materials.'
      write(iunit,'()')
      write(iunit,'(a)') repeat('=', 80)
      write(iunit,'()')
    enddo
    if (sig%noffdiag/=0) then
      write(6,776)
      write(8,776)
    endif
    ! FHJ: Flush these units otherwise the warnings that follow will break the
    ! formating of the messages we want to write to unit 6
    FLUSH(0)
    FLUSH(6)
    FLUSH(8)
  endif
776 format( &
      4x,"n = band index of bra wavefunction", &
      /,4x,"m = band index of ket wavefunction", &
      /,4x,"l = band index of energy eigenvalue",/, &
      /,1x,"< psi_n(k) |      X      | psi_m(k) >" &
      /,1x,"< psi_n(k) | SX(Eo_l(k)) | psi_m(k) >" &
      /,1x,"< psi_n(k) | CH(Eo_l(k)) | psi_m(k) >" &
      /,1x,"< psi_n(k) |     Vxc     | psi_m(k) >" &
      /,1x,"< psi_n(k) | T + V_ion + V_H | psi_m(k) >",/)
  if (peinf%inode==0) then
    FLUSH(0)
    FLUSH(6)
    FLUSH(8)
    if (eqp1_warns(1)) then
      write(0,'(/a)') "WARNING: |Eqp0 - Eo| > finite_difference_spacing. Linear extrapolation for eqp1"
      write(0,'(a)') "may be inaccurate. You should test the validity of eqp1 by rerunning the"
      write(0,'(a)') "calculation with the self energy evaluated at the eqp0 energies. For that,"
      write(0,'(a)') "use the eqp_outer.dat file, created with eqp.py script and point WFN_outer to"
      write(0,'(a/)') "WFN_inner, if you were not already using WFN_outer."
    endif
    if (eqp1_warns(2)) then
      write(0,'(/a)') "WARNING: Some solutions to Dyson`s equation for Eqp1 were extrapolated."
      write(0,'(a/)') "         PROCEED WITH CAUTION!"
    endif
    if (eqp1_warns(3)) then
      write(0,'(/a)') "WARNING: There are multiple solutions to Dyson`s equation for some values of Eqp1,"
      write(0,'(a)') "         so we keep those closer to Eqp0."
      write(0,'(a/)') "         PROCEED WITH **EXTREME** CAUTION!"
    endif
    if (eqp1_warns(4)) then
      write(0,'(/a)') "WARNING: Could not find any solution to Dyson`s equation for some values of Eqp1,"
      write(0,'(a)') "         so we are using Eqp0 for those states."
      write(0,'(a/)') "         PROCEED WITH **EXTREME** CAUTION! SOME VALUES FOR EQP1 ARE JUST *WRONG*!"
    endif
  endif
  if (peinf%inode.eq.0) then
    if (imagvxcflag) write(0,677) 'Vxc'
    if (imagxflag) write(0,677) 'exchange'
  endif
677 format( &
      1x,"WARNING: ",a," diagonal matrix elements have large imaginary part.",/, &
      3x,"This may indicate a problem with your calculation.",/)
!--------- Clean Up and Finish -------------------------------------------------
  call destroy_fftw_plans()
  call kp%free()
  if (sig%spl_tck%n>0) then
    if(associated(sig%spl_tck%t))then;deallocate(sig%spl_tck%t);nullify(sig%spl_tck%t);endif
    if(associated(sig%spl_tck%c))then;deallocate(sig%spl_tck%c);nullify(sig%spl_tck%c);endif
  endif
  if (sig%spl_tck_outer%n>0) then
    if(associated(sig%spl_tck_outer%t))then;deallocate(sig%spl_tck_outer%t);nullify(sig%spl_tck_outer%t);endif
    if(associated(sig%spl_tck_outer%c))then;deallocate(sig%spl_tck_outer%c);nullify(sig%spl_tck_outer%c);endif
  endif
  if(associated(wfnk%isrtk))then;deallocate(wfnk%isrtk);nullify(wfnk%isrtk);endif
  if(associated(wfnk%ek))then;deallocate(wfnk%ek);nullify(wfnk%ek);endif
  if(associated(wfnk%elda))then;deallocate(wfnk%elda);nullify(wfnk%elda);endif
  if(associated(sig%qpt))then;deallocate(sig%qpt);nullify(sig%qpt);endif
  if (sig%ndiag.ne.0) then
    if(associated(sig%diag))then;deallocate(sig%diag);nullify(sig%diag);endif
  end if
  if (sig%noffdiag.gt.0) then
    if(associated(sig%off1))then;deallocate(sig%off1);nullify(sig%off1);endif
    if(associated(sig%off2))then;deallocate(sig%off2);nullify(sig%off2);endif
    if(associated(sig%off3))then;deallocate(sig%off3);nullify(sig%off3);endif
    if(associated(sig%offmap))then;deallocate(sig%offmap);nullify(sig%offmap);endif
    if(sig%elph) then
      if(associated(sig%off_ep))then;deallocate(sig%off_ep);nullify(sig%off_ep);endif
      if(associated(sig%offmap_ep))then;deallocate(sig%offmap_ep);nullify(sig%offmap_ep);endif
    endif
  endif
  if(associated(sig%indkn))then;deallocate(sig%indkn);nullify(sig%indkn);endif
  if(.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX) then
    if(associated(sig%vxc))then;deallocate(sig%vxc);nullify(sig%vxc);endif
  end if
  if(associated(gvec%components))then;deallocate(gvec%components);nullify(gvec%components);endif
  if(associated(gvec%index_vec))then;deallocate(gvec%index_vec);nullify(gvec%index_vec);endif
  if(sig%freq_dep.eq.1) then
    if(associated(wpg%rho))then;deallocate(wpg%rho);nullify(wpg%rho);endif
  endif
  if(associated(wfnkq%isrtkq))then;deallocate(wfnkq%isrtkq);nullify(wfnkq%isrtkq);endif
  if(associated(wfnkq%ekq))then;deallocate(wfnkq%ekq);nullify(wfnkq%ekq);endif
  if(associated(peinf%indext))then;deallocate(peinf%indext);nullify(peinf%indext);endif
  if(associated(peinf%indext_dist))then;deallocate(peinf%indext_dist);nullify(peinf%indext_dist);endif
  if(associated(peinf%ntband_dist))then;deallocate(peinf%ntband_dist);nullify(peinf%ntband_dist);endif
  if(allocated(alda))then;deallocate(alda);endif
  if(allocated(akih))then;deallocate(akih);endif ! ZL: for KIH
  if(allocated(ax))then;deallocate(ax);endif
  if(allocated(achcor))then;deallocate(achcor);endif
  if(allocated(asig_imag))then;deallocate(asig_imag);endif
  if (sig%freq_dep/=0 .and. sig%exact_ch==1) then
    if(associated(achtcor_n1))then;deallocate(achtcor_n1);nullify(achtcor_n1);endif
  endif
  if(allocated(achcor_n1))then;deallocate(achcor_n1);endif
  if (sig%freq_dep.eq.-1.or.sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
    if(allocated(asx))then;deallocate(asx);endif
    if(allocated(ach))then;deallocate(ach);endif
    if(allocated(asig))then;deallocate(asig);endif
    if(associated(acht_n1))then;deallocate(acht_n1);nullify(acht_n1);endif
    if(allocated(ach_n1))then;deallocate(ach_n1);endif
    if(allocated(enew))then;deallocate(enew);endif
    if(allocated(efsto))then;deallocate(efsto);endif
  endif
  if (sig%freq_dep.eq.2) then
    ! ZL: EP not implemented for this frequency dependence
    if(allocated(asxDyn))then;deallocate(asxDyn);endif
    if(allocated(achDyn))then;deallocate(achDyn);endif
    if(allocated(achDyn_cor))then;deallocate(achDyn_cor);endif
    if(allocated(achDyn_corb))then;deallocate(achDyn_corb);endif
    if(allocated(ach2Dyn))then;deallocate(ach2Dyn);endif
    if(allocated(asigDyn))then;deallocate(asigDyn);endif
    if(associated(achtD_n1))then;deallocate(achtD_n1);nullify(achtD_n1);endif
    if(allocated(achD_n1))then;deallocate(achD_n1);endif
    if(allocated(efstoDyn))then;deallocate(efstoDyn);endif
    if(allocated(enewDyn))then;deallocate(enewDyn);endif
    if(allocated(enewDyn_nosr))then;deallocate(enewDyn_nosr);endif
    if(allocated(neqp1))then;deallocate(neqp1);endif
    if(allocated(neqp1_nosr))then;deallocate(neqp1_nosr);endif
    if(associated(asxtDyn))then;deallocate(asxtDyn);nullify(asxtDyn);endif
    if(associated(achtDyn))then;deallocate(achtDyn);nullify(achtDyn);endif
    if(associated(achtDyn_cor))then;deallocate(achtDyn_cor);nullify(achtDyn_cor);endif
    if(associated(achtDyn_corb))then;deallocate(achtDyn_corb);nullify(achtDyn_corb);endif
    if(associated(ach2tDyn))then;deallocate(ach2tDyn);nullify(ach2tDyn);endif
  endif
  if (sig%freq_dep.eq.1.or.sig%freq_dep.eq.3) then
    if(associated(acht))then;deallocate(acht);nullify(acht);endif
    if(associated(asxt))then;deallocate(asxt);nullify(asxt);endif
  endif
  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
    if(associated(eps))then;deallocate(eps);nullify(eps);endif
  endif
  if (sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
    if(associated(epsR))then;deallocate(epsR);nullify(epsR);endif
    if (sig%need_advanced) then
      if(associated(epsA))then;deallocate(epsA);nullify(epsA);endif
    endif
  endif
  if(allocated(isrtrq))then;deallocate(isrtrq);endif
  if (sig%freq_dep.eq.0.or.sig%exact_ch.eq.1) then
    if(associated(isrtrqi))then;deallocate(isrtrqi);nullify(isrtrqi);endif
  endif
  if(allocated(ekin))then;deallocate(ekin);endif
  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
    if(allocated(isrtq))then;deallocate(isrtq);endif
    if(allocated(isrtqi))then;deallocate(isrtqi);endif
    if(allocated(ind))then;deallocate(ind);endif
    if(allocated(indinv))then;deallocate(indinv);endif
    if(allocated(ph))then;deallocate(ph);endif
    if (peinf%inode .eq. 0) then
      if(associated(epsmpi%isrtq))then;deallocate(epsmpi%isrtq);nullify(epsmpi%isrtq);endif
      if(associated(epsmpi%isrtqi))then;deallocate(epsmpi%isrtqi);nullify(epsmpi%isrtqi);endif
    endif
    if(associated(epsmpi%qk))then;deallocate(epsmpi%qk);nullify(epsmpi%qk);endif
    if(associated(epsmpi%nmtx))then;deallocate(epsmpi%nmtx);nullify(epsmpi%nmtx);endif
    if(sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
      if(associated(epsmpi%eps))then;deallocate(epsmpi%eps);nullify(epsmpi%eps);endif
    else
      if(associated(epsmpi%epsR))then;deallocate(epsmpi%epsR);nullify(epsmpi%epsR);endif
      if (sig%need_advanced) then
        if(associated(epsmpi%epsA))then;deallocate(epsmpi%epsA);nullify(epsmpi%epsA);endif
      endif
      if (sig%do_sigma_subspace) then
       if(associated(sig%epssub%eps_sub_info))then;deallocate(sig%epssub%eps_sub_info);nullify(sig%epssub%eps_sub_info);endif
       if(associated(sig%epssub%eigenvec_sub))then;deallocate(sig%epssub%eigenvec_sub);nullify(sig%epssub%eigenvec_sub);endif
       if(associated(sig%epssub%eps_sub))then;deallocate(sig%epssub%eps_sub);nullify(sig%epssub%eps_sub);endif
       if(associated(sig%epssub%eps_wings_rows))then;deallocate(sig%epssub%eps_wings_rows);nullify(sig%epssub%eps_wings_rows);endif
       if(associated(sig%epssub%eps_wings_cols))then;deallocate(sig%epssub%eps_wings_cols);nullify(sig%epssub%eps_wings_cols);endif
       !XXX if(associated(sig%epssub%eps_wings_correction_rows))then;deallocate(sig%epssub%eps_wings_correction_rows);nullify(sig%epssub%eps_wings_correction_rows);endif
       !XXX if(associated(sig%epssub%eps_wings_correction_cols))then;deallocate(sig%epssub%eps_wings_correction_cols);nullify(sig%epssub%eps_wings_correction_cols);endif
       !MDB name too long, workaround
       if(associated(sig%epssub%eps_wings_correction_rows))then
         deallocate(sig%epssub%eps_wings_correction_rows)
         nullify(sig%epssub%eps_wings_correction_rows)
       endif
       if(associated(sig%epssub%eps_wings_correction_cols))then
         deallocate(sig%epssub%eps_wings_correction_cols)
         nullify(sig%epssub%eps_wings_correction_cols)
       endif
       if(associated(sig%epssub%vcoul_sub))then;deallocate(sig%epssub%vcoul_sub);nullify(sig%epssub%vcoul_sub);endif
      end if
    endif
  endif
  if(associated(wfnkqmpi%nkptotal))then;deallocate(wfnkqmpi%nkptotal);nullify(wfnkqmpi%nkptotal);endif
  if(associated(wfnkqmpi%isort))then;deallocate(wfnkqmpi%isort);nullify(wfnkqmpi%isort);endif
  if(associated(wfnkqmpi%el))then;deallocate(wfnkqmpi%el);nullify(wfnkqmpi%el);endif
  if(associated(wfnkqmpi%qk))then;deallocate(wfnkqmpi%qk);nullify(wfnkqmpi%qk);endif
  if(associated(wfnkqmpi%band_index))then;deallocate(wfnkqmpi%band_index);nullify(wfnkqmpi%band_index);endif
  if(associated(wfnkqmpi%cg))then;deallocate(wfnkqmpi%cg);nullify(wfnkqmpi%cg);endif
  if(sig%nkn.gt.1) then
    if(associated(wfnkmpi%nkptotal))then;deallocate(wfnkmpi%nkptotal);nullify(wfnkmpi%nkptotal);endif
    if(associated(wfnkmpi%isort))then;deallocate(wfnkmpi%isort);nullify(wfnkmpi%isort);endif
    if(associated(wfnkmpi%qk))then;deallocate(wfnkmpi%qk);nullify(wfnkmpi%qk);endif
    if(associated(wfnkmpi%el))then;deallocate(wfnkmpi%el);nullify(wfnkmpi%el);endif
    if(associated(wfnkmpi%elda))then;deallocate(wfnkmpi%elda);nullify(wfnkmpi%elda);endif
    if(associated(wfnkmpi%cg))then;deallocate(wfnkmpi%cg);nullify(wfnkmpi%cg);endif
  endif
  if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
    if(associated(epsmpi%inv_igp_index))then;deallocate(epsmpi%inv_igp_index);nullify(epsmpi%inv_igp_index);endif
  endif
  call show_references()
!----------------------------
! Time Accounting
  call timing%stop(timing%total)
  call timing%print(common_timing)
!----------------------------
! Close files and finish
  call close_file(55) ! file sigma.inp
  if(peinf%inode == 0) then
    call close_file(8) ! file sigma_hp.log
    call close_file(30)
    call close_file(31)
    if (sig%coul_mod_flag .and. (.not. sig%use_vxc2dat)) then
      call close_file(121) ! file vxc2.dat
    elseif ((.not. sig%use_xdat) .and. xflag .and. (.not. sig%coul_mod_flag)) then
      call close_file(119) ! file x.dat
    endif
    if(.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX) then
      if(.not.sig%use_kihdat) then ! ZL: only if not using KIH
        call close_file(120) ! file vxc.dat
      endif
    end if
    if (.not.(sig%freq_dep .eq. 0 .and. sig%exact_ch .eq. 1) .and. .not. (sig%freq_dep == -1)) then
      call close_file(127) ! file ch_converge.dat
    endif
    if (sig%iwritecoul .eq. 1) then
      call close_file(19) ! file vcoul
    endif
  endif
  if (peinf%inode.eq.0 .and. (sig%freq_dep.eq.2 .or. (sig%fdf.eq.-3 .and. sig%freq_dep.eq.1))) then
    write(8000, '(/,a)')'# Please refer to Sigma/README for more information about this file.'
    call close_file(8000) ! file spectrum.dat
  endif
  call write_memory_usage()
end program sigma
