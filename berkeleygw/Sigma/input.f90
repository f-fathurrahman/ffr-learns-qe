!===============================================================================
!
! Routines:
!
! (1) input() Originally By ? Last Modified 10/5/2009 (gsm)
!
! Sets up various data structures. Reads in and distributes wavefunctions.
! Stores wavefunctions in structures distributed in memory.
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
module input_m
  use global_m
  use checkbz_m
  use createpools_m
  use eqpcor_m
  use fftw_m
  use fullbz_m
  use input_utils_m
  use misc_m
  use read_rho_vxc_m
  use scissors_m
  use sort_m
  use wfn_rho_vxc_io_m
  use io_utils_m
  use epsread_hdf5_m
  use timing_m, only: timing => sigma_timing
  use inread_m
  use genwf_mpi_m
  use wfn_io_hdf5_m
  implicit none
  private
  public :: input
contains
subroutine input(crys,gvec,syms,kp,wpg,sig,wfnk,iunit_c,iunit_k,fnc,fnk,wfnkqmpi,wfnkmpi, wfnk_phonq, wfnk_phonq_mpi, ep_read_in)
  type (crystal), intent(out) :: crys
  type (gspace), intent(out) :: gvec
  type (symmetry), intent(out) :: syms
  type (kpoints), intent(out) :: kp
  type (wpgen), intent(out) :: wpg
  type (siginfo), intent(inout) :: sig ! ZL modifies
  type (wfnkstates), intent(out) :: wfnk, wfnk_phonq ! ZL adds wfnk_phonq
  integer, intent(in) :: iunit_c, iunit_k
  character*20, intent(in) :: fnc, fnk
  type (wfnkqmpiinfo), intent(out) :: wfnkqmpi
  type (wfnkmpiinfo), intent(out) :: wfnkmpi, wfnk_phonq_mpi ! ZL adds wfnk_phonq_mpi
  character :: fncor*32, tmpfn*16
  character :: tmpstr*100,tmpstr1*16,tmpstr2*16
  character :: filenameh5*80
  integer :: i,ierr,ii,j,jj,k,kk
  integer :: ikn, irk, nq
  integer :: ncoul,ncoulb,ncouls,nmtx,nmtx_l,nmtx_col
  integer, allocatable :: isrt(:), isrti(:)
  integer :: ipe
  integer :: Nrod,Nplane,Nfft(3),dNfft(3),dkmax(3),nmpinode
  integer :: ntband,npools,ndiag_max,noffdiag_max,ntband_max
  real(DP) :: scalarsize,qk(3),vcell,qtot,omega_plasma,rhog0
  real(DP) :: mem,fmem,rmem,rmem2,scale,dscale
  integer, allocatable :: gvecktmp(:,:)
  logical :: do_ch_sum, skip_checkbz
  type(grid) :: gr
  real(DP) :: del
  integer :: g0(3), i_looper ! ZL: to map EP sig%indk_phonq_g0
  character(len=3) :: sheader
  integer :: iflavor, freq_dep, error
  integer :: ng_eps, nq_eps, nmtxmax_eps, qgrid_eps(3), nfreq_eps, nfreq_imag_eps
  real(DP) :: ecuts
  logical :: file_exists
  ! ZL: electron phonon (EP) control
  logical, intent(in), optional :: ep_read_in
  logical :: ep_read
  logical :: do_phonq ! ZL: if do_phonq = .false., wfnk_phonq and wfnk_phonq_mpi will not
                     ! be used, nor allocated or evaluated
 
  ! ZL: set up EP
  do_phonq = .false.
  ep_read = .false.
!---------------
! Read sigma.inp
! write(6,*) 'before inread'
  if (.not.ep_read) then
    call inread(sig) ! ZL: in inread(), if nphonq .eq. 1, program exits.
                     ! ZL: so after inread(), simply use sig%elph
  else
    if (peinf%inode.eq.0) write(6,*) 'Skip inread sigma for dWFN.'
  endif
! write(6,*) 'after inread'
  ! ZL: conditions for generating wfnk_phonq, wfnk_phonq_mpi
!------------------------
! Initialize HDF5
  if(sig%nkn .lt. 1) then
    if(peinf%inode.eq.0) write(0,*) 'sig%nkn < 1'
    call die('sig%nkn')
  endif
!---------------------------------
! Determine the available memory
  call procmem(mem,nmpinode)
  if(peinf%inode.eq.0) then
    write(6,998) mem/1024.0d0**2
  endif
998 format(1x,'Memory available: ',f0.1,' MB per PE')
  fmem=mem/8.0d0
!-------------------------------
! (gsm) Estimate the required memory
  if(peinf%inode.eq.0) then
! determine number of frequencies and number of q-points
    if (sig%freq_dep.eq.-1) then
      nq=sig%nq
      sig%nFreq=1 ! ZL: for EP, this is fine reevaluate
    endif
    if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1.or.sig%freq_dep.eq.2.or.sig%freq_dep.eq.3) then
        call open_file(10,file='eps0mat',form='unformatted',status='old')
        read(10)
        read(10) i, sig%nFreq
        read(10)
        read(10)
        read(10)
        read(10)
        read(10) ecuts
        read(10) nq
        call close_file(10)
        call open_file(11,file='epsmat',form='unformatted',status='old',iostat=ierr)
        if (ierr.eq.0) then
          read(11)
          read(11)
          read(11)
          read(11)
          read(11)
          read(11)
          read(11)
          read(11) nq_eps
          nq = nq + nq_eps
          call close_file(11)
        endif
    endif
  endif
! FHJ: Read header of WFN file (and WFNq, if requred), perform consistency
! checks, and figure out the number of bands and ncrit.
  sheader = 'WFN'
  iflavor = 0
  if(sig%wfn_hdf5) then
  else
    if(peinf%inode == 0) then
      if(.not.ep_read) then
        write(6,'(a)') ' Reading header of WFN_inner'
        call open_file(25,file='WFN_inner',form='unformatted',status='old')
      else
        ! ZL: open dWFN, the perturbed wavefunction
        write(6,'(a)') ' Reading header of dWFN'
        call open_file(25,file='dWFN',form='unformatted',status='old')
      endif
    endif
    call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys)
  endif
  call check_trunc_kpts(sig%icutv, kp) ! This does not affect EP
  if (sig%ecutb<TOL_ZERO) then
    sig%ecutb = kp%ecutwfc
  endif
  if (sig%freq_dep/=-1 .and. sig%ecuts<TOL_ZERO) then
    sig%ecuts = ecuts
  endif
  if (sig%ntband==0 .and. (.not.ep_read)) sig%ntband = kp%mnband - 1
! DVF: Here we re-define the number of bands in the dynamical self energy
! summation to exclude the deep core states. This ensures that the pools will be
! set up properly. This achieves most of what is needed to exclude core states.
! There is a little more work in the i/o routines to make sure you`re getting the
! higher valence states.
! ZL: we disable this operation if doing elph to avoid possible mismaps of band indices
  if ((sig%ncore_excl.ne.0) .and. sig%elph) then
    call die('Electron-Phonon calculation is not consistent with sig%ncore_excl.ne.0 .')
  endif
  sig%ntband=sig%ntband-sig%ncore_excl
  allocate(kp%elda (sig%ntband,kp%nrk,kp%nspin))
  kp%el(:,:,:) = kp%el(:,:,:) - sig%avgpot / ryd
! DVF : add sig%ncore_excl to make sure we are getting the right energies. kp%el is setup in
! wfn_rho_vxc_io.f90 routine (see .p.f for something comprehensible) and includes core states.
  kp%elda(1:sig%ntband, 1:kp%nrk, 1:kp%nspin)=ryd*kp%el(1+sig%ncore_excl:sig%ntband+sig%ncore_excl, 1:kp%nrk, 1:kp%nspin)
  call scissors_shift(kp, sig%scis, sig%spl_tck)
  ! If QP corrections requested, read the corrected QP energies from file (in eV)
  if(sig%eqp_corrections .and. (.not.ep_read)) then
    fncor='eqp.dat'
! DVF : add sig%ncore_excl to make sure we are getting the right energies. The eqp.dat file
! will have band indices referenced to the total number of states, including the core states.
    call eqpcor(fncor,peinf%inode,peinf%npes,kp,1+sig%ncore_excl,sig%ntband+sig%ncore_excl,0,0,kp%el,kp%el,kp%el,1,0)
    ! note: want in Ry since conversion occurs below
  endif
  if (peinf%inode==0 .and. (.not.ep_read)) then
    call find_efermi(sig%rfermi, sig%efermi, sig%efermi_input, kp, sig%ntband+sig%ncore_excl, 1+sig%ncore_excl, &
      "unshifted grid", should_search=.true., should_update=.true., write7=.false.)
  endif
  if (sig%nvband==0 .and. (.not.ep_read)) then
    sig%nvband = minval(kp%ifmax)
    sig%ncrit = maxval(kp%ifmax) - minval(kp%ifmax)
  endif
! DVF: Here we re-define the number of bands in the bare exchange
! summation to exclude the deep core states.
! ZL: this has been disabled already for EP
  sig%nvband=sig%nvband-sig%ncore_excl
! DVF: again, this is to exclude the core states.
  if ((sig%freq_dep==0 .and. sig%exact_ch==1) .or. sig%freq_dep==-1) then
    if (sig%ntband + sig%ncore_excl > max(maxval(sig%diag(1:sig%ndiag)), sig%nvband + sig%ncore_excl + sig%ncrit)) then
      ! Do not waste time and memory with the unoccupied bands which will not be used
      sig%ntband = max(maxval(sig%diag(1:sig%ndiag))-sig%ncore_excl, sig%nvband + sig%ncrit)
      if (peinf%inode==0 .and. sig%freq_dep==0) then
        write(0,'(a)') "WARNING: Resetting number_bands to number of bands actually needed for exact static CH."
      else if (peinf%inode==0) then
        write(0,'(a)') "WARNING: Resetting number_bands to number of bands actually needed for exchange."
      endif
    endif
  endif
  if(peinf%inode == 0) then
    if(.not.ep_read) then
      call logit('Reading WFN_inner -- gvec info')
    else
      call logit('Reading dWFN -- gvec info')
    endif
  endif
  allocate(gvec%components (3, gvec%ng))
  if(sig%wfn_hdf5) then
  else
      call read_binary_gvectors(25, gvec%ng, gvec%ng, gvec%components)
  endif
!-----------------------
! sort G-vectors with respect to their kinetic energy
  allocate(isrt (gvec%ng))
  allocate(isrti (gvec%ng))
  allocate(gvec%ekin (gvec%ng))
  call kinetic_energies(gvec, crys%bdot, gvec%ekin)
  call sortrx(gvec%ng, gvec%ekin, isrt, gvec = gvec%components)
  ncouls = gcutoff(gvec%ng, gvec%ekin, isrt, sig%ecuts)
  ncoulb = gcutoff(gvec%ng, gvec%ekin, isrt, sig%ecutb)
  if(associated(gvec%ekin))then;deallocate(gvec%ekin);nullify(gvec%ekin);endif
  if (peinf%inode==0) then
    write(6,*)
    write(6,'(a)') ' Calculation parameters:'
    write(6,'(1x,a,f0.2)') '- Cutoff of the bare Coulomb interaction (Ry): ', sig%ecutb
    write(6,'(1x,a,f0.2)') '- Cutoff of the screened Coulomb interaction (Ry): ', sig%ecuts
    write(6,'(1x,a,i0)') '- Number of G-vectors up to the bare int. cutoff: ', ncoulb
    write(6,'(1x,a,i0)') '- Number of G-vectors up to the screened int. cutoff: ', ncouls
    write(6,'(1x,a,i0)') '- Total number of bands in the calculation: ', sig%ntband
    write(6,'(1x,a,i0)') '- Number of fully occupied valence bands: ', sig%nvband
    write(6,'(1x,a,i0)') '- Number of partially occ. conduction bands: ', sig%ncrit
    write(6,*)
  endif
  if(sig%ecuts > sig%ecutb) then
    call die("The screened_coulomb_cutoff cannot be greater than the bare_coulomb_cutoff.", &
      only_root_writes = .true.)
  endif
  if(sig%ntband > kp%mnband) then
    if(.not.ep_read) then
      call die("number_bands is larger than are available in WFN_inner", only_root_writes = .true.)
    else
      call die("number_bands is larger than are available in dWFN", only_root_writes = .true.)
    endif
  endif
  do ii = 1, sig%ndiag
    if (sig%diag(ii)>sig%ntband) then
      write(0,'(a,3(i0,1x))') 'ERROR: the band_index_max cannot be greater than number_bands: ', &
        ii, sig%diag(ii), sig%ntband
      call die('band_index_max cannot be greater than number_bands', only_root_writes = .true.)
    endif
  enddo
! FHJ: do vanilla WFN reading stuff
  allocate(gvecktmp (3,gvec%ng))
  gvecktmp(:,:)=gvec%components(:,:)
  do i=1,gvec%ng
    gvec%components(:,i)=gvecktmp(:,isrt(i))
    isrti(isrt(i)) = i
  enddo
  if(allocated(gvecktmp))then;deallocate(gvecktmp);endif
  call gvec_index(gvec)
!------------------------------------------------------------------------------
! Pools, distribution and memory estimation
!------------------------------------------------------------------------------
  if(peinf%inode == 0) then
    ncoul=max(ncouls,ncoulb)
    nmtx=ncouls
! divide bands over processors (this is repeated below)
    ntband=min(sig%ntband,kp%mnband)
    if (peinf%npools .le. 0 .or. peinf%npools .gt. peinf%npes) then
      call createpools(sig%ndiag,ntband,peinf%npes,npools,ndiag_max,ntband_max)
    else
      npools = peinf%npools
      if (mod(sig%ndiag,npools).eq.0) then
        ndiag_max=sig%ndiag/npools
      else
        ndiag_max=sig%ndiag/npools+1
      endif
      if (mod(ntband,peinf%npes/npools).eq.0) then
        ntband_max=ntband/(peinf%npes/npools)
      else
        ntband_max=ntband/(peinf%npes/npools)+1
      endif
    endif
    if (sig%noffdiag.gt.0) then
      if (mod(sig%noffdiag,npools).eq.0) then
        noffdiag_max=sig%noffdiag/npools
      else
        noffdiag_max=sig%noffdiag/npools+1
      endif
    endif
    nmtx_l=int(sqrt(dble(nmtx)**2/dble(peinf%npes/npools)))
    nmtx_col=int(dble(nmtx)/dble(peinf%npes/npools))+1
    scalarsize = sizeof_scalar()
! required memory
    rmem=0.0d0
! arrays eps and epstemp in program main
    if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
      rmem=rmem+(dble(nmtx_l)**2+dble(nmtx))*scalarsize
      rmem=rmem+(dble(nmtx_col)*nmtx)*scalarsize
    endif
! wtilde_array, imagloop_ig, imaginary_flag in mtxel_cor
    if (sig%freq_dep.eq.1) then
      rmem=rmem+(dble(nmtx_col)*nmtx)*16
      rmem=rmem+(dble(nmtx_col)*nmtx)*4
      rmem=rmem+(dble(nmtx_col)*nmtx)*1
    endif
! arrays epsR, epsRtemp, epsA, and epsAtemp in program main
    if (sig%freq_dep.eq.2 .or.sig%freq_dep.eq.3) then
      rmem=rmem+(dble(nmtx_l)**2+dble(nmtx)) &
        *dble(sig%nFreq)*2.0d0*scalarsize
    endif
! array epsmpi%eps in subroutine epscopy
    if (sig%freq_dep.eq.0.or.sig%freq_dep.eq.1) then
      rmem=rmem+(dble(nmtx_col)*nmtx)*dble(nq)*scalarsize
    endif
! arrays epsmpi%epsR and epsmpi%epsA in subroutine epscopy
    if (sig%freq_dep.eq.2 .or.sig%freq_dep.eq.3) then
      rmem=rmem+(nmtx_col*nmtx)*dble(nq) &
        *dble(sig%nFreq)*2.0d0*scalarsize
    endif
! array aqs in program main
    rmem=rmem+dble(ntband_max)*dble(ncoul)*scalarsize
! array aqsaug in program main
    if (sig%noffdiag.gt.0) then
      rmem=rmem+dble(ntband_max)*dble(ncoul) &
        *dble(sig%ndiag)*dble(sig%nspin)*scalarsize
    endif
! array aqsch in program main
    if (sig%exact_ch.eq.1) then
      rmem=rmem+dble(ncoul)*scalarsize
    endif
! arrays aqsaugchd and aqsaugcho in program main
    if (sig%exact_ch.eq.1.and.nq.gt.1) then
      rmem=rmem+dble(ncoul)*dble(ndiag_max) &
        *dble(sig%nspin)*scalarsize
      if (sig%noffdiag.gt.0) then
        rmem=rmem+dble(ncoul)*dble(noffdiag_max)* &
          dble(sig%nspin)*scalarsize
      endif
    endif
! array wfnk%zk in subroutine input
    rmem=rmem+dble(ndiag_max)*dble(kp%ngkmax)* &
      dble(sig%nspin*kp%nspinor)*scalarsize
! array wfnkq%zkq in subroutine genwf
    rmem=rmem+dble(ntband_max)*dble(kp%ngkmax) &
      *dble(sig%nspin*kp%nspinor)*scalarsize
! arrays zc in subroutine input and zin in subroutine genwf
    rmem=rmem+dble(kp%ngkmax)*dble(kp%nspin*kp%nspinor)*scalarsize
! array wfnkoff%zk in program main and subroutines mtxel_vxc and mtxel_cor
    if ((sig%exact_ch.eq.1.and.sig%noffdiag.gt.0).or. &
      (.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX)) then
      rmem=rmem+2.0d0*dble(kp%ngkmax)*scalarsize
    endif
! array gvec%index_vec in input
    rmem=rmem+dble(gvec%nFFTgridpts)*4.0d0
! arrays fftbox1 and fftbox2 in subroutines mtxel and mtxel_occ
    call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)
    rmem=rmem+dble(Nfft(1)*Nfft(2)*Nfft(3))*32.0d0
! arrays wfnkqmpi%isort and wfnkqmpi%cg in subroutine input
    rmem=rmem+dble(kp%ngkmax)*dble(kp%nrk)*4.0d0+ &
      dble(kp%ngkmax)*dble(ntband_max)*dble(sig%nspin*kp%nspinor)* &
      dble(kp%nrk)*dble(scalarsize)
! arrays wfnkmpi%isort and wfnkmpi%cg in subroutine input
    rmem=rmem+dble(kp%ngkmax)*dble(sig%nkn)*4.0d0+ &
      dble((ndiag_max*kp%ngkmax)/(peinf%npes/npools))* &
      dble(sig%nspin*kp%nspinor)*dble(sig%nkn)* &
      dble(scalarsize)
989 format(1x,'Memory required for execution: ',f0.1,' MB per PE')
    write(6,989) rmem/1024.0d0**2
! random numbers
    rmem=0.0D0
    if (sig%icutv/=TRUNC_BOX) then
! arrays ran, qran, and qran2
! (ran is deallocated before qran2 is allocated)
      rmem=rmem+6.0D0*dble(nmc)*8.0D0
    endif
! various truncation schemes
    rmem2=0.0d0
! cell wire truncation
    if (sig%icutv==TRUNC_WIRE) then
      dkmax(1) = gvec%FFTgrid(1) * n_in_wire
      dkmax(2) = gvec%FFTgrid(2) * n_in_wire
      dkmax(3) = 1
      call setup_FFT_sizes(dkmax,dNfft,dscale)
! array fftbox_2D
      rmem2=rmem2+dble(dNfft(1))*dble(dNfft(2))*16.0d0
! array inv_indx
      rmem2=rmem2+dble(Nfft(1))*dble(Nfft(2))*dble(Nfft(3))* &
        4.0d0
! array qran
      rmem2=rmem2+3.0D0*dble(nmc)*8.0D0
    endif
! cell box truncation (parallel version only)
    if (sig%icutv==TRUNC_BOX) then
      dkmax(1) = gvec%FFTgrid(1) * n_in_box
      dkmax(2) = gvec%FFTgrid(2) * n_in_box
      dkmax(3) = gvec%FFTgrid(3) * n_in_box
      call setup_FFT_sizes(dkmax,dNfft,dscale)
      if (mod(dNfft(3),peinf%npes) == 0) then
        Nplane = dNfft(3)/peinf%npes
      else
        Nplane = dNfft(3)/peinf%npes+1
      endif
      if (mod(dNfft(1)*dNfft(2),peinf%npes) == 0) then
        Nrod = (dNfft(1)*dNfft(2))/peinf%npes
      else
        Nrod = (dNfft(1)*dNfft(2))/peinf%npes+1
      endif
! array fftbox_2D
      rmem2=rmem2+dble(dNfft(1))*dble(dNfft(2))*dble(Nplane)* &
        16.0d0
! array fftbox_1D
      rmem2=rmem2+dble(dNfft(3))*dble(Nrod)*16.0d0
! array dummy
! rmem2=rmem2+dble(dNfft(1))*dble(dNfft(2))*16.0d0
! arrays dummy1 and dummy2
      rmem2=rmem2+dble(Nrod)*dble(peinf%npes+1)*16.0d0
! array inv_indx
      rmem2=rmem2+dble(Nfft(1))*dble(Nfft(2))*dble(Nfft(3))* &
        4.0d0
    endif
    if (rmem2 .gt. rmem) rmem = rmem2
988 format(1x,'Memory required for vcoul: ',f0.1,' MB per PE')
    write(6,988) rmem/1024.0d0**2
    write(6,*)
  endif
  if(peinf%inode == 0) then
!----------------------------------------------
! Compute cell volume from wave number metric
    call get_volume(vcell,crys%bdot)
    if (abs(crys%celvol-vcell).gt.TOL_Small) then
      call die('volume mismatch')
    endif
! (gsm) check consistency of spin indices
    do k=1,sig%nspin
      if (sig%spin_index(k).lt.1.or.sig%spin_index(k).gt.kp%nspin) &
        call die('inconsistent spin indices')
    enddo
    if(sig%ntband.gt.kp%mnband) then
      write(tmpstr1,660) sig%ntband
      write(tmpstr2,660) kp%mnband
      write(0,666) TRUNC(tmpstr1), TRUNC(tmpstr2)
      if(.not.ep_read) then
        call die('More bands specified in sigma.inp than available in WFN_inner.')
      else
        call die('More bands specified in sigma.inp than available in dWFN.')
      endif
660 format(i16)
666 format(1x,'The total number of bands (',a,') specified in sigma.inp',/, &
        3x,'is larger than the number of bands (',a,') available in WFN_inner.',/)
    endif
    do_ch_sum = .not.((sig%freq_dep == -1) .or. (sig%freq_dep == 0 .and. sig%exact_ch == 1))
    if(sig%ntband.eq.kp%mnband) then
      if(.not.ep_read) then
        call die("You must provide one more band in WFN_inner than used in sigma.inp number_bands in order to assess degeneracy.")
      else
        call die("You must provide one more band in dWFN than used in sigma.inp number_bands in order to assess degeneracy.")
      endif
    endif
!--------------------------------------------------------------------
! SIB: Find the k-points in sig%kpt in the list kp%rk (die if not found).
! sig%indkn holds the index of the k-points in sig%kpt in kp%rk, i.e.
! kp%rk(sig%indkn(ikn))=sig%kpt(ikn)
    allocate(sig%indkn (sig%nkn))
    do ikn=1,sig%nkn
      sig%indkn(ikn)=0
      qk(:)=sig%kpt(:,ikn)
      do irk=1,kp%nrk
        if(all(abs(kp%rk(1:3,irk)-qk(1:3)) .lt. TOL_Small)) sig%indkn(ikn)=irk
      enddo
      ! ZL: print information
      if(do_phonq .and. peinf%inode.eq.0) write(6,*) "sig%indkn(",ikn,") =", sig%indkn(ikn)
      if(sig%indkn(ikn) .eq. 0) then
        if(.not.ep_read) then
          write(0,'(a,3f10.5,a)') 'Could not find the k-point ', (qk(i),i=1,3), ' among those read from WFN_inner :'
        else
          write(0,'(a,3f10.5,a)') 'Could not find the k-point ', (qk(i),i=1,3), ' among those read from dWFN :'
        endif
        write(0,'(3f10.5)') ((kp%rk(i,irk),i=1,3),irk=1,kp%nrk)
        call die('k-point in sigma.inp k_points block not available.')
      endif
    enddo
    ! ZL: EP, when reading regular WFN_inner, NOT dWFN, we map k+phonq points
    if(do_phonq .and. (.not. ep_read) .and. (sig%nphonq .eq. 1)) then ! ZL: only nphonq=1 supported for now
      allocate(sig%indk_phonq (sig%nkn * sig%nphonq))
      allocate(sig%indk_phonq_g0 (3, sig%nkn * sig%nphonq))
      do ikn=1,sig%nkn*sig%nphonq
        sig%indk_phonq(ikn)=0
        qk(:)=sig%k_phonq(:,ikn)
        irk_loop_g0: do irk=1,kp%nrk
          do i_looper=1,3
            ! ZL: this idea is taken from genwf_mpi.f90
            ! the most interesting part is that
            ! it uses the fact that g0(3) is INTEGER
            ! therefore it will automatically drop the fractional part
            ! and keep only the integer part
            ! This is also guaranteed by adding TOL_small so that 0.999999
            ! is indeed 1
            ! Then we just need to compare if del is close to integer
            ! NOTE: gmap uses a different idea, it maps between old and new wfn
            ! but here we map to master G-list
            del = qk(i_looper) - kp%rk(i_looper,irk)
            if (del .ge. 0.0d0) g0(i_looper) = del + TOL_Small
            if (del .lt. 0.0d0) g0(i_looper) = del - TOL_Small
            if (abs(del-g0(i_looper)) .gt. TOL_Small) cycle irk_loop_g0
          enddo
          sig%indk_phonq(ikn)=irk
          sig%indk_phonq_g0(:,ikn)=g0(:)
          exit irk_loop_g0
        enddo irk_loop_g0
        ! ZL: print information
        if(peinf%inode .eq. 0) then
          write(6,*) "sig%indk_phonq(",ikn,") =", sig%indk_phonq(ikn)
          write(6,*) "sig%indk_phonq_g0(",ikn,") =", sig%indk_phonq_g0(:,ikn)
        endif
        if(sig%indk_phonq(ikn) .eq. 0) then
          write(0,'(a,3f10.5,a)') 'Could not find the k_phonq-point ', (qk(i),i=1,3), ' among those in WFN_inner :'
          write(0,'(3f10.5)') ((kp%rk(i,irk),i=1,3),irk=1,kp%nrk)
          call die('k_phonq-point from sigma.inp k_points+phonq not available.')
        endif
      enddo
    endif ! EP k+phonq mapping
    if(do_ch_sum .and. (.not.ep_read)) then
      if(any(abs(kp%el(sig%ntband+sig%ncore_excl, 1:kp%nrk, 1:kp%nspin) - &
                 kp%el(sig%ntband+sig%ncore_excl+1, 1:kp%nrk, 1:kp%nspin)) .lt. sig%tol/RYD)) then
        if(sig%degeneracy_check_override) then
          write(0,'(a)') &
            "WARNING: Selected number of bands for CH sum (number_bands) breaks degenerate subspace. " // &
            "Run degeneracy_check.x for allowable numbers."
          write(0,*)
        else
          write(0,'(a)') &
            "Run degeneracy_check.x for allowable numbers, or use keyword " // &
            "degeneracy_check_override to run anyway (at your peril!)."
          call die("Selected number of bands for CH sum (number_bands) breaks degenerate subspace.")
        endif
      endif
    endif
  endif
!----------------------------------------------------------------
! Check consistency
  if (kp%mnband < maxval(sig%diag) .and. peinf%inode == 0) then
    write(0,*) 'The highest requested band is ', maxval(sig%diag),' but WFN_inner contains only ', kp%mnband,' bands.'
    call die('too few bands')
  endif
  ! qgridsym has no meaning when there is only one q-point (Gamma)
  ! setting this to false means warnings about its applicability will not be triggered
  if(all(abs(kp%rk(1:3,1:kp%nrk)) < TOL_Zero)) sig%qgridsym = .false.
  if(peinf%inode == 0) then
    ! ZL: assess_degeneracies will allocate kp%degeneracy
    ! kp%el is still in Ry here, but sig%tol is in eV.
    call assess_degeneracies(kp, kp%el(sig%ntband+sig%ncore_excl+1, :, :), sig%ntband, sig%efermi, sig%tol/RYD, sig = sig, &
                             ncore_excl=sig%ncore_excl)
    if(.not.ep_read) then
      call calc_qtot(kp, crys%celvol, sig%efermi, qtot, omega_plasma, write7=.false.)
    endif
  endif
  call timing%start(timing%fullbz)
  gr%nr = kp%nrk
  allocate(gr%r (3, gr%nr))
  gr%r = kp%rk
  call fullbz(crys,syms,gr,syms%ntran,skip_checkbz,wigner_seitz=.false.,paranoid=.true.)
  if(.not.ep_read) then
    tmpfn='WFN_inner'
  else
    tmpfn='dWFN'
  endif
  if (.not. skip_checkbz) then
    call checkbz(gr%nf,gr%f,kp%kgrid,kp%shift,crys%bdot, &
      tmpfn,'k',.false.,sig%freplacebz,sig%fwritebz)
  endif
  call dealloc_grid(gr)
  call timing%stop(timing%fullbz)
  ! For discussion of how q-symmetry may (and may not) be used with degenerate states,
  ! see Hybertsen and Louie, Phys. Rev. B 35, 5585 (1987) Appendix A
  if(sig%qgridsym .and. sig%noffdiag > 0) then
    if(peinf%inode == 0) then
      write(0,'(a)') "WARNING: Cannot calculate offdiagonal elements unless no_symmetries_q_grid is set."
      write(0,'(a)') "This flag is being reset to enable the calculation."
    endif
    sig%qgridsym = .false.
  endif
  ! ZL: for EP, turn off q-symmetry. Interestingly, this coincides with off-diag
  if(sig%qgridsym .and. ep_read) then
    if(peinf%inode == 0) then
      write(0,'(a)') "WARNING: Cannot calculate electron phonon unless no_symmetries_q_grid is set."
      write(0,'(a)') "This flag is being reset to enable the calculation."
    endif
    sig%qgridsym = .false.
  endif
  if(peinf%inode == 0) then
    if(sig%qgridsym) then
      write(6,'(1x,a)') 'Q-grid symmetries are being used.'
    else
      write(6,'(1x,a)') 'Q-grid symmetries are not being used.'
    endif
  endif
  kp%el(:,:,:)=kp%el(:,:,:)*ryd
 !----------------------------------------------------------------
 ! Distribute data
!----------------------------------------------------------------
! CHP: Set the init_freq_eval relative to the Fermi energy in case
! of full frequency. This part should NOT be put before the
! Fermi energy is set.
  if ((sig%freq_dep==2 .or. (sig%freq_dep==1.and.sig%fdf==-3)) .and. sig%freq_grid_shift==1) then
    ! FHJ: only do this for FE-based freq. grids
    sig%freqevalmin = sig%freqevalmin + sig%efermi
  endif
!--------------------------------------------------------------------
! Read in the exchange-correlation potential and store in array vxc
  if(.not.sig%use_vxcdat .and. .not.sig%sigma_correction .and. .not. sig%is_EXX .and. .not.sig%use_kihdat) then
    ! ZL: include condition that when not running with kih.dat
    call read_vxc(sig, gvec, kp, syms, crys, isrti, isrt, 1, tolerant_value=sig%tolerant_value)
  endif ! not using vxc.dat
  if(.not. sig%use_vxc2dat) then
    ! This is for hybrid functional like calculation (one shot)
    if (sig%freq_dep.eq.-1 .and. ((1.0d0 - sig%coulomb_mod%long_range_frac_fock > TOL_SMALL) .or. &
                  (1.0d0 - sig%coulomb_mod%short_range_frac_fock > TOL_SMALL))) &
      call read_vxc(sig, gvec, kp, syms, crys, isrti, isrt, 2, tolerant_value=sig%tolerant_value)
  endif ! not using vxc2.dat
  ! ZL: read dVXC when reading dWFN
  if (sig%elph .and. ep_read) then
    call read_vxc(sig, gvec, kp, syms, crys, isrti, isrt, 3, tolerant_value=sig%tolerant_value)
  endif
!--------------------------------------------------------------------
! Read in the charge density and store in array rho (formerly known as CD95)
! CD95 Ref: http://www.nature.com/nature/journal/v471/n7337/full/nature09897.html
  if(sig%freq_dep.eq.1) then
    call timing%start(timing%input_read)
    call read_rho(wpg, gvec, kp, syms, crys, isrti, isrt, 'WFN_inner', tolerant_value=sig%tolerant_value) ! ZL: only read from WFN_inner, NOT dWFN, for RHO
    rhog0 = sum(wpg%nelec(1:kp%nspin))
    if(peinf%inode == 0 .and. abs(rhog0 - qtot) > TOL_Small) then
      write(0,'(a,g12.6,a,g12.6)') 'WARNING: RHO has total charge ', rhog0, ' but WFN_inner occupations give charge ', qtot
    endif
    call timing%stop(timing%input_read)
  endif ! sig%freq_dep
!------------------------------------------------------------------------
! Divide bands over processors
!
! sig%ntband number of total (valence and conduction) bands
! peinf%npools number of parallel sigma calculations
! peinf%ndiag_max maximum number of diag sigma calculations
! peinf%noffdiag_max maximum number of offdiag sigma calculations
! peinf%ntband_max maximum number of total bands per node
!
! ntband_node number of total bands per current node
! nvband_node number of valence bands per current node
! peinf%indext(itb) indices of total bands belonging to current node
!
! peinf%index_diag band index for diag sigma calculation
! peinf%flag_diag flag for storing diag sigma calculation
! peinf%index_offdiag band index for offdiag sigma calculation
! peinf%flag_offdiag flag for storing offdiag sigma calculation
!
! flags are needed in case of uneven distribution, nodes still do the
! calculation because they share epsilon and wavefunctions and need to
! participate in global communications, but the result is not stored
  if (peinf%npools .le. 0 .or. peinf%npools .gt. peinf%npes) then
    call createpools(sig%ndiag,sig%ntband,peinf%npes,npools,ndiag_max,ntband_max)
    peinf%npools = npools
    peinf%ndiag_max = ndiag_max
    peinf%ntband_max = ntband_max
  else
    if (mod(sig%ndiag,peinf%npools).eq.0) then
      peinf%ndiag_max=sig%ndiag/peinf%npools
    else
      peinf%ndiag_max=sig%ndiag/peinf%npools+1
    endif
    if (mod(sig%ntband,peinf%npes/peinf%npools).eq.0) then
      peinf%ntband_max=sig%ntband/(peinf%npes/peinf%npools)
    else
      peinf%ntband_max=sig%ntband/(peinf%npes/peinf%npools)+1
    endif
  endif
! JRD Create Separate Comms for Each Pool
  peinf%npes_pool = peinf%npes/peinf%npools
  peinf%my_pool=peinf%inode/peinf%npes_pool
  peinf%pool_rank = mod(peinf%inode,peinf%npes_pool)
  if (sig%noffdiag.gt.0) then
    if (mod(sig%noffdiag,peinf%npools).eq.0) then
      peinf%noffdiag_max=sig%noffdiag/peinf%npools
    else
      peinf%noffdiag_max=sig%noffdiag/peinf%npools+1
    endif
  endif
  allocate(peinf%index_diag (peinf%ndiag_max))
  allocate(peinf%flag_diag (peinf%ndiag_max))
  peinf%index_diag=1
  peinf%flag_diag=.false.
  do ii=1,peinf%ndiag_max*peinf%npools
    jj=(ii-1)/peinf%npools+1 ! which number in the pool it is
    kk=mod(ii-1,peinf%npools) ! which pool it is in
    if (peinf%inode/(peinf%npes/peinf%npools).eq.kk) then
      if (ii.le.sig%ndiag) then
        peinf%index_diag(jj)=ii
        peinf%flag_diag(jj)=.true.
      else
        peinf%index_diag(jj)=1
        peinf%flag_diag(jj)=.false.
      endif
    endif
  enddo
  if (sig%noffdiag.gt.0) then
    allocate(peinf%index_offdiag (peinf%noffdiag_max))
    allocate(peinf%flag_offdiag (peinf%noffdiag_max))
    peinf%index_offdiag=1
    peinf%flag_offdiag=.false.
    do ii=1,peinf%noffdiag_max*peinf%npools
      jj=(ii-1)/peinf%npools+1
      kk=mod(ii-1,peinf%npools)
      if (peinf%inode/(peinf%npes/peinf%npools).eq.kk) then
        if (ii.le.sig%noffdiag) then
          peinf%index_offdiag(jj)=ii
          peinf%flag_offdiag(jj)=.true.
        else
          peinf%index_offdiag(jj)=1
          peinf%flag_offdiag(jj)=.false.
        endif
      endif
    enddo
  endif
  allocate(peinf%indext (peinf%ntband_max))
!JRD Now is over procs per pool
  allocate(peinf%indext_dist (peinf%ntband_max,peinf%npes_pool))
  allocate(peinf%ntband_dist (peinf%npes_pool))
  peinf%ntband_node=0
  peinf%ntband_dist=0
  peinf%nvband_node=0
  peinf%indext=0
  peinf%indext_dist=0
  do ii=1,sig%ntband
      ! FHJ: The bands are assigned to the processors in a sequential fashion.
      ! This is more efficient for the parallel IO. Eg: for 2 bands/processor,
      ! v1->PE0, v2->PE0, v3->PE1, etc.
      ipe = (ii-1)/peinf%ntband_max
      if (peinf%pool_rank.eq.ipe) then
        peinf%ntband_node=peinf%ntband_node+1
        if(ii.le.(sig%nvband+sig%ncrit)) then
          peinf%nvband_node=peinf%nvband_node+1
        endif
        peinf%indext(peinf%ntband_node)=ii
      endif
      peinf%ntband_dist(ipe+1)=peinf%ntband_dist(ipe+1)+1
      peinf%indext_dist(peinf%ntband_dist(ipe+1),ipe+1)=ii
  enddo
!------------------------------------------------------------------------
! Report distribution of bands over processors
  if (peinf%inode==0) then
    write(6,*)
    write(6,'(1x,a)') 'Parallelization report:'
701 format(1x,'- Using ',i0,' processor(s), ' ,i0,' pool(s), ',i0,' processor(s) per pool.')
    write(6,701) peinf%npes, peinf%npools, &
      peinf%npes/peinf%npools
702 format(1x,' - Note: distribution is not ideal; ',i0,' processor(s) is/are idle.')
    if (mod(peinf%npes,peinf%npools).ne.0) &
      write(6,702) peinf%npes-(peinf%npes/peinf%npools)* &
      peinf%npools
    if (mod(sig%ndiag,peinf%npools).eq.0) then
703 format(1x,'- Each pool is computing ',i0,' diagonal sigma matrix element(s).')
      write(6,703) peinf%ndiag_max
    else
704 format(1x,'- Each pool is computing ',i0,' to ',i0,' diagonal sigma matrix element(s).')
      write(6,704) peinf%ndiag_max-1, peinf%ndiag_max
705 format(1x,'- Note: distribution is not ideal because the number of diagonal sigma',/,&
      3x,'matrix elements (',i0,') is not divisible by the number of pools (',i0,').')
      write(6,705) sig%ndiag, peinf%npools
    endif
    if (mod(sig%noffdiag,peinf%npools).eq.0) then
706 format(1x,'- Each pool is computing ',i0,' off-diagonal sigma matrix element(s).')
      write(6,706) peinf%noffdiag_max
    else
707 format(1x,'- Each pool is computing ',i0,' to ',i0,' off-diagonal sigma matrix element(s).')
      write(6,707) peinf%noffdiag_max-1, peinf%noffdiag_max
708 format(1x,'- Note: distribution is not ideal because the number of off-diagonal sigma',/,&
      3x,'matrix elements (',i0,') is not divisible by the number of pools (',i0,').')
      write(6,708) sig%noffdiag, peinf%npools
    endif
    if (mod(sig%ntband,peinf%npes/peinf%npools).eq.0) then
709 format(1x,'- Each processor is holding ',i0, ' band(s).')
      write(6,709) peinf%ntband_max
    else
710 format(1x,'- Each pool is holding ',i0,' to ',i0,' band(s).')
      write(6,710) minval(peinf%ntband_dist), maxval(peinf%ntband_dist)
711 format(1x,'- Note: distribution is not ideal because the total number of bands',/,&
      3x,'(',i0,') is not divisible by the number of processors per pool (',i0,').')
      write(6,711) sig%ntband, peinf%npes/peinf%npools
    endif
    write(6,*)
  endif
!-----------------------------------------------------------
!
! LOOP OVER K-POINT GRID AND READ IN WAVEFUNCTIONS
!
!-----------------------------------------------------------
  if(sig%wfn_hdf5) then
  else
    ! ZL: we call read_wavefunctions() only once in input()
    if(.not.ep_read) then
      if(.not.do_phonq) then
        call read_wavefunctions(kp, gvec, sig, wfnk, iunit_c, iunit_k, fnc, fnk, wfnkqmpi, wfnkmpi)
      else
        call read_wavefunctions(kp, gvec, sig, wfnk, iunit_c, iunit_k, fnc, fnk, wfnkqmpi, wfnkmpi,&
                                ep_read_in=.false., do_phonq_in=do_phonq, wfnk_phonq=wfnk_phonq, wfnk_phonq_mpi=wfnk_phonq_mpi)
      endif
    else
      call read_wavefunctions(kp, gvec, sig, wfnk, iunit_c, iunit_k, fnc, fnk, wfnkqmpi, wfnkmpi, ep_read_in = ep_read)
    endif
    if(peinf%inode.eq.0) then
      call close_file(25)
    endif
  endif
  if(peinf%inode.eq.0) then
! Write out information about self-energy calculation
    call scissors_write(6, sig%scis)
    call scissors_write(6, sig%scis_outer, "outer")
    write(6,'(a)')
    if(sig%freq_dep.eq.1) then
      do k=1,sig%nspin
        write(6,160) sig%spin_index(k),wpg%nelec &
          (sig%spin_index(k)),sqrt(wpg%wpsq(sig%spin_index(k)))
160 format(/,1x,'Data for sum rule:',/,&
          1x,'- rho(0,',i0,') = ',f0.6,' electrons',/,&
          1x,'- wp = ',f0.6,' eV',/)
      enddo
    endif ! sig%freq_dep
    write(6,'(1x,a,i0)') 'Number of bands to compute diagonal self-energy matrix elements: ', sig%ndiag
    write(6,'(1x,a)') 'Bands:'
    write(6,'(1x,"- ",i0)') sig%diag(:)
    write(6,'(1x,a,i0)') 'Number of off-diagonal bands to compute self-energy matrix elements: ', sig%noffdiag
    if (sig%noffdiag>0) then
      write(6,'(1x,a)') 'Off-diagonal bands:'
      do k=1,sig%noffdiag
174 format(1x,'- offmap(:, ',i6,') =',3i6)
        write(6,174) k, (sig%offmap(k,ii), ii = 1, 3)
      enddo
    endif
    write(6,*)
  endif ! node 0
  if(allocated(isrt))then;deallocate(isrt);endif
  if(allocated(isrti))then;deallocate(isrti);endif
  if(allocated(kp%w))then;deallocate(kp%w);endif
  if(allocated(kp%el))then;deallocate(kp%el);endif
  if(allocated(kp%elda))then;deallocate(kp%elda);endif
 
  return
end subroutine input
subroutine read_wavefunctions(kp,gvec,sig,wfnk,iunit_c,iunit_k,fnc,fnk,wfnkqmpi,wfnkmpi,read_hdf5, &
                              ep_read_in,do_phonq_in,wfnk_phonq,wfnk_phonq_mpi)
  type (kpoints), intent(in) :: kp
  type (gspace), intent(in) :: gvec
  type (siginfo), intent(in) :: sig
  type (wfnkstates), intent(out) :: wfnk
  integer, intent(in) :: iunit_c, iunit_k
  character*20, intent(in) :: fnc, fnk
  type (wfnkqmpiinfo), intent(out) :: wfnkqmpi
  type (wfnkmpiinfo), intent(out) :: wfnkmpi
  ! hdf5 wfn
  logical, optional, intent(in) :: read_hdf5
  ! ZL: for EP
  type (wfnkstates), optional, intent(out) :: wfnk_phonq ! similar to wfnk
  type (wfnkmpiinfo), optional, intent(out) :: wfnk_phonq_mpi ! similar to wfnkmpi
  character :: tmpstr*100,tmpstr1*16,tmpstr2*16
  integer :: i,ii,j,k
  integer :: ikn, irk
  integer :: istore,nstore,inum,g1,g2,iknstore, istore_ep, iknstore_ep ! ZL adds istore_ep, iknstore_ep
  integer :: ndvmax, ndvmax_l
  integer, allocatable :: isort(:)
  real(DP) :: qk(3)
  real(DP), allocatable :: zc(:,:)
  logical :: dont_read
  type(gspace) :: gvec_kpt
  type(progress_info) :: prog_info !< a user-friendly progress report
  ! ZL: electron phonon (EP)
  logical, intent(in), optional :: ep_read_in, do_phonq_in
  logical :: ep_read, do_phonq
  integer :: induse ! index in use, for indkn or indk_phonq
  integer :: tmpngk
  ! WFN hdf5 variables
  logical :: my_read_hdf5
  real(DP), allocatable :: zc_all(:,:,:)
  integer :: error
  integer :: ngktot
  integer :: istart, ib_first, ib_last
  integer, allocatable :: components(:,:)
  logical :: i_got_this
 
  my_read_hdf5 = .false.
  if ( present(read_hdf5) ) my_read_hdf5 = read_hdf5
  ! ZL: set up EP
  ep_read = .false.
  if(present(ep_read_in)) then
    ep_read = ep_read_in
  endif
  do_phonq = .false.
  if(present(do_phonq_in)) then
    do_phonq = do_phonq_in
  endif
  ! ZL: although it should never happen, in case someone mis-uses, double check here
  if(do_phonq .and. ep_read) then
    call die("Doing k+phonq points (designed for outer) is NOT allowed for dWFN (only used as inner).", &
              only_root_writes = .true.)
  endif
  if (do_phonq) then
    if ((.not.present(wfnk_phonq)) .or. (.not.present(wfnk_phonq_mpi))) then
      call die("Doing k_phonq wfn needs wfklkrq and wfnk_phonq_mpi as arguments!", only_root_writes = .true.)
    endif
  else
    if (present(wfnk_phonq) .or. present(wfnk_phonq_mpi)) then
      call die("You are not doing k_phonq wfn, then do not pass wfnk_phonq or wfnk_phonq_mpi.", only_root_writes=.true.)
    endif
  endif
  ! ZL: done EP setup
  allocate(isort (gvec%ng))
  allocate(wfnkqmpi%nkptotal (kp%nrk))
  allocate(wfnkqmpi%isort (kp%ngkmax,kp%nrk))
  allocate(wfnkqmpi%band_index (peinf%ntband_max,kp%nrk))
  allocate(wfnkqmpi%qk (3,kp%nrk))
  allocate(wfnkqmpi%el (sig%ntband,sig%nspin,kp%nrk))
  ! ZL: ntband_max is total number of bands distributed to each processor
  ! wfnkqmpi%cg is distributed over bands
  allocate(wfnkqmpi%cg (kp%ngkmax,peinf%ntband_max,sig%nspin*kp%nspinor,kp%nrk))
  if (sig%nkn.gt.1) then ! ZL: same logic for EP
    ! it is using ngkmax, so same for EP, wfnkmpi%cg is distributed over ndv
    ndvmax=peinf%ndiag_max*kp%ngkmax
    if (mod(ndvmax,peinf%npes/peinf%npools).eq.0) then
      ndvmax_l=ndvmax/(peinf%npes/peinf%npools)
    else
      ndvmax_l=ndvmax/(peinf%npes/peinf%npools)+1
    endif
    allocate(wfnkmpi%nkptotal (sig%nkn))
    allocate(wfnkmpi%isort (kp%ngkmax,sig%nkn))
    allocate(wfnkmpi%qk (3,sig%nkn))
    allocate(wfnkmpi%el (sig%ntband,sig%nspin,sig%nkn))
    allocate(wfnkmpi%elda (sig%ntband,sig%nspin,sig%nkn))
    allocate(wfnkmpi%cg (ndvmax_l,sig%nspin*kp%nspinor,sig%nkn))
    ! ZL: allocate wfnk_phonq_mpi
    ! if sig%nkn=1, then no mpi, only use wfnk_phonq which stores k+phonq
    if (do_phonq) then
      ! ZL: the size should in principle be sig%nkn*sig%nphonq,
      ! since we consider only nphonq=1, we omit it here
      allocate(wfnk_phonq_mpi%nkptotal (sig%nkn))
      allocate(wfnk_phonq_mpi%isort (kp%ngkmax,sig%nkn))
      allocate(wfnk_phonq_mpi%qk (3,sig%nkn))
      allocate(wfnk_phonq_mpi%el (sig%ntband,sig%nspin,sig%nkn))
      allocate(wfnk_phonq_mpi%elda (sig%ntband,sig%nspin,sig%nkn))
      allocate(wfnk_phonq_mpi%cg (ndvmax_l,sig%nspin*kp%nspinor,sig%nkn))
    endif
  endif
!-----------------------------------
! Read in wavefunction information
  allocate(wfnk%isrtk (gvec%ng))
  allocate(wfnk%ek (sig%ntband,sig%nspin))
  allocate(wfnk%elda (sig%ntband,sig%nspin))
  ! ZL: allocate wfnk_phonq
  if (do_phonq) then
    allocate(wfnk_phonq%isrtk (gvec%ng))
    allocate(wfnk_phonq%ek (sig%ntband,sig%nspin))
    allocate(wfnk_phonq%elda (sig%ntband,sig%nspin))
  endif
  if(.not.ep_read) then
    call progress_init(prog_info, 'reading wavefunctions (WFN_inner)', 'state', kp%nrk*sig%ntband)
  else
    call progress_init(prog_info, 'reading wavefunctions (dWFN)', 'state', kp%nrk*sig%ntband)
  endif
  ! Actual HDF5 WFN reading
  if ( my_read_hdf5 ) then
    if ( peinf%inode == 0 ) then
      write(6,*) 'Reading band block from WFN_inner.h5'
    end if
    ngktot = SUM(kp%ngk)
    allocate(gvec_kpt%components (3, ngktot))
    allocate(zc_all (ngktot, kp%nspin*kp%nspinor , peinf%ntband_node))
  end if
  ! ZL: outmost loop, irk over kp%nrk
  ! keep track of the current band block for each k-point
  istart = 1
  do irk=1,kp%nrk
    if(.not.ep_read) then
      write(tmpstr,*) 'Reading WFN_inner -> cond/val wfns irk=',irk
    else
      write(tmpstr,*) 'Reading dWFN -> cond/val wfns irk=',irk
    endif
    call logit(tmpstr)
    qk(:)=kp%rk(:,irk)
!----------------------------
! Read in and sort gvectors
    if ( my_read_hdf5 ) then
      ! we already have all components, just make a copy
      allocate(components (3,kp%ngk(irk)))
      components(:,:) = gvec_kpt%components(:,istart:istart+kp%ngk(irk)-1)
      do i = 1, kp%ngk(irk)
        call findvector(isort(i), components(:, i), gvec)
        if (isort(i) .eq. 0) then
          write(0,*) 'could not find gvec', kp%ngk(irk), i, components(1:3, i)
          call die('findvector')
        endif
      enddo
      if(allocated(components))then;deallocate(components);endif
    else
      allocate(gvec_kpt%components (3, kp%ngk(irk)))
      call timing%start(timing%input_read)
      call read_binary_gvectors(25, kp%ngk(irk), kp%ngk(irk), gvec_kpt%components)
      call timing%stop(timing%input_read)
      do i = 1, kp%ngk(irk)
        call findvector(isort(i), gvec_kpt%components(:, i), gvec)
        if (isort(i) .eq. 0) then
          write(0,*) 'could not find gvec', kp%ngk(irk), i, gvec_kpt%components(1:3, i)
          call die('findvector')
        endif
      enddo
      if(associated(gvec_kpt%components))then;deallocate(gvec_kpt%components);nullify(gvec_kpt%components);endif
    end if ! my_read_hdf5
!--------------------------------------------------------
! Determine if Sigma must be computed for this k-point.
! If so, store the bands and wavefunctions on file iunit_k.
! If there is only one k-point, store directly in wfnk.
    istore=0
    if(do_phonq) istore_ep=0 ! ZL: initialize
    do ikn=1,sig%nkn ! ZL: here we take advantage that we implement for one phonq point for now
                     ! otherwize, need another outer loop: do ikn=1,sig%nkn*sig%nphonq
      if(sig%indkn(ikn).eq.irk) then ! ZL: KEY if-statement: here it judges whether store wfnkmpi
        istore=1 ! ZL: for a given irk point, only one ind(ikn) index can match
        iknstore=ikn ! ZL: and this ind(ikn) is recorded as iknstore
        wfnk%nkpt=kp%ngk(irk)
        wfnk%ndv=peinf%ndiag_max*kp%ngk(irk)
        wfnk%k(:)=qk(:)
        wfnk%isrtk(:)=isort(:)
        do k=1,sig%nspin
          wfnk%ek(1:sig%ntband,k)= &
            kp%el(1+ sig%ncore_excl:sig%ntband+ sig%ncore_excl,irk,sig%spin_index(k))
          wfnk%elda(1:sig%ntband,k)= &
            kp%elda(1:sig%ntband,irk,sig%spin_index(k))
        enddo
        allocate(wfnk%zk (wfnk%ndv,sig%nspin*kp%nspinor))
        wfnk%zk=0.0d0
      endif
      ! ZL: for do_phonq
      if(do_phonq) then
        if (sig%indk_phonq(ikn).eq.irk) then
          istore_ep=1
          iknstore_ep=ikn
          wfnk_phonq%nkpt=kp%ngk(irk)
          wfnk_phonq%ndv=peinf%ndiag_max*kp%ngk(irk)
          wfnk_phonq%k(:)=qk(:) + sig%indk_phonq_g0(:,ikn)
          wfnk_phonq%isrtk(:)=isort(:)
          call shift_phonq_g0(sig%indk_phonq_g0(:,ikn), wfnk_phonq, gvec)
          do k=1,sig%nspin
            wfnk_phonq%ek(1:sig%ntband,k)= &
              kp%el(1+ sig%ncore_excl:sig%ntband+ sig%ncore_excl,irk,sig%spin_index(k))
            wfnk_phonq%elda(1:sig%ntband,k)= &
              kp%elda(1:sig%ntband,irk,sig%spin_index(k))
          enddo
          allocate(wfnk_phonq%zk (wfnk_phonq%ndv,sig%nspin*kp%nspinor))
          wfnk_phonq%zk=0.0d0
        endif
      endif
    enddo ! ikn
    ! ZL: nkptotal is ngk
    wfnkqmpi%nkptotal(irk) = kp%ngk(irk)
    wfnkqmpi%isort(1:kp%ngk(irk),irk) = isort(1:kp%ngk(irk))
    do k = 1, sig%nspin
      wfnkqmpi%el(1:sig%ntband,k,irk) = &
        kp%el(1+ sig%ncore_excl:sig%ntband+ sig%ncore_excl,irk,sig%spin_index(k))
    enddo
    wfnkqmpi%qk(1:3,irk) = qk(1:3)
!-------------------------------------------------------------------------
! SIB: Read wave functions from file WFN_inner (ZL: or dWFN) (unit=25) and have
! the appropriate processor write it to iunit_c after checking norm.
! The wavefunctions for bands where Sigma matrix elements are
! requested are stored in wfnk%zk for later writing to unit iunit_k.
! If band index is greater than sig%ntband, we will actually
! not do anything with this band (see code below).
! We still have to read it though, in order
! to advance the file to get to the data for the next k-point.
    allocate(zc (kp%ngk(irk), kp%nspin*kp%nspinor))
    inum=0
    do i=1,kp%mnband
      call timing%start(timing%input_read)
! DVF: don`t read deep core states.
      dont_read = (i > sig%ntband+sig%ncore_excl .or. i <= sig%ncore_excl)
      i_got_this = .false.
      if ( my_read_hdf5 ) then
        ! hdf5 version all bands have already be read and stored in zc_all
        if ( .not.dont_read ) then
          if(i-sig%ncore_excl >= ib_first .and. i-sig%ncore_excl <= ib_last) then
            i_got_this = .true.
            do k=1,sig%nspin*kp%nspinor
              zc(:,k) = zc_all( istart:(istart + kp%ngk(irk)-1), k , i-sig%ncore_excl-ib_first+1)
            end do
          end if
        end if
      else
        ! binary version
        call read_binary_data(25, kp%ngk(irk), kp%ngk(irk), kp%nspin*kp%nspinor, zc, dont_read = dont_read)
      end if
      call timing%stop(timing%input_read)
      if (.not.dont_read) then
        call progress_step(prog_info, sig%ntband*(irk-1) + i)
        if(peinf%inode.eq.0) then
          do k=1,sig%nspin
            if(.not.ep_read) then
              call checknorm('WFN_inner',i,irk,kp%ngk(irk),k,kp%nspinor,zc(:,:))
            endif
          enddo
        endif
        nstore=0 ! ZL: for bands
        do j=1,peinf%ndiag_max
          if (.not.peinf%flag_diag(j)) cycle
          if (i==sig%diag(peinf%index_diag(j))) nstore=j
        enddo
        if ( my_read_hdf5 ) then
          ! all reduce such all pe in the pools have this current (outer) WFN
          if ( (istore.eq.1).and.(nstore.ne.0)) then
            if ( .not. i_got_this ) zc = 0.0d0
          end if
        end if
        ! ZL: if this is the kpoint in outer, save in wfnk
        if((istore.eq.1).and.(nstore.ne.0)) then
          do k=1,sig%nspin*kp%nspinor
            do j=1,kp%ngk(irk)
              wfnk%zk((nstore-1)*kp%ngk(irk)+j,k) = zc(j,sig%spin_index(k))
            enddo
          enddo
        endif
        ! ZL: for do_phonq
        if(do_phonq.and.(istore_ep.eq.1).and.(nstore.ne.0)) then
          do k=1,sig%nspin*kp%nspinor
            do j=1,kp%ngk(irk)
              wfnk_phonq%zk((nstore-1)*kp%ngk(irk)+j,k) = zc(j,sig%spin_index(k))
            enddo
          enddo
        endif
! DVF : i is indexed including the deep core states (since it runs over all the bands
! in the wavefunction), while the indexing arrays
! are not. We have to subtract sig%ncore_excl from i to get the right states.
        j=0 ! ZL: index to decide if save this band
        do ii=1,peinf%ntband_node
          if (peinf%indext(ii)==i-sig%ncore_excl) j=1
        enddo
        call timing%start(timing%input_write)
        if (j.eq.1) then
          inum=inum+1
          wfnkqmpi%band_index(inum,irk)=i-sig%ncore_excl
          do k=1,sig%nspin*kp%nspinor
            wfnkqmpi%cg(1:kp%ngk(irk),inum,k,irk)= &
              zc(1:kp%ngk(irk),sig%spin_index(k))
          enddo
        endif
        call timing%stop(timing%input_write)
      else
        ! FHJ: the following lines were introduced in r6294 and are supposed to
        ! be a shortcut if we are past the last band of the last k-point. However,
        ! in light of a previous bug (#223), this feature is commented out for now.
        !! FHJ: shortcut if this is past the last band of the last k-point
        !if (irk==kp%nrk) exit
      endif
    enddo ! i (loop over bands) ZL: kp%mnband
    if(allocated(zc))then;deallocate(zc);endif
    ! ZL: this is the needed kpoint for wfnkmpi, copy from wfnk to wfnkmpi
    if((istore.eq.1).and.(sig%nkn.gt.1)) then ! ZL: here it stores wfnkmpi
      call timing%start(timing%input_write)
      ikn=iknstore
      wfnkmpi%nkptotal(ikn)=kp%ngk(irk)
      wfnkmpi%isort(1:kp%ngk(irk),ikn)=wfnk%isrtk(1:kp%ngk(irk))
      wfnkmpi%qk(1:3,ikn)=qk(1:3)
      wfnkmpi%el(1:sig%ntband,1:sig%nspin,ikn)= &
        wfnk%ek(1:sig%ntband,1:sig%nspin)
      wfnkmpi%elda(1:sig%ntband,1:sig%nspin,ikn)= &
        wfnk%elda(1:sig%ntband,1:sig%nspin)
      do k=1,sig%nspin*kp%nspinor
        wfnkmpi%cg(1:wfnk%ndv,k,ikn)=wfnk%zk(1:wfnk%ndv,k)
      enddo
      call timing%stop(timing%input_write)
      ! ZL: if sig%nkn.gt.1, clean wfnk%zk, otherwise (only one kpoint in outer) we keep it
      if(associated(wfnk%zk))then;deallocate(wfnk%zk);nullify(wfnk%zk);endif
    endif
    ! ZL: for do_phonq
    if(do_phonq.and.(istore_ep.eq.1).and.(sig%nkn.gt.1)) then
      call timing%start(timing%input_write)
      ikn=iknstore_ep
      wfnk_phonq_mpi%nkptotal(ikn)=kp%ngk(irk)
      wfnk_phonq_mpi%isort(1:kp%ngk(irk),ikn)=wfnk_phonq%isrtk(1:kp%ngk(irk))
      wfnk_phonq_mpi%qk(1:3,ikn)=wfnk_phonq%k(:)
      wfnk_phonq_mpi%el(1:sig%ntband,1:sig%nspin,ikn)= &
        wfnk_phonq%ek(1:sig%ntband,1:sig%nspin)
      wfnk_phonq_mpi%elda(1:sig%ntband,1:sig%nspin,ikn)= &
        wfnk_phonq%elda(1:sig%ntband,1:sig%nspin)
      do k=1,sig%nspin*kp%nspinor
        wfnk_phonq_mpi%cg(1:wfnk_phonq%ndv,k,ikn)=wfnk_phonq%zk(1:wfnk_phonq%ndv,k)
      enddo
      call timing%stop(timing%input_write)
      ! ZL: if sig%nkn.gt.1, remove wfnk%zk, otherwise (only one kpoint in outer) we keep it
      if(associated(wfnk_phonq%zk))then;deallocate(wfnk_phonq%zk);nullify(wfnk_phonq%zk);endif
    endif
    istart = istart + kp%ngk(irk)
  enddo ! irk (loop over k-points)
  call progress_free(prog_info)
  if(allocated(isort))then;deallocate(isort);endif
 
end subroutine read_wavefunctions
end module input_m
