!==============================================================================
!
! Routines:
!
! (1) epsilon Originally by (MSH) Last Edited 5/1/2008 (JRD)
!
! Send comments/bugs to jdeslip@berkeley.edu
!
! See README file for more details
!
!==============================================================================
program epsilon
  use algos_epsilon_m
  use blas_m
  use chi_convergence_m
  use chi_summation_m
  use epsinv_m
  use epsread_hdf5_m
  use epswrite_hdf5_m
  use fftw_m
  use fullbz_m
  use genwf_eps_m
  use global_m
  use gmap_m
  use input_m
  use input_q_m
  use input_utils_m
  use io_utils_m
  use irrbz_m
  use mtxel_m
  use mtxelmultiply_m
  use read_matrix_m
  use references_m
  use rqstar_m
  use scalapack_m
  use sort_m
  use subgrp_m
  use tile_m
  use timing_m, only: common_timing, timing => epsilon_timing
  use vcoul_generator_m
  use write_matrix_m
  use write_program_header_m
  implicit none
  type (kpoints) :: kp
  type (kpoints) :: kpq
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (crystal) :: crys
  type (polarizability) :: pol
  type (valence_wfns) :: vwfn
  type (conduction_wfns) :: cwfn
  type (scalapack) :: scal
  type (int_wavefunction) :: intwfnv
  type (int_wavefunction) :: intwfnvq
  type (int_wavefunction) :: intwfnc
  type(chi_summator_t) :: chi_summator
  type(chi_converger_t) :: chi_converger
  type(wfn_FFT_comm_t) :: wfn_FFT_comm
!-----------------------
! Arrays for kpoints (fullbz, ...etc)
  integer :: nrk
  integer :: indst(48)
  integer, allocatable :: indrk(:),neq(:)
  type(grid) :: gr
!-----------------------
! Arrays for polarizability matrix
  integer :: nstar
  logical :: is_q0, use_WFNq, qpt_done
  integer :: npcmax,nprmax,ivin,neqmax
  integer, allocatable :: ind(:),indt(:,:,:)
  integer, allocatable :: nprdtemp(:),npcdtemp(:)
  real(DP) :: qq(3) !< The current q-pt under consideration
  real(DP) :: omega_plasma, kfact
  real(DP), allocatable :: ekin(:), kweights(:), kvols(:)
  real(DP), allocatable :: ph(:)
  real(DP), allocatable :: pht(:,:,:)
  integer, allocatable :: nst(:)
!-----------------------
! Variables used with HDF5
  integer :: my_iq ! If this is not a q->0 point, write to (iq-nq0) q-point
  character(len=13) :: filename_chi_hdf5='', filename_eps_hdf5='', filename_out_hdf5=''
  character :: tmpstr*100
  character :: filename*13
  integer :: initial_access = 0
  integer :: i,j,k,n,irk,iv,itran,it
  integer :: ncount,ix,jx,kgq(3)
  integer :: iunit,iq
  integer :: ig_l, ig_g, ipe
  integer :: ia, ib, ic, id, ie, if
  integer :: nmtx_t
  real(DP) :: tsec(2)
  real(DP) :: fact,rk(3)
  integer :: ispin
  character*24 :: routnam(120)
  integer, allocatable :: routsrt(:)
  integer :: Ntime
  integer :: iunit_v, iunit_c
  real(DP) :: dnpcr,dnpcrmax
  real(DP) :: E_rpa
  integer :: jj
  type(progress_info) :: prog_info !< a user-friendly progress report
!--------------- Begin Program ---------------------------------------
  call peinfo_init()
!----------------------
! Initialize random numbers
  peinf%jobtypeeval = 0
!--------------------
! Initialize timer
  call timing%init()
  call common_timing%init()
  call timing%start(timing%total)
!------------------------
! Initialize files
  call timing%start(timing%job_setup)
  call open_file(55,file='epsilon.inp',form='formatted',status='old')
  if(peinf%inode .eq. 0) then
    call open_file(7,file='epsilon.log',form='formatted',status='replace')
  endif
  call write_program_header('Epsilon', .true.)
!----------- Call Input: Read crystal data from unit 25 ---------------
! read parameters and q-points from unit 55 (input file)
! initialize unit 10 (output for polarizability matrix)
  call timing%stop(timing%job_setup)
  call timing%start(timing%input)
  call input(kp,crys,syms,gvec,pol,cwfn,vwfn,intwfnv,intwfnc,omega_plasma,gr)
  ! FHJ: consistency check
  if (pol%os_opt_ffts/=0 .and. (kp%nspin*kp%nspinor/=1 .or. pol%ncrit/=0)) then
    call die('Cannot use os_opt_fft/=0 for metals or spin-polarized calculations', &
      only_root_writes=.true.)
  endif
  ! ZL: add warning for broken time reversal symmetry
  if (kp%nspin*kp%nspinor/=1) then
    if (peinf%inode .eq. 0) then
      write(6,*) '==================================== WARNING ===================================='
      write(6,*) 'You have enabled spin polarization or full spinor calculations.'
      write(6,*) 'Current formalism and implementation adopted by BerkeleyGW implicitly assumes'
      write(6,*) 'time-reversal symmetry (TRS), e.g. in the standard Adler-Wiser RPA polarizability.'
      write(6,*) 'Running magnetic/broken-TRS systems is fundamentally against these assumptions.'
      write(6,*) 'One needs to be aware of and justify its direct usage. Use at your own risk.'
      write(6,*) '================================================================================='
    endif
  endif
  !
  if (kp%nspinor == 2) then
    call require_reference(REF_Barker2018)
  endif
  !
  pol%FFTgrid = gvec%FFTgrid
  if (pol%min_fftgrid) then
    if (pol%skip_chi) then
      ! FHJ: We don`t have WFNs, so we can`t do FFTs.
      ! Let`s set it to zero to avoid developers` errors.
      pol%FFTgrid(:) = 0
    else
      ! FHJ: Figure our the FFT grid that holds the WFNs
      call get_wfn_fftgrid(pol, gvec, kp, intwfnv)
    endif
  endif
  if(.not. pol%skip_chi .and. peinf%inode == 0) then
    call open_file(17,file='chi_converge.dat',form='formatted',status='replace')
  endif
  if (pol%iwritecoul .eq. 1) then
    if (peinf%inode .eq. 0) then
      call open_file(19,file='vcoul',form='formatted',status='replace')
    endif
  endif
! CHP: saves the (G=0,G`=0) component of (retarded) epsilon inverse
  if(peinf%inode==0 .and. pol%freq_dep>0 .and. .not.pol%skip_epsilon) then
    call open_file(51,file='EpsInvDyn',form='formatted',status='replace')
    call open_file(52,file='EpsDyn',form='formatted',status='replace')
  endif
  allocate(vwfn%isort (gvec%ng))
  allocate(cwfn%isort (gvec%ng))
  if (pol%nq0>0) then
    ! FHJ: no q->0 point can have all coordinates set to zero
    if (pol%icutv==TRUNC_NONE .and. any(all(abs(pol%qpt(:,:pol%nq0))<TOL_ZERO,dim=1))) then
      call die('No truncation and zero q0', only_root_writes=.true.)
    endif
  endif
  call timing%stop(timing%input)
!------------------------
! Initialize GPU
  call output_algos()
  if (chi_summation_algo == OPENACC_ALGO .or. mtxel_algo == OPENACC_ALGO) then
    call initialize_gpu(OPENACC_ALGO)
  end if
  if (chi_summation_algo == OMP_TARGET_ALGO .or. mtxel_algo == OMP_TARGET_ALGO) then
    call initialize_gpu(OMP_TARGET_ALGO)
  end if
!-------------- Read wavefunctions for (k+q) points ---------------------
! SIB: The input_q routine looks the same as the input routine but
! if reads from a file called WFNq instead of WFN. Presumably
! these are the shifted (by "q") wave functions.
  call timing%start(timing%input_q)
  if (pol%need_WFNq) then
    if (peinf%inode .eq. 0) then
      write(6,*) 'You have a slightly shifted q0 vector and a semiconductor.'
      write(6,*) 'So, reading from WFNq.'
    endif
    call input_q(gvec,kpq,cwfn,vwfn,pol,intwfnvq)
  elseif (pol%lin_denominator>TOL_Zero) then
  endif
  call timing%stop(timing%input_q)
!-------------- GENERATE FULL BZ ----------------------------------------
! SIB: fullbz() takes the kpoints in kp%components(1:3,kp%nrk) and applies all
! the symmetries in syms to them. The resulting set of unique vectors
! are in gr%f(1:3,gr%nf) (gr%nf of them).
  if (peinf%inode .eq. 0) then
    write(6,'(1x,a)') 'Summary of the WFN files:'
    write(6,'(1x,a,i0)') '- Number of k-points in WFN: ', kp%nrk
    if (pol%need_WFNq) then
      write(6,'(1x,a,i0)') '- Number of k-points in WFNq: ', kpq%nrk
    endif
    write(6,'(1x,a,i0)') '- Number of k-points in the full BZ of WFN: ', gr%nf
    if (peinf%verb_high) then
      write(6,'(1x,a,i0,a)') '- Full list of k-points generated with ',syms%ntran,' symmetries:'
      write(6,'(/7x,a,6x)') 'k-point rk (irr. BZ)'
      write(6,'(1x,32("-"))')
      write(6,'(3(1x,f10.6))') (kp%rk(:,iq), iq=1, kp%nrk)
      write(6,'(/7x,a,6x)') 'k-point fk (full BZ)'
      write(6,'(1x,32("-"))')
      write(6,'(3(1x,f10.6))') (gr%f(:,iq), iq=1, gr%nf)
    endif
    write(6,'()')
  endif
  allocate(ekin (gvec%ng))
  allocate(scal%nprd (peinf%npes_freqgrp))
  allocate(scal%npcd (peinf%npes_freqgrp))
  if (pol%os_opt_ffts==2) then
    ! FHJ: FFTs opt. level 2 => do all the FFTs using all the processor, int_wfn arrays
    call genwf_FFT_Isend(wfn_FFT_comm,crys,gvec,syms,kp,kpq,vwfn,pol,cwfn,intwfnv,intwfnvq,intwfnc)
    !call genwf_FFT(crys,gvec,syms,kp,kpq,vwfn,pol,cwfn,intwfnv,intwfnvq,intwfnc,need_WFNq)
  endif
!----------- LOOP over q points for which chi and eps are calculated -----
  do iq=1,pol%nq
    if (pol%stop_after_qpt>-1 .and. iq>pol%stop_after_qpt) then
      if (peinf%inode==0) write(6,'(/,1x,a,/)') 'stop_after_qpt: emulating a sudden application stop.'
      FLUSH(6)
      FLUSH(0)
      call sleep(1d0)
      stop
    endif
    ! SIB: qq(1:3) is the current q vector under consideration
    qq(:)=pol%qpt(:,iq)
    if(peinf%inode.eq.0) then
      call print_dealing_with(iq, pol%nq, qq, 'q')
    endif
    is_q0 = iq<=pol%nq0
    use_WFNq = (is_q0.and.pol%need_WFNq).or.pol%patched_sampling.or.(pol%qflags(iq)==-1)
    if (peinf%inode==0) then
      if (is_q0) then
        write(6,'(1x,a)') 'This is the special q->0 point.'
      else
        write(6,'(1x,a)') 'This is a regular non-zero q-point.'
      endif
    endif
    my_iq = iq
    if (.not.is_q0) my_iq = iq - pol%nq0
    if (is_q0) then
      filename_eps_hdf5 = 'eps0mat.h5'
      filename_chi_hdf5 = 'chi0mat.h5'
    else
      filename_eps_hdf5 = 'epsmat.h5'
      filename_chi_hdf5 = 'chimat.h5'
    endif
    if (pol%skip_epsilon) then
      filename_out_hdf5 = filename_chi_hdf5 ! Write to/restart chimat file
    else
      filename_out_hdf5 = filename_eps_hdf5 ! Write to/restart epsmat file
    endif
!--------------------
! Sort kinetic energies and determine number of matrix elements
! Index of ordered kinetic energies in array isrtx
! SIB: pol%isrtx has the indices for sorted ekin
    call timing%start(timing%q_loop_setup)
    allocate(pol%isrtx (gvec%ng))
    allocate(pol%isrtxi (gvec%ng))
      ! FHJ: Need to (re)compute kinetic energies and sorting arays
      call timing%start(timing%gvec)
      ! Calculate kinetic energies |q+G|^2
      if (is_q0) then
        call kinetic_energies(gvec, crys%bdot, ekin)
      else
        call kinetic_energies(gvec, crys%bdot, ekin, qvec = qq)
      endif
      call logit('sorting gvec')
      call sortrx(gvec%ng, ekin, pol%isrtx, gvec = gvec%components)
      call timing%stop(timing%gvec)
      ! Compute inverse array to isrtx
      do i=1,gvec%ng
        pol%isrtxi(pol%isrtx(i))=i
      enddo
    call timing%stop(timing%q_loop_setup)
! SIB: pol%nmtx becomes the number of matrix elements to be computed;
! the matrix is computed if its ekin is < ecuts
    call timing%start(timing%init_cutoff)
    pol%nmtx = gcutoff(gvec%ng, ekin, pol%isrtx, pol%ecuts)
    if(peinf%inode.eq.0) then
      write(6, '(1x,a,i0)') 'Rank of the polarizability matrix (nmtx): ', pol%nmtx
    endif
    ! FHJ: Do we want to use the economical fftgrid/box? If so, we pad the WFN FFTbox
    ! by the box that holds nmtx gvectors, which is much smaller than the WFN fftbox.
    if (pol%min_fftgrid.and.pol%os_opt_ffts<2) call get_eps_fftgrid(pol, gvec)
    call timing%stop(timing%init_cutoff)
    call timing%start(timing%init_scalapack)
    allocate(pol%irow (pol%nmtx))
    pol%irow=0
! JRD: Determine size of distributed matrices
! DVF : when using parallel frequencies and the number of processors is not divisible
! by the number of frequencies done in parallel, there are excess processors that we
! don`t want to include in our distributed matrix operations. So, for these processors
! we call blacs_setup and then reset the values of the ScaLAPACK variables
! to harmless values that won`t affect any of the global variables obatined via MPI_Allreduce.
    if(pol%nfreq_group .gt. 1) then
      if (peinf%inode .lt. peinf%npes) then
        call blacs_setup(scal, pol%nmtx, .false.,peinf%npes_freqgrp,pol%nfreq_group)
      else
        call blacs_setup(scal, pol%nmtx, .false.,peinf%npes_freqgrp,pol%nfreq_group,peinf%npes_orig-peinf%npes)
        scal%npr=0
        scal%npc=0
        scal%nbl=1
        scal%nprow=1
        scal%npcol=1
        scal%myprow=1000000 ! DVF: A very large number so that we never take anything from
        scal%mypcol=1000000 ! these processors in the garbage frequency/mtxel groups
                             ! See where these variables are used in epsinv.f90 to see
                             ! what I mean.
      endif
    else
      call blacs_setup(scal, pol%nmtx, .true.)
    endif
    scal%nprd = scal%npr
    scal%npcd = scal%npc
! DVF: Get the maximum number of columns/rows owned by one of the processors. This is so you can
! allocate arrays of the right size. For what determines how many rows and columns a procesor
! has, see blacs_setup routine in Common directory and google the numroc routine of scalapack
! Numroc is a nifty little routine
    dnpcr = scal%npc*1D0*scal%npr
! write(*,*) "before npr allreduce",peinf%inode
    npcmax = scal%npc
    nprmax = scal%npr
    ! FHJ: I think there`s actually a but in this report..
    if (dnpcr>(pol%nmtx*1.5D0*pol%nmtx/(1D0*peinf%npes)) .and. &
      peinf%inode==0 .and. peinf%verb_high) then
      write(6,'(/1x,a)') 'NOTE: ScaLAPACK layout is not well load-balanced.'
      write(6,'(1x,a,f12.0/)') 'Max number of elements owned by one task is:', dnpcr
    endif
! Allocate scalapack arrays needed later. See scalapack in Common/scalapack.f90 to see what
! these arrays hold
    allocate(scal%isrtxcol (scal%npc))
    allocate(scal%isrtxrow (scal%npr))
    allocate(scal%imycol (scal%npc))
    allocate(scal%imyrow (scal%npr))
    allocate(scal%imycolinv (pol%nmtx))
    allocate(scal%imyrowinv (pol%nmtx))
    ! FHJ: FIXME - DVF should be sprayed hard for overwriting peinf%npes!!
    allocate(scal%imycold (npcmax,peinf%npes_orig))
    allocate(scal%imyrowd (nprmax,peinf%npes_orig))
    scal%imycold = 0
    scal%imyrowd = 0
    ! scal%imyrow(ig_l) = indxl2g(ig_l, scal%nbl, scal%myprow, 0, scal%nprow)
    ! ipe = indxg2p(ig_g, scal%nbl, 0, 0, scal%npcol)
    ! scal%imyrowd(ig_l, inode+1) = indxl2g(ig_l, scal%nbl, ipe+1, 0, scal%nprow)
    ! FHJ: FIXME - DVF should be sprayed hard for overwriting peinf%npes!!
    if (peinf%inode<peinf%npes_orig) then
      ! FHJ: FIXME - these indexing arrays are completely useless and should be
      ! removed. Let`s not reinvent BLACS, plz!
      do ig_l = 1, scal%npr
        ig_g = indxl2g(ig_l, scal%nbl, scal%myprow, 0, scal%nprow)
        scal%isrtxrow(ig_l) = pol%isrtx(ig_g)
        scal%imyrow(ig_l) = ig_g
        scal%imyrowd(ig_l, peinf%inode+1) = ig_g
        scal%imyrowinv(ig_g) = ig_l
      enddo
      do ig_l = 1, scal%npc
        ig_g = indxl2g(ig_l, scal%nbl, scal%mypcol, 0, scal%npcol)
        scal%isrtxcol(ig_l) = pol%isrtx(ig_g)
        scal%imycol(ig_l) = ig_g
        scal%imycold(ig_l, peinf%inode+1) = ig_g
        scal%imycolinv(ig_g) = ig_l
      enddo
      !do ig_g = 1, pol%nmtx
      ! ig_l = indxl2g(ig_g, scal%nbl, 0, 0, scal%npcol)
      ! ipcol = indxg2p(ig_g, scal%nbl, 0, 0, scal%npcol)
      ! scal%imycold(ig_l,ipe+1) = ig_g
      ! ipe = indxg2p(ig_g, scal%nbl, 0, 0, scal%nprow)
      ! ig_l = indxl2g(ig_g, scal%nbl, 0, 0, scal%nprow)
      ! scal%imyrowd(ig_l,ipe+1) = ig_g
      !enddo
    endif
    call timing%stop(timing%init_scalapack)
    call timing%start(timing%init_arrays)
!----------------------
! Determine subgroup which leaves qq invariant
!
! SIB: figures out which symmetries acting on qq result in qq + integer
! entries. syms%ntranq is their number, syms%indsub are their indices
! (pointing to syms%mtrx), and syms%kgzero(1:3,:) are the integer entries.
    call timing%start(timing%subgrp)
    call subgrp(qq,syms)
    if (pol%patched_sampling) then
      syms%ntranq = 1
    endif
    call timing%stop(timing%subgrp)
!-----------------------
! Determine independent elements of polarizability matrix
!
! SIB: independent means due to symmetries. This initializes
! the polarizability matrix pol%chi to zero (for independent entries)
! and figure out phases due to symmetries for dependent ones,
! and points dependent ones to the entries they depend on (pol%kxi indices)
! JRD: we don`t do this anymore
! call logit('calling indep')
! call indep(nind,gvec,syms,pol,kp%nspin)
!
! JRD: Testing what if we set pol%kxi to zero
! pol%kxi = 0
! pol%chi = 0D0
! nind=pol%nmtx*(pol%nmtx+1)/2
!
!----------------------
! Reduce the k-points to the irr. bz with respect to q
!
! SIB: figure out k-points in irr. BZ (based on symmetries for current q)
! nrk is # of irr. points, indrk are their indices in the full zone,
! and neq is the number of equiv. points for an irr. point.
! (Full zone vectors come in through gr%f(1:3,1:gr%nf).)
    call timing%start(timing%irrbz)
    allocate(indrk (gr%nf))
    allocate(neq (gr%nf))
    call irrbz(syms,gr%nf,gr%f,nrk,neq,indrk)
    call timing%stop(timing%irrbz)
    neqmax = maxval(neq(1:nrk))
! write(6,*) peinf%inode, 'neqmax', neq(1), neqmax
!---------------------------
! Output points in irr. bz
    if(peinf%inode.eq.0) then
      write(6,'(1x,a,i0)') 'Number of k-points in the irreducible BZ(q) (nrk): ', nrk
      if (peinf%verb_medium) then
        write(6,'(/6x,a,5x,a)') 'k-point rk (irr. BZ)', '#eq/fBZ'
        write(6,'(1x,29("-"),1x,7("-"))')
        write(6,'(3(1x,f9.6),1x,i7)') (gr%f(1:3,indrk(j)), neq(j), j=1,nrk)
      endif
      write(7,70) (qq(i),i=1,3),pol%nmtx,nrk
70 format(/ /,5x,'q=',3f7.4,2x,'nmtx=',i8,2x,'nrk=',i3)
    endif
    if (pol%patched_sampling) then
      fact=4.0d0/(product(kp%kgrid)*crys%celvol*kp%nspin*kp%nspinor)
    else
      fact=4.0d0/(dble(gr%nf)*crys%celvol*kp%nspin*kp%nspinor)
    endif
    if (pol%freq_dep .eq. 0) then
      allocate(pol%chi (scal%npr,scal%npc,kp%nspin))
      pol%chi=0.0d0
    endif
    ! some check for subspace truncation method
    IF(pol%subspace) THEN
      IF(.NOT.(pol%freq_dep==2 .AND. pol%freq_dep_method==2)) THEN
        call die('Subspace truncation implemented only for freq_dep=2 and freq_dep_method=2')
      END IF
      ! set this flag to false for the meantime, this will be used to decide
      ! if regenerate the full chi or proceed in the calculation of epsilon
      ! directly within the subspace
      pol%need_full_chi = .FALSE.
    END IF
    ! allocate pol%chiRDyn later
    IF(.NOT. pol%subspace) THEN
      if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
        allocate(pol%chiRDyn (scal%npr,scal%npc,pol%nfreq_in_group,kp%nspin))
        pol%chiRDyn=(0.0,0.0)
      endif
    END IF
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      allocate(pol%chiRDyn (scal%npr,scal%npc,pol%nfreq_in_group,kp%nspin))
      pol%chiRDyn=(0.0,0.0)
      allocate(pol%chiTDyn (pol%os_nsfreq_para,scal%npr,scal%npc,kp%nspin))
      pol%chiTDyn=(0.0,0.0)
    endif
    if (pol%freq_dep .eq. 3) then
      allocate(pol%chiRDyn (scal%npr,scal%npc,pol%nfreq_in_group,kp%nspin))
      pol%chiRDyn=(0.0,0.0)
    endif
    if (.not. pol%skip_chi) then
!------------------------------------
! SIB: allocate space
! write(6,*) peinf%inode, 'allocating pht', neqmax, pol%nmtx, nrk
      allocate(ind (pol%nmtx))
      allocate(ph (pol%nmtx))
      allocate(pht (pol%nmtx,neqmax,nrk))
      allocate(indt (pol%nmtx,neqmax,nrk))
      ind=0
      allocate(nst (nrk))
! JRD: Possible Memory Hazard. We can speed this up by possibly
! only allocating number of bands on current proc and doing send/recvs
      if(pol%nfreq_group .gt. 1) then
        if(pol%gcomm .eq. -1) then
          allocate(pol%gme (pol%nmtx,peinf%ncownmax,peinf%nvownmax,kp%nspin,nrk,pol%nfreq_group))
!disabled PARALLEL DO collapse(6)
          do ia = 1, pol%nfreq_group
            do ib = 1, nrk
              do ic = 1 , kp%nspin
                do id = 1, peinf%nvownmax
                  do ie = 1, peinf%ncownmax
                    do if = 1, pol%nmtx
                      pol%gme(if,ie,id,ic,ib,ia)=0.0d0
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        else
          allocate(pol%gme (pol%nmtx,peinf%ncownmax,peinf%nvownmax,kp%nspin,nrk,1))
!disabled PARALLEL DO collapse(6)
          do ia = 1, 1
            do ib = 1, nrk
              do ic = 1 , kp%nspin
                do id = 1, peinf%nvownmax
                  do ie = 1, peinf%ncownmax
                    do if = 1, pol%nmtx
                      pol%gme(if,ie,id,ic,ib,ia)=0.0d0
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        endif
      else
        allocate(pol%gme (pol%nmtx,peinf%ncownactual,peinf%nvownactual,kp%nspin,nrk,pol%nfreq_group))
!disabled PARALLEL DO collapse(6)
        do ia = 1, pol%nfreq_group
          do ib = 1, nrk
            do ic = 1 , kp%nspin
              do id = 1, peinf%nvownactual
                do ie = 1, peinf%ncownactual
                  do if = 1, pol%nmtx
                    pol%gme(if,ie,id,ic,ib,ia)=0.0d0
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
      if (pol%freq_dep .eq. 2 .or. pol%freq_dep .eq. 3) then
        if(pol%nfreq_group .eq. 1) then
          allocate(pol%edenDyn (peinf%nvownactual,peinf%ncownactual,kp%nspin,nrk,pol%nfreq_group))
        else
          allocate(pol%edenDyn (peinf%nvownmax,peinf%ncownmax,kp%nspin,nrk,pol%nfreq_group))
        endif
      endif
      call timing%stop(timing%init_arrays)
!--------- LOOP OVER K-POINTS IN SET RK ---------------------------------
      ! FHJ: this is to generate nice output / time estimate
      call progress_init(prog_info, 'calculation of matrix elements', 'transition', nrk*peinf%nvownmax)
! SIB: loop over points in irreducible zone
      do irk=1,nrk
        rk(:)=gr%f(:,indrk(irk)) ! rk(:) is the current irr. k-point
! Regenerate star of rk,store the index of the rotation
! SIB: Star is the set of vectors generated by applying all
! subgroup operations for the current q-vector to the k-point rk.
        call timing%start(timing%rqstar)
        call rqstar(syms,nstar,indst,rk)
        if(nstar.ne.neq(irk)) then
          write(0,*) 'nstar?',irk,nstar,neq(irk)
          call die('nstar mismatch')
        endif
        call timing%stop(timing%rqstar)
! JRD: loop over transfs which generate the star of rk for gmap
        nst(irk) = nstar
        call timing%start(timing%gmap)
        do it=1,nstar
! Map g-vectors in polarizability to r**(-1)(g-gq)
! note that gmap requires index of transf in full group
! whereas indst gives index in subgroup
          itran = syms%indsub(indst(it))
          kgq(:) = -syms%kgzero(:,indst(it)) ! note minus sign
          call gmap(gvec,syms,pol%nmtx,itran,kgq,pol%isrtx,pol%isrtxi,ind,ph,.true.)
          pht(:,it,irk) = ph(:)
          indt(:,it,irk) = ind(:)
! debug Statement here
        enddo
        call timing%stop(timing%gmap)
!--------- loop over occupied states -------------------------------------
! SIB: loop over valence states (iv,ispin) where iv is the band index.
        if (pol%os_opt_ffts==2) then
          if (.not.wfn_FFT_comm%done) call genwf_FFT_Wait(wfn_FFT_comm)
          !FHJ: TODO free me later!
          call genwf_lvl2(kp,kpq,vwfn,pol,cwfn)
        endif
        do iv=1,peinf%nvownmax
          ! FHJ : friendly output / running time estimate
          call progress_step(prog_info)
          if (peinf%verb_debug .and. peinf%inode==0) then
            write(6,*) 'Doing iv', iv,'of', peinf%nvownmax
          endif
          if (pol%os_opt_ffts/=2) then
            call genwf_gen(syms,gvec,crys,kp,kpq,irk,rk,qq,vwfn,pol,cwfn,use_WFNq,intwfnv,intwfnvq,intwfnc,iv)
            if (pol%os_opt_ffts==1) then
              ! FHJ: FFT opt. level 1: precalculates FFTs of the conduction bands
              ! each kpt at a time.
              if (iv==1) then
                call mtxel_init_FFT_cond(gvec,pol,cwfn,kp)
              endif
            endif
          endif
! SIB: compute matrix elements and energy denominators for (iv,ispin)
! with all other conduction bands.
          do ispin=1,kp%nspin
            write(tmpstr,'(a,i2,a,i4,a)') "is =", ispin, " iv = ", iv, " calling mtxel"
            call logit(tmpstr)
            if ( iv .le. peinf%nvownactual) then
              call timing%start(timing%mtxel)
              ivin=peinf%invindexv(iv)
              kfact = 1d0
              call mtxel(ivin,gvec,vwfn,cwfn,pol,ispin,irk,kp,kpq,peinf%rank_mtxel,kfact)
              call timing%stop(timing%mtxel)
            endif
          enddo ! ispin
          if (pol%os_opt_ffts<2) then
            ! FHJ: opt. lvl 2 doesn`t even own the WFNs..
            if ( iv .le. peinf%nvownactual) then
              if(associated(vwfn%zv))then;deallocate(vwfn%zv);nullify(vwfn%zv);endif
            endif
          endif
        enddo ! iv
        if (peinf%nvownactual>0) then
          if (pol%os_opt_ffts<2) then
            if(associated(cwfn%zc))then;deallocate(cwfn%zc);nullify(cwfn%zc);endif
          endif
          if (pol%os_opt_ffts==1) then
            ! FHJ: destroy FFTs of conduction bands
            call mtxel_free_FFT_cond(cwfn)
          endif
          if(associated(vwfn%ev))then;deallocate(vwfn%ev);nullify(vwfn%ev);endif
          if(associated(cwfn%ec))then;deallocate(cwfn%ec);nullify(cwfn%ec);endif
        endif
      enddo ! irk
      call progress_free(prog_info)
!------------------------------------------------------------------
! DVF: if requested, test convergence of chi with conduction bands
      if (peinf%inode .eq. 0 .and. pol%freq_dep .eq. 0 .and. pol%fullConvLog .ne. -1) then
        write(6,'(1x,a,i0,a)') 'Preparing simple convergence tests with ', &
          cwfn%nband-vwfn%nband, ' unoccupied bands.'
      endif
      call timing%start(timing%converge_tests)
      if (pol%freq_dep .eq. 0 .and. pol%fullConvLog .ne. -1) then
        call create_chi_converger(chi_converger,vwfn%nband,cwfn%nband)
        if (peinf%verb_debug .and. peinf%inode==0) then
          write(6,'(/,1x,"Starting Convergence Tests")')
        endif
        call chi_convergence_test(pol,pht,indt,kp,nrk,nst,vwfn%nband,cwfn%nband,fact,chi_converger)
        if(peinf%inode .eq. 0) then
          call chi_convergence_print(pol,iq,vwfn%nband,cwfn%nband,chi_converger)
        endif
        call free_chi_converger(chi_converger)
! write(6,*) 'End Convergence Writing'
      endif ! pol%freq_dep .eq. 0
      call timing%stop(timing%converge_tests)
!-----------------------------------------------------------------------
! Construct part of chi that this proc owns by summing the pol%gme matrices
      if (peinf%verb_debug .and. peinf%inode==0) write(6,'(/,1x,"Doing chi Summation")')
      call timing%start(timing%chi_sum_total)
      call create_chi_summator(chi_summator, pol, scal, fact, kp%nspin)
        if (pol%gcomm .eq. 0) then
          call chi_summation_comm_elements(chi_summator,&
                                       pol,scal,kp,vwfn,cwfn,&
                                       nst,nrk,indt,pht)
        else
          call chi_summation_comm_matrix(chi_summator,&
                                         pol,scal,kp,kpq,vwfn,&
                                         nst,nrk,indt,pht)
        endif
      call free_chi_summator(chi_summator, pol)
      call timing%stop(timing%chi_sum_total)
      if (peinf%verb_debug .and. peinf%inode==0) write(6,'(1x,a)') "Done Polarizability"
! Done ChiSum
!-----------------------------------------------------------------------
! Deallocate some arrays no longer needed
      if(allocated(pht))then;deallocate(pht);endif
      if(allocated(indt))then;deallocate(indt);endif
      if(allocated(ind))then;deallocate(ind);endif
      if(allocated(ph))then;deallocate(ph);endif
      if(allocated(nst))then;deallocate(nst);endif
      if(associated(pol%gme))then;deallocate(pol%gme);nullify(pol%gme);endif
      if (pol%freq_dep .eq. 2 .or. pol%freq_dep .eq. 3) then
        if(associated(pol%edenDyn))then;deallocate(pol%edenDyn);nullify(pol%edenDyn);endif
      endif
    else ! pol%skip_chi
!DVF: read chi from chi if this is specified
      if (peinf%inode==0) write(6,'(/1x,a)') 'Reading polarizability matrix from file'
      if (is_q0) then
        iunit=10
        if(peinf%inode.eq.0) then
          write(6,'(1x,a)') 'Reading from file chi0mat.'
          call open_file(unit=10,file='chi0mat',form='unformatted',status='old')
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10) nmtx_t
        endif
        do ispin=1,kp%nspin
          if (pol%freq_dep .eq. 0) then
            call read_matrix_d(scal,pol%chi(:,:,ispin),pol%nmtx,iunit)
          endif ! pol%freq_dep .eq. 0
          if (pol%freq_dep .eq. 2) then
            call read_matrix_f(scal, pol%nFreq, pol%nfreq_in_group, &
              pol%chiRDyn(:,:,:,ispin), pol%nmtx, pol%nfreq_group, iunit)
          endif ! pol%freq_dep .eq. 2!
        enddo ! ispin
      else ! is_q0
        iunit=11
        if(peinf%inode.eq.0) then
          write(6,'(1x,a)') 'Reading from file chimat.'
          if (initial_access .eq. 0) then
            call open_file(unit=11,file='chimat',form='unformatted',status='old')
            read(11)
            read(11)
            read(11)
            read(11)
            read(11)
            read(11)
            read(11)
            read(11)
            read(11)
            read(11)
          endif
          read(11)
          read(11) nmtx_t
! write(6,*) 'nmtx_t for chimat', nmtx_t
        endif
        do ispin=1,kp%nspin
          if (pol%freq_dep .eq. 0) then
            call read_matrix_d(scal,pol%chi(:,:,ispin),pol%nmtx,iunit)
          endif ! pol%freq_dep .eq. 0
          if (pol%freq_dep .eq. 2) then
            call read_matrix_f(scal, pol%nFreq, pol%nfreq_in_group, &
              pol%chiRDyn(:,:,:,ispin), pol%nmtx, pol%nfreq_group, iunit)
          endif ! pol%freq_dep .eq. 2
        enddo ! ispin
        initial_access = 1
      endif ! is_q0
      if (peinf%inode==0) write(6,'(1x,a/)') 'Ok'
    endif ! pol%skip_chi
!-----------------------------------------------------------------------
! JRD: Now write out elements that Proc 1 owns
    do ispin = 1, kp%nspin
      if (pol%freq_dep.eq.0 .and. peinf%inode.eq.0) then
        write(7,940) ispin,kp%nspin
        do i=1,scal%npr
          ix=scal%isrtxrow(i)
          do j=1,scal%npc
! JRD: Diagonal, subdiagonal and wings only
            jx=scal%isrtxcol(j)
            if(i.eq.j .or. (gvec%components(1,ix) .eq. 0 .and. gvec%components(2,ix) .eq. 0 .and. gvec%components(3,ix) .eq. 0)) &
              write(7,950) (gvec%components(k,ix),k=1,3),ekin(ix),(gvec%components(k,jx),k=1,3),ekin(jx), &
                            pol%chi(i,j,ispin)
          enddo
        enddo
      endif ! pol%freq_dep.eq.0 .and. peinf%inode.eq.0
      if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2) .and. peinf%inode.eq.0) then
        IF(pol%subspace .AND. (.NOT.pol%need_full_chi) ) THEN
          ! write only for omega = 0
          ! write(7,940) ispin,kp%nspin
          IF(ispin == 1) THEN
            write(7,*)'frq=',0
            do i=1,scal%npr
              ix=scal%isrtxrow(i)
              do j=1,scal%npc
                jx=scal%isrtxcol(j)
                if(i.eq.j) &
                  write(7,950) (gvec%components(k,ix),k=1,3),ekin(ix),(gvec%components(k,jx),k=1,3),ekin(jx), &
                                pol%chiRDyn_sym_omega0(i,j)
              enddo
            enddo
          END IF
          WRITE(7,*)
          WRITE(7,*) 'Eigenvalues symmetrized chi at omega = 0'
          DO i = 1, pol%nmtx
            WRITE(7,*) i, pol%eigenval_omega0(i)
          END DO
          WRITE(7,*)
        ELSE
          write(7,940) ispin,kp%nspin
          do jj=1,pol%nfreq_in_group
            write(7,*)'frq=',jj
! JRD XXX the i and j loops are out of order here....
            do i=1,scal%npr
              ix=scal%isrtxrow(i)
              do j=1,scal%npc
 !!!JRD: Diagonal and subdiagonal only
                jx=scal%isrtxcol(j)
                if(i.eq.j) &
                  write(7,950) (gvec%components(k,ix),k=1,3),ekin(ix),(gvec%components(k,jx),k=1,3),ekin(jx), &
                                pol%chiRDyn(i,j,jj,ispin)
              enddo
            enddo
          enddo
        END IF ! IF(pol%subspace
      endif ! (pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2) .and. peinf%inode.eq.0
      if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1).and. peinf%inode.eq.0) then
        write(7,940) ispin,kp%nspin
         do jj=1,pol%nfreq_in_group
           write(7,*) 'frq=',jj
! JRD XXX the i and j loops are out of order here
           do i=1,scal%npr
            ix=scal%isrtxrow(i)
            do j=1,scal%npc
              jx=scal%isrtxcol(j)
              if(i.eq.j) &
                write(7,950) (gvec%components(k,ix),k=1,3),ekin(ix),(gvec%components(k,jx),k=1,3),ekin(jx), &
                              pol%chiRDyn(i,j,jj,ispin)
            enddo
          enddo
        enddo
      endif ! (pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1) .and. peinf%inode.eq.0
940 format(/,10x,' independent matrix elements of chi', 7x,'spin index= ',1i1,1x,1i1,/,/,&
        10x,'g',10x,'g**2',10x,'gp',10x,'gp**2',10x,'chi(g,gp)')
! if last value is real, only one of the f13.8 will be used.
950 format(3i5,f10.5,5x,3i5,f10.5,5x,2f15.10)
    enddo ! ispin (loop over spins)
! write(6,*) 'End Element Writing'
!--------- write polarizability matrix and crystal info to file ---------
    if (pol%skip_epsilon) call write_chi(trim(filename_chi_hdf5), iq, my_iq, &
      qq, crys, gvec, pol, ekin, scal, is_q0)
! Use pol%chi(j,1) as sum over spin components
! JRD: Why was proc 0 the only one doing this?!
    if (pol%freq_dep .eq. 0) then
      if(kp%nspin.eq.2) pol%chi(:,:,1)=pol%chi(:,:,1)+pol%chi(:,:,2)
    endif ! pol%freq_dep .eq. 0
    if (pol%freq_dep .eq. 2) then
      if(kp%nspin.eq.2) pol%chiRDyn(:,:,:,1)=pol%chiRDyn(:,:,:,1)+pol%chiRDyn(:,:,:,2)
    endif ! pol%freq_dep .eq. 2
    if (.not. pol%skip_epsilon) then
      call timing%start(timing%epsinv_total)
      call logit('Calling epsinv')
      if(pol%do_rpa) then
        call epsinv(gvec,pol,ekin,qq,is_q0,crys,scal,kp,omega_plasma,iq,E_rpa)
        pol%E_rpa_qp(iq) = E_rpa
      else
        call epsinv(gvec,pol,ekin,qq,is_q0,crys,scal,kp,omega_plasma,iq)
      endif
      if (pol%subspace .and. (.not.pol%need_full_chi)) then
        if(allocated(pol%chiRDyn_sym_omega0))then;deallocate(pol%chiRDyn_sym_omega0);endif
        if(allocated(pol%eigenvect_omega0))then;deallocate(pol%eigenvect_omega0);endif
        if(allocated(pol%eigenval_omega0))then;deallocate(pol%eigenval_omega0);endif
        if(allocated(pol%vcoul_sub))then;deallocate(pol%vcoul_sub);endif
      endif
      call logit('Finished epsinv')
      call timing%stop(timing%epsinv_total)
    endif ! pol%skip_epsilon
    if (pol%freq_dep .eq. 0) then
      if(associated(pol%chi))then;deallocate(pol%chi);nullify(pol%chi);endif
    endif
    if (pol%freq_dep .eq. 2) then
      if(associated(pol%chiRDyn))then;deallocate(pol%chiRDyn);nullify(pol%chiRDyn);endif
      if(pol%freq_dep_method .eq. 1) then
        if(associated(pol%chiTDyn))then;deallocate(pol%chiTDyn);nullify(pol%chiTDyn);endif
      endif
    endif
    if (pol%freq_dep .eq. 3) then
      if(associated(pol%chiRDyn))then;deallocate(pol%chiRDyn);nullify(pol%chiRDyn);endif
    endif
    if(allocated(indrk))then;deallocate(indrk);endif
    if(allocated(neq))then;deallocate(neq);endif
    if(associated(pol%isrtx))then;deallocate(pol%isrtx);nullify(pol%isrtx);endif
    if(associated(pol%isrtxi))then;deallocate(pol%isrtxi);nullify(pol%isrtxi);endif
    if(associated(pol%irow))then;deallocate(pol%irow);nullify(pol%irow);endif
    if(associated(scal%isrtxcol))then;deallocate(scal%isrtxcol);nullify(scal%isrtxcol);endif
    if(associated(scal%isrtxrow))then;deallocate(scal%isrtxrow);nullify(scal%isrtxrow);endif
    if(associated(scal%imycol))then;deallocate(scal%imycol);nullify(scal%imycol);endif
    if(associated(scal%imyrow))then;deallocate(scal%imyrow);nullify(scal%imyrow);endif
    if(associated(scal%imycolinv))then;deallocate(scal%imycolinv);nullify(scal%imycolinv);endif
    if(associated(scal%imyrowinv))then;deallocate(scal%imyrowinv);nullify(scal%imyrowinv);endif
    if(associated(scal%imycold))then;deallocate(scal%imycold);nullify(scal%imycold);endif
    if(associated(scal%imyrowd))then;deallocate(scal%imyrowd);nullify(scal%imyrowd);endif
  enddo ! iq (loop over q points)
  if(associated(scal%nprd))then;deallocate(scal%nprd);nullify(scal%nprd);endif
  if(associated(scal%npcd))then;deallocate(scal%npcd);nullify(scal%npcd);endif
  call dealloc_grid(gr)
! End q point loop!
!-------------------------------------------------------------------
! XXXXXXXXXXXXXXXX Do BZ sum to get RPA correlation energy
  if(pol%do_rpa) then
    E_rpa = 0.0D+00
    do iq = 1, pol%nq
      E_rpa = E_rpa + pol%qw_rpa(iq) * pol%E_rpa_qp(iq)
    enddo
    if(peinf%inode.eq.0) then
      write(*,*) "RPA energy (Ry) = ", E_rpa
      write(*,*) "RPA energy (Ha) = ", E_rpa/2.0D+00
      write(*,*) "RPA energy (eV) = ", E_rpa*ryd
    endif
  endif
! XXXXXXXXXXXXXXXX
!----------- Clean House -------------------------------------------
  call logit('Cleaning up')
  if (.not. pol%skip_epsilon) then
    call destroy_qran()
  endif
  if(.not. pol%skip_chi) then
    if(peinf%inode == 0) call close_file(17) ! file chi_converge.dat
    call destroy_fftw_plans()
  endif
  if (pol%iwritecoul .eq. 1) then
    if (peinf%inode .eq. 0) then
      call close_file(19) ! file vcoul
    endif
  endif
  if (peinf%inode==0 .and. pol%freq_dep==2 .and. .not.pol%skip_epsilon) then
    call close_file(51) !file EpsInvDyn
    call close_file(52) !file EpsDyn
  endif
  if(allocated(ekin))then;deallocate(ekin);endif
  call kp%free()
  if(pol%need_WFNq) then
    call kpq%free()
  endif
  if(associated(gvec%components))then;deallocate(gvec%components);nullify(gvec%components);endif
  if(associated(gvec%index_vec))then;deallocate(gvec%index_vec);nullify(gvec%index_vec);endif
  if(associated(pol%qpt))then;deallocate(pol%qpt);nullify(pol%qpt);endif
  if(associated(vwfn%isort))then;deallocate(vwfn%isort);nullify(vwfn%isort);endif
  if(associated(cwfn%isort))then;deallocate(cwfn%isort);nullify(cwfn%isort);endif
  if(associated(peinf%global_nvown))then;deallocate(peinf%global_nvown);nullify(peinf%global_nvown);endif
  if(associated(peinf%global_ncown))then;deallocate(peinf%global_ncown);nullify(peinf%global_ncown);endif
  if(associated(peinf%indexc))then;deallocate(peinf%indexc);nullify(peinf%indexc);endif
  if(associated(peinf%indexv))then;deallocate(peinf%indexv);nullify(peinf%indexv);endif
  if(associated(peinf%global_indexv))then;deallocate(peinf%global_indexv);nullify(peinf%global_indexv);endif
  if(associated(peinf%invindexv))then;deallocate(peinf%invindexv);nullify(peinf%invindexv);endif
  if(associated(peinf%invindexc))then;deallocate(peinf%invindexc);nullify(peinf%invindexc);endif
  if(associated(peinf%doiownv))then;deallocate(peinf%doiownv);nullify(peinf%doiownv);endif
  if(associated(peinf%doiownc))then;deallocate(peinf%doiownc);nullify(peinf%doiownc);endif
  if(associated(peinf%does_it_ownv))then;deallocate(peinf%does_it_ownv);nullify(peinf%does_it_ownv);endif
  if(associated(peinf%does_it_ownc))then;deallocate(peinf%does_it_ownc);nullify(peinf%does_it_ownc);endif
  if(associated(peinf%global_pairowner))then;deallocate(peinf%global_pairowner);nullify(peinf%global_pairowner);endif
  if(associated(pol%dFreqGrid))then;deallocate(pol%dFreqGrid);nullify(pol%dFreqGrid);endif
  if(associated(pol%dFreqBrd))then;deallocate(pol%dFreqBrd);nullify(pol%dFreqBrd);endif
  if(peinf%inode.eq.0) then
    if(associated(pol%nmtx_of_q))then;deallocate(pol%nmtx_of_q);nullify(pol%nmtx_of_q);endif
    call close_file(7) ! epsilon.log
    if (pol%skip_epsilon.and..not.pol%use_hdf5) then
      if (pol%nq0>0) call close_file(10) ! chi0mat
      if (pol%nq1>0) call close_file(11) ! chimat
    else
      if (.not.pol%use_hdf5) then
        if (pol%nq0>0) call close_file(12) ! eps0mat
        if (pol%nq1>0) call close_file(13) ! epsmat
      endif
    endif ! pol%skip_epsilon
  endif
  if(pol%subspace) then
    if(allocated(pol%neigen_of_q))then;deallocate(pol%neigen_of_q);endif
  end if
  if (.not. pol%skip_chi) then
    call free_wfns(pol, intwfnv, intwfnvq, intwfnc, .true.)
  endif
  call show_references()
!------------- Print Timing Info -----------------------------------------
  call logit('Calculating Timing Info')
  call timing%stop(timing%total)
  call timing%print(common_timing)
  if(allocated(routsrt))then;deallocate(routsrt);endif
  call write_memory_usage()
!-------------------------------
! JIM: Close HDF interface
contains
  !> Computes v(q+G) and writes it to an HDF5 file.
  !! We only use this routine when writing chi, i.e., skipping epsilon, as
  !! epsinv.f90 does this is we are actually writing epsilon inverse.
  subroutine write_vcoul_hdf5(filename, iq, my_iq, qq, crys, gvec, pol)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: iq, my_iq
    real(DP), intent(in) :: qq(3)
    type(crystal), intent(in) :: crys
    type(gspace), intent(in) :: gvec
    type(polarizability), intent(in) :: pol
   
  end subroutine write_vcoul_hdf5
  !> Write the polarizability matrix and v(q+G) to file.
  subroutine write_chi(filename, iq, my_iq, qq, crys, gvec, pol, ekin, scal, is_q0)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: iq, my_iq
    real(DP), intent(in) :: qq(3)
    type(crystal), intent(in) :: crys
    type(gspace), intent(in) :: gvec
    type(polarizability), intent(in) :: pol
    real(DP), intent(in) :: ekin(:)
    type(scalapack), intent(in) :: scal
    logical, intent(in) :: is_q0
    integer :: isorti(gvec%ng), ii, ispin, iunit, np
   
    if (peinf%inode==0) write(6,'(/1x,a)') 'Writing polarizability matrix to file'
      iunit = 11
      if (is_q0) iunit = 10
      if (peinf%inode==0) then
        write(iunit) syms%ntranq, (((syms%mtrx(i,j,syms%indsub(n)),i=1,3),j=1,3), &
          (syms%tnp(k,syms%indsub(n)),syms%kgzero(k,n),k=1,3),n=1,syms%ntranq)
        np = pol%nmtx*(pol%nmtx+1)/2
        write(iunit) pol%nmtx, np, (pol%isrtx(i),ekin(i),i=1,gvec%ng), (pol%irow(i),i=1,pol%nmtx)
      endif
      do ispin = 1, kp%nspin
        select case (pol%freq_dep)
          case (0)
            call write_matrix_d(scal, pol%chi(:,:,ispin), pol%nmtx, iunit)
          case (2,3)
            call write_matrix_f(scal, pol%nFreq, pol%chiRDyn(:,:,:,ispin), &
            pol%nmtx, iunit, pol%nfreq_group)
        endselect
      enddo ! ispin
    if (peinf%inode==0) write(6,'(1x,a/)') 'Ok'
   
  end subroutine write_chi
end program epsilon
