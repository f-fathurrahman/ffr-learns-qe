!===============================================================================
!
! Routines:
!
! (1) input_outer() Originally By gsm Last Modified 10/5/2009 (gsm)
!
! Adapted from subroutine input. Reads in and distributes outer
! wavefunctions. Stores wavefunctions in structures distributed in memory.
!
!===============================================================================

module input_outer_m
  use eqpcor_m
  use global_m
  use input_utils_m
  use io_utils_m
  use misc_m
  use scissors_m
  use wfn_rho_vxc_io_m
  implicit none
  private
  public :: &
    input_outer
contains
subroutine input_outer(crys,gvec,syms,kp,sig,wfnk,iunit_k,fnk,wfnkmpi)
  type (crystal), intent(in) :: crys
  type (symmetry), intent(in) :: syms
  type (gspace), intent(in) :: gvec
  type (kpoints), intent(in) :: kp
  type (siginfo), intent(inout) :: sig
  type (wfnkstates), intent(inout) :: wfnk !< do not lose the allocations from input
  integer, intent(in) :: iunit_k
  character*20, intent(in) :: fnk
  type (wfnkmpiinfo), intent(inout) :: wfnkmpi !< do not lose the allocations from input
  character :: fncor*32
  integer :: i,j,k,ikn,irk,nbnmin,nbnmax, &
    iouter,istore,nstore,nstore_sum,g1,g2,iknstore, &
    ntband_outer
  integer, allocatable :: isort(:)
  real(DP) :: qk(3)
  real(DP), allocatable :: zc(:,:)
  logical :: dont_read
  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvec_outer, gvec_kpt
  type(kpoints) :: kp_outer
  type(crystal) :: crys_outer
  type(symmetry) :: syms_outer
  logical, allocatable :: kpt_outer_required(:)
  type(progress_info) :: prog_info !< a user-friendly progress report
 
  call logit('Opening WFN_outer')
!-------------------------------------------
! Try to open WFN_outer file
! Return to caller if error is encountered
  if(peinf%inode.eq.0) call open_file(25,file='WFN_outer',form='unformatted',status='old',iostat=iouter)
  sig%wfn_outer_present = (iouter == 0)
  if (iouter.ne.0) then
   
    return
  endif ! iouter
  sheader = 'WFN'
  iflavor = 0
  call read_binary_header_type(25, sheader, iflavor, kp_outer, gvec_outer, syms_outer, crys_outer, warn = .false.)
  call check_trunc_kpts(sig%icutv, kp_outer)
  call check_header('WFN_inner', kp, gvec, syms, crys, 'WFN_outer', kp_outer, gvec_outer, syms_outer, crys_outer, is_wfn = .true., tolerant=sig%tolerant_value)
  ! DAS: there is no fundamental reason for the condition below, but seems to be assumed in the code
  ! and it would require some rewriting to make that kp_outer need only contain the kpts in sig%kpt (the ones in the input file).
  if(kp_outer%nrk /= kp%nrk) call die('different number of k-points in WFN_outer')
  if(kp_outer%mnband < maxval(sig%diag)) call die('not enough bands in WFN_outer')
  ntband_outer = min(kp_outer%mnband, sig%ntband)
  kp_outer%el(:,:,:) = kp_outer%el(:,:,:) - sig%avgpot_outer / ryd
  allocate(kp_outer%elda (ntband_outer,kp_outer%nrk,kp_outer%nspin))
  kp_outer%elda(1:ntband_outer,1:kp_outer%nrk,1:kp_outer%nspin)=ryd*kp_outer%el(1:ntband_outer,1:kp_outer%nrk,1:kp_outer%nspin)
  call scissors_shift(kp_outer, sig%scis_outer, sig%spl_tck_outer)
  kp_outer%el(:,:,:)=ryd*kp_outer%el(:,:,:)
!---------------------------------------------------------------
! If quasi-particle corrections requested, read in qp energies
  if(sig%eqp_outer_corrections) then
    fncor='eqp_outer.dat'
    nbnmin=sig%ntband+1
    nbnmax=0
    do i=1,sig%ndiag
      if (sig%diag(i).lt.nbnmin) nbnmin=sig%diag(i)
      if (sig%diag(i).gt.nbnmax) nbnmax=sig%diag(i)
    enddo
    ! FHJ: Determine which k-points from WFN_outer/eqp_outer.dat we need.
    allocate(kpt_outer_required (kp_outer%nrk))
    kpt_outer_required(:) = .false.
    do ikn = 1, sig%nkn
      kpt_outer_required(sig%indkn(ikn)) = .true.
    enddo
    call eqpcor(fncor,peinf%inode,peinf%npes,kp_outer,nbnmin,nbnmax,0,0,&
      kp_outer%el,kp_outer%el,kp_outer%el,0,0,kpt_required=kpt_outer_required)
    if(allocated(kpt_outer_required))then;deallocate(kpt_outer_required);endif
  endif
!-------------------
! Read in g-vectors
  allocate(gvec_outer%components (3, gvec%ng))
  call read_binary_gvectors(25, gvec%ng, gvec%ng, gvec_outer%components)
  if(associated(gvec_outer%components))then;deallocate(gvec_outer%components);nullify(gvec_outer%components);endif
  allocate(isort (gvec%ng))
!---------------------
! Loop over k-points
  call progress_init(prog_info, 'reading wavefunctions (WFN_outer)', 'state', kp_outer%nrk*ntband_outer)
  do irk=1,kp_outer%nrk
    qk(:)=kp_outer%rk(:,irk)
!--------------------------------------------
! Read in and sort g-vectors for k-point irk
    allocate(gvec_kpt%components (3, kp_outer%ngk(irk)))
    call read_binary_gvectors(25, kp_outer%ngk(irk), kp_outer%ngk(irk), gvec_kpt%components)
    do i = 1, kp_outer%ngk(irk)
      call findvector(isort(i), gvec_kpt%components(:, i), gvec)
      if (isort(i) == 0) call die('outer wfn: could not find G-vector')
    enddo
    if(associated(gvec_kpt%components))then;deallocate(gvec_kpt%components);nullify(gvec_kpt%components);endif
!--------------------------------------------------------
! Determine if Sigma must be computed for this k-point.
! If so, store the bands and wavefunctions on file iunit_k.
! If there is only one k-point, store directly in wfnk.
    istore=0
    do ikn=1,sig%nkn
      if(sig%indkn(ikn).eq.irk) then
        istore=1
        iknstore=ikn
        wfnk%nkpt=kp_outer%ngk(irk)
        wfnk%ndv=peinf%ndiag_max*kp_outer%ngk(irk)
        wfnk%k(:)=qk(:)
        wfnk%isrtk(:)=isort(:)
        do k=1,sig%nspin
          wfnk%ek(1:ntband_outer,k)= &
            kp_outer%el(1:ntband_outer,irk,sig%spin_index(k))
          wfnk%elda(1:ntband_outer,k)= &
            kp_outer%elda(1:ntband_outer,irk,sig%spin_index(k))
        enddo
        allocate(wfnk%zk (wfnk%ndv,sig%nspin*kp_outer%nspinor))
        wfnk%zk=0.0d0
      endif
    enddo ! ikn
!---------------------------------------------------------------
! Read the wavefunctions, check the norm, and store in wfnk%zk
    allocate(zc (kp_outer%ngk(irk),kp_outer%nspin*kp_outer%nspinor))
    do i=1,kp_outer%mnband
! Check whether i-th band is needed for diagonal Sigma matrix elements.
! We don`t need to check for off-diagonal Sigma matrix elements separately
! because whenever off-diagonals are requested the corresponding diagonals
! are also included in the calculation as enforced in Sigma/inread.f90.
      if (i.le.ntband_outer) then
        call progress_step(prog_info, ntband_outer*(irk-1) + i)
        nstore=0
        do j=1,peinf%ndiag_max
          if (.not.peinf%flag_diag(j)) cycle
          if (i==sig%diag(peinf%index_diag(j))) nstore=j
        enddo
        nstore_sum=nstore
      else
        nstore_sum=0
      endif
      dont_read = (istore == 0 .OR. nstore_sum == 0)
      call read_binary_data(25, kp_outer%ngk(irk), kp_outer%ngk(irk), kp_outer%nspin*kp_outer%nspinor, zc, dont_read = dont_read)
      if(.NOT. dont_read) then
        if(peinf%inode.eq.0) then
          do k=1,sig%nspin
            call checknorm('WFN_outer',i,irk,kp_outer%ngk(irk),k,kp_outer%nspinor,zc(:,:))
          enddo
        endif
        if((istore.eq.1).and.(nstore.ne.0)) then
          do k=1,sig%nspin*kp_outer%nspinor
            do j=1,kp_outer%ngk(irk)
              wfnk%zk((nstore-1)*kp_outer%ngk(irk)+j,k) = zc(j,sig%spin_index(k))
            enddo
          enddo
        endif
      else
        ! FHJ: the following lines were introduced in r6294 and are supposed to
        ! be a shortcut if we are past the last band of the last k-point. However,
        ! in light of a previous bug (#223), this feature is commented out for now.
        !! FHJ: shortcut if this is past the last band of the last k-point
        !if (irk==kp_outer%nrk) exit
      endif ! .NOT. dont_read
    enddo ! i (loop over bands)
    if(allocated(zc))then;deallocate(zc);endif
!---------------------------------------------------------
! Write the wavefunctions stored in wfnk%zk to file iunit_k
    if((istore.eq.1).and.(sig%nkn.gt.1)) then
      ikn=iknstore
      wfnkmpi%nkptotal(ikn)=kp_outer%ngk(irk)
      wfnkmpi%isort(1:kp_outer%ngk(irk),ikn)=wfnk%isrtk(1:kp_outer%ngk(irk))
      wfnkmpi%qk(1:3,ikn)=qk(1:3)
      wfnkmpi%el(1:ntband_outer,1:sig%nspin,ikn)=wfnk%ek(1:ntband_outer,1:sig%nspin)
      wfnkmpi%elda(1:ntband_outer,1:sig%nspin,ikn)=wfnk%elda(1:ntband_outer,1:sig%nspin)
      do k=1,sig%nspin*kp_outer%nspinor
        wfnkmpi%cg(1:wfnk%ndv,k,ikn)=wfnk%zk(1:wfnk%ndv,k)
      enddo
      if(associated(wfnk%zk))then;deallocate(wfnk%zk);nullify(wfnk%zk);endif
    endif
  enddo ! irk (loop over k-points)
  call progress_free(prog_info)
  if(peinf%inode.eq.0) then
    call close_file(25) ! WFN_outer
  endif
  if(allocated(isort))then;deallocate(isort);endif
  call kp_outer%free()
 
  return
end subroutine input_outer
end module input_outer_m
