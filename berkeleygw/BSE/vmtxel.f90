!==============================================================================
!
! Module vmtxel_m
!
! Originally by GKA (2018)
!
! Objects holding dipole operator.
!
!==============================================================================

module vmtxel_m
  use global_m
  use mtxel_optical_m
  use misc_m, only: bse_index
  implicit none
public
! =========================================================================== !
!> dipole operator (velocity)
type vmtxel_t
  ! Objects
  type (crystal) :: crys
  type (gspace) :: gvec
  ! Dimensions
  integer :: ns !< Number of spins
  integer :: nk !< Number of k-points
  integer :: nband !< Number of bands at k (conduction bands in BSE)
  integer :: mband !< Number of bands at k+q (valence bands in BSE)
  integer :: nmat !< Flat dimension (ns * nk * nband * mband)
  integer :: npol !< Number of polarization
  ! Options
  integer :: band_ordering !< 0 = counting bands from fermi levels, with
                            !! conduction bands going up in energy
                            !! valence bands down in energy.
                            !! 1 = counting all bands from bottom-up,
                            !! starting at the first valence band
                            !! and going up in energy.
  integer :: opr !< 0 = use velocity operator
                  !! 1 = use momentum operator
  logical :: has_velocity = .false. !< Whether the velocity operator has
                                      !! been computed along with the
                                      !! dipole operator
  logical :: use_hdf5 = .true. !< Use hdf5 for I/O
  ! Parallel info
  logical :: is_master
  ! Arrays
  !> The list of k-points
  !> kpt(3,nkpt)
  real(dp), allocatable :: kpt(:,:)
  !> The polarization vector.
  !> pol(3,npol)
  real(dp), allocatable :: pol(:,:)
  !> Dipole operator, for a single k-point, in band index basis.
  !! s1k(nband,mband,ns,npol)
  real(DP), allocatable :: s1k(:,:,:,:)
  !> Dipole operator, for all k-point, in flat BSE indices basis.
  !! s1(nmat,npol)
  real(DP), allocatable :: s1(:,:)
  !> Dipole operator, for all k-point, in band index basis.
  !! At the moment, this is treated as a temporary array for io only.
  !! rnm(nband,mband,nk,ns,npol)
  real(DP), allocatable :: rnm(:,:,:,:,:)
  !> Velocity operator, for all k-point, in band index basis.
  !! vnm(nband,mband,nk,ns,npol)
  real(DP), allocatable :: vnm(:,:,:,:,:)
  contains
  ! Core procedures
  procedure :: init => vmtxel_init
  procedure :: init_from_xctinfo => vmtxel_init_from_xctinfo
  procedure :: alloc => vmtxel_alloc
  procedure :: free => vmtxel_free
  ! Computation / communication
  procedure :: compute_ik_vmtxel
  procedure :: reduce => reduce_vmtxel
  procedure :: band_to_flat => vmtxel_band_to_flat
  procedure :: flat_to_band => vmtxel_flat_to_band
  ! I/O
  procedure :: write_vmtxel
  procedure :: read_vmtxel
  procedure :: write_vmtxel_bin
  procedure :: read_vmtxel_bin
end type vmtxel_t
! =========================================================================== !
contains
! =========================================================================== !
!> Initialize object, setting manually dimensions and options
subroutine vmtxel_init(this, ns, nk, nband, mband, opr, npol, &
                       band_ordering, with_velocity)
  class(vmtxel_t), intent(inout) :: this
  integer, intent(in) :: ns !< Number of spins
  integer, intent(in) :: nk !< Number of k-points
  integer, intent(in) :: nband !< Number of bands at k (conduction)
  integer, intent(in) :: mband !< Number of bands at k+q (valence)
  integer, intent(in),optional :: opr !< 0 = use velocity operator
                                          !! 1 = use momentum operator
  integer, intent(in), optional :: npol !< Number of polarization
  integer, intent(in), optional :: band_ordering !< 0 = from fermi levels
                                                  !! 1 = bottom-up
  logical, intent(in), optional :: with_velocity !< Do compute the velocity
                                                  !! operator
  !real(dp), intent(in), optional :: pol(:,:) !< Polarization vectors
 
  ! Get parallel info
  this%is_master = (peinf%inode.eq.0)
  this%ns = ns
  this%nk = nk
  this%nband = nband
  this%mband = mband
  this%nmat = this%ns * this%nk * this%nband * this%mband
  if (present(opr)) then
    this%opr = opr
  else
    this%opr = 1
  end if
  if (present(npol)) then
    this%npol = npol
  else
    if (this%opr .eq. 1) then
      this%npol = 3
    else
      this%npol = 1
    end if
  end if
  if (present(band_ordering)) then
    this%band_ordering = band_ordering
  else
    this%band_ordering = 0
  end if
  if (present(with_velocity)) then
    this%has_velocity = with_velocity
  else
    this%has_velocity = .false.
  end if
 
end subroutine vmtxel_init
! =========================================================================== !
!> Initialize object from an xctinfo object, copying dimensions and options
!! as well as the polarization vector
subroutine vmtxel_init_from_xctinfo(this, xct, opr)
  class(vmtxel_t), intent(inout) :: this
  type(xctinfo), intent(in) :: xct
  integer, intent(in),optional :: opr !< 0 = use velocity operator
                                          !! 1 = use momentum operator
  integer :: opr_ = 0
 
  if (present(opr)) opr_ = opr
  call this%init(xct%nspin, xct%nkpt_fi, xct%ncb_fi, xct%nvb_fi, &
                 opr=opr_, npol=xct%npol)
  ! Set polarization vector
  if (this%npol .eq. 1) then
    allocate(this%pol (3, this%npol))
    this%pol(:,1) = xct%pol
  end if
 
end subroutine vmtxel_init_from_xctinfo
! =========================================================================== !
!> Allocate arrays
subroutine vmtxel_alloc(this, rnm)
  class(vmtxel_t), intent(inout) :: this
  logical, intent(in), optional :: rnm
  integer :: ipol
 
  allocate(this%kpt (3, this%nk))
  allocate(this%s1 (this%nmat, this%npol))
  allocate(this%s1k (this%nband, this%mband, this%ns, this%npol))
  this%s1 = 0.0d0
  this%s1k = 0.0d0
  if (.not. allocated(this%pol)) then
    allocate(this%pol (3, this%npol))
    do ipol=1,this%npol
      this%pol(:,ipol) = 0.0d0
      this%pol(ipol,ipol) = 1.0d0
    end do
  end if
  if (present(rnm)) then
    if (rnm) then
      allocate(this%rnm (this%nband, this%mband, this%nk, this%ns, this%npol))
    end if
  end if
  if (this%has_velocity) then
    allocate(this%vnm (this%nband, this%mband, this%nk, this%ns, this%npol))
  end if
 
end subroutine vmtxel_alloc
! =========================================================================== !
!> Free memory
subroutine vmtxel_free(this)
  class(vmtxel_t), intent(inout) :: this
 
  if(allocated(this%kpt))then;deallocate(this%kpt);endif
  if(allocated(this%pol))then;deallocate(this%pol);endif
  if(allocated(this%s1k))then;deallocate(this%s1k);endif
  if(allocated(this%s1))then;deallocate(this%s1);endif
  if(allocated(this%rnm))then;deallocate(this%rnm);endif
 
end subroutine vmtxel_free
! =========================================================================== !
!> Share matrix elements among all PEs
subroutine reduce_vmtxel(this)
  class(vmtxel_t), intent(inout) :: this
  real(DP), allocatable :: dummy(:,:)
 
  allocate(dummy (this%nmat, this%npol))
  dummy = this%s1
  if(allocated(dummy))then;deallocate(dummy);endif
 
end subroutine reduce_vmtxel
! =========================================================================== !
!> Transform the dipole operator from flat BSE indices to band indices
subroutine vmtxel_flat_to_band(this)
  class(vmtxel_t), intent(inout) :: this
  type(xctinfo) :: xct_
  integer :: ipol, is, ik, ic, iv
 
  if (.not. allocated(this%rnm)) then
    allocate(this%rnm (this%nband, this%mband, this%nk, this%ns, this%npol))
  end if
  ! GKA: FIXME Do this more elegantly
  xct_%nspin = this%ns
  xct_%ncb_fi = this%nband
  xct_%nvb_fi = this%mband
  do is=1,this%ns
    do ik=1,this%nk
      do iv=1,this%mband
        do ic=1,this%nband
          this%rnm(ic,iv,ik,is,:) = this%s1(bse_index(ik, ic, iv, is, xct_),:)
        enddo
      enddo
    enddo
  enddo
 
end subroutine vmtxel_flat_to_band
!> Transform the dipole operator from band indices to flat BSE indices
subroutine vmtxel_band_to_flat(this)
  class(vmtxel_t), intent(inout) :: this
  type(xctinfo) :: xct_
  integer :: ipol, is, ik, ic, iv
 
  if (.not. allocated(this%s1)) then
    allocate(this%s1 (this%nmat, this%npol))
  end if
  ! GKA: FIXME Do this more elegantly
  xct_%nspin = this%ns
  xct_%ncb_fi = this%nband
  xct_%nvb_fi = this%mband
  do is=1,this%ns
    do ik=1,this%nk
      do iv=1,this%mband
        do ic=1,this%nband
          this%s1(bse_index(ik, ic, iv, is, xct_),:) = this%rnm(ic,iv,ik,is,:)
        enddo
      enddo
    enddo
  enddo
 
end subroutine vmtxel_band_to_flat
! =========================================================================== !
!> Compute the dipole operator
subroutine compute_ik_vmtxel(this, ik, wfnc_fi, wfnvq_fi, gvec, qshift, &
                             crys, eqp)
  class(vmtxel_t), intent(inout) :: this
  integer, intent(in) :: ik
  type (gspace), intent(in) :: gvec
  type (wavefunction), intent(in) :: wfnc_fi
  type (wavefunction), intent(in) :: wfnvq_fi
  real(DP), intent(in), optional :: qshift
  type (crystal), intent(in), optional :: crys
  type (eqpinfo), intent(in), optional :: eqp
  type(xctinfo) :: xct_
  integer :: ipol, is, ic, iv
  real(DP) :: de
 
  if (this%npol==1) then
    if (this%opr.eq.0) then
      call mtxel_v(wfnc_fi,wfnvq_fi,gvec,qshift,this%nband,this%mband, &
                   this%s1k(:,:,:,1))
    elseif (this%opr.eq.1) then
      call mtxel_m(crys,wfnc_fi,wfnvq_fi,gvec,eqp,this%pol(:,1), &
                   this%nband,this%mband,this%s1k(:,:,:,1),ik,.true.)
    endif
  else
    do ipol=1,3
      this%pol(:,ipol) = 0.0d0
      this%pol(ipol,ipol) = 1.0d0
      if (this%opr.eq.0) then
        call mtxel_v(wfnc_fi,wfnvq_fi,gvec,qshift,this%nband,this%mband, &
                     this%s1k(:,:,:,ipol))
      else
        call mtxel_m(crys,wfnc_fi,wfnvq_fi,gvec,eqp,this%pol(:,ipol), &
                     this%nband,this%mband,this%s1k(:,:,:,ipol),ik,.true.)
      endif
    enddo
  endif
  ! GKA: FIXME Do this more elegantly
  xct_%nspin = this%ns
  xct_%ncb_fi = this%nband
  xct_%nvb_fi = this%mband
  do is=1,this%ns
    do ic=1,this%nband
      do iv=1,this%mband
        this%s1(bse_index(ik, ic, iv, is, xct_),:) = this%s1k(ic,iv,is,:)
      enddo
    enddo
  enddo
  ! Also compute the velocity operator, along with the dipole operator
  if (this%has_velocity .and. present(eqp)) then
    do is=1,this%ns
      do ic=1,this%nband
        do iv=1,this%mband
          de = (eqp%eclda(ic,ik,is) - eqp%evlda(iv,ik,is))
          this%vnm(ic,iv,ik,is,:) = this%s1k(ic,iv,is,:) * de
        enddo
      enddo
    enddo
  end if
 
end subroutine compute_ik_vmtxel
! =========================================================================== !
! Writing routines
! =========================================================================== !
!> Write the file, using the format determined with use_hdf5
subroutine write_vmtxel(this)
  class(vmtxel_t), intent(inout) :: this
 
    call this%write_vmtxel_bin()
 
end subroutine write_vmtxel
! =========================================================================== !
!> Write binary file
subroutine write_vmtxel_bin(this)
  class(vmtxel_t), intent(inout) :: this
  integer :: ipol
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
 
  if (this%is_master) then
    write(6,'(1x,a)') 'Writing matrix elements into vmtxel'
    do ipol=1,this%npol
      if (this%npol==1) then
        fname = 'vmtxel'
      else
        fname = 'vmtxel_'//suffix(ipol)
      endif
      call open_file(16, file=trim(fname), form='unformatted', status='replace')
      write(16) this%nk, this%nband, this%mband, this%ns, this%opr
      write(16) this%s1(:,ipol)
      call close_file(16)
    enddo
  endif
 
end subroutine write_vmtxel_bin
! =========================================================================== !
! =========================================================================== !
! Reading routines
! =========================================================================== !
!> Read the file, using the format determined with use_hdf5
subroutine read_vmtxel(this)
  class(vmtxel_t), intent(inout) :: this
 
    call this%read_vmtxel_bin()
 
end subroutine read_vmtxel
! =========================================================================== !
!> Read binary file
subroutine read_vmtxel_bin(this)
  class(vmtxel_t), intent(inout) :: this
  integer :: ii,ipol
  integer :: ic,iv,ik,is
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
 
  if (this%is_master) then
    write(6,'(1x,a)') 'Reading matrix elements from vmtxel'
    do ipol=1,this%npol
      if (this%npol==1) then
        fname = 'vmtxel'
      else
        fname = 'vmtxel_'//suffix(ipol)
      endif
      call open_file(16, file=trim(fname), form='unformatted', status='old')
      read(16) ik,ic,iv,is,ii
      if (ik.ne.this%nk.or.ic.ne.this%nband.or.iv.ne.this%mband &
        .or.is.ne.this%ns.or.ii.ne.this%opr) then
        write(0,'(a,5i6)') 'read  : ', ik,ic,iv,is,ii
        write(0,'(a,5i6)') 'needed: ', this%nk,this%nband,this%mband,this%ns,this%opr
        call die('parameter mismatch in vmtxel')
      endif
      read(16) this%s1(:,ipol)
      call close_file(16)
    enddo
  endif
 
end subroutine read_vmtxel_bin
! =========================================================================== !
end module vmtxel_m
