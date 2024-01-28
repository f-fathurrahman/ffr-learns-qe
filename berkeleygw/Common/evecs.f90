!==============================================================================
!
! Module evecs_m
!
! Created by GKA (2018).
!
! Objects holding the exciton wavefunctions, with corresponding IO routines.
!
!==============================================================================
module evecs_m
use global_m
use message_m
implicit none
public
! =========================================================================== !
!> BSE eigenvectors and eigenvalues
type evecs_t
  ! Dimensions
  integer :: ns !< Number of spins
  integer :: nv !< Number of valence bands
  integer :: nc !< Number of conduction bands
  integer :: nk !< Number of k-points
  integer :: nmat !< Size of the BSE hamiltonian (maximum number of eigenvalues)
                      !! We expect nmat=ns*nk*nc*nv
  integer :: usize !< Size of the eigenvectors (nmat or 2*nmat)
  integer :: neig !< Total number of eigenvectors
  integer :: meig=1 !< Number of eigenvectors held locally (array dimension)
  integer :: nq=1 !< Number of q-points (hard-coded)
  integer :: krnl=1 !< Kernel type
                      !! krnl = 0 --> spin triplet kernel, direct kernel only (only allow for nspin = 1)
                      !! 1 --> spin singlet kernel (default)
                      !! 2 --> local-fields + RPA, exchange kernel only
                      !! 3 --> spinor kernel
  logical :: tda = .true. !< Whether Tamm-Dancoff approximation is used
  logical :: has_Avc = .false. !< Whether the object holds eigenvectors as Avc
  logical :: has_flat_Avc = .false. !< Whether the object holds the flat eigenvectors
                                     !! Note: eigenvectors are named u when they are stored
                                     !! as flat arrays, and Avc when they are stored as (s,v,c,k).
  type(mf_header_t) :: mf_header !< The mean field header
  logical :: has_mf_header = .false. !< Whether the mean field header is initialized
  !> q-points, or center-of-mass momentum of exciton (qpt=0 for optical excitons).
  !! qpts(3,nq)
  real(DP), allocatable :: qpts(:,:)
  !> k-points
  !! kpts(3,nk)
  real(DP), allocatable :: kpts(:,:)
  !> BSE eigenvalues. This array is not distributed. Each node has all eigenvalues
  !! evals(neig)
  real(DP), allocatable :: evals(:)
  !> Right eigenvectors (flat)
  !! u_r(usize,meig)
  real(DP), allocatable :: u_r(:,:)
  !> Left eigenvectors (flat)
  !! u_l(usize,meig)
  real(DP), allocatable :: u_l(:,:)
  ! GKA: At the moment, the number of q-points (nq) is hard-coded to 1, and
  ! for simplicity, the e-h coefficient arrays do not have an nq dimension.
  ! In the hdf5 file, however, they do have this extra dimension.
  !> Right e-h coefficients
  !! Avc(ns,nv,nc,nk,meig)
  real(DP), allocatable :: Avc(:,:,:,:,:)
  !> Left e-h coefficients
  !! Avc_l(ns,nv,nc,nk,meig)
  real(DP), allocatable :: Avc_l(:,:,:,:,:)
  !> Right e-h coefficients, deexcitation component
  !! Bvc_r(ns,nv,nc,nk,meig)
  real(DP), allocatable :: Bvc_r(:,:,:,:,:)
  !> Left e-h coefficients, deexcitation component
  !! Bvc_l(ns,nv,nc,nk,meig)
  real(DP), allocatable :: Bvc_l(:,:,:,:,:)
  ! Parallel info
  logical :: is_master = .true. !< Whether we are on the master node
  logical :: is_distributed = .false. !< Whether each node owns a different set of eigenvalues
  integer :: npes=1 !< Number of processors
  integer :: inode=0 !< Rank of this node
  integer :: comm=0 !< MPI communicator
  integer :: mpi_info=0 !< MPI info
  integer :: global_meig !< The maximum value of meig accross processors
  integer, allocatable :: global_ieigs(:) !< The global eigenvalue indices of those owned locally
  integer, allocatable :: who_owns(:) !< The index of the node holding all eigenvector
  integer, allocatable :: local_ieigs(:) !< The local index on the node of all eigenvector
                                          !! This is a global array.
  ! IO parameters
  integer :: file_version=1 !< HDF5 file version
  logical :: is_open = .false. !< Whether the file is open for sequential reading
  integer :: nwrite=0 !< Number of eigenvectors to be written
  integer :: ieig_read=0 !< Reading counter
  integer :: uur = 11 !< Unit number for right eigenvectors
  integer :: uul = 12 !< Unit number for left eigenvectors
  character(30) :: default_fname_bin='eigenvectors'
  character(30) :: default_fname_hdf5='eigenvectors.h5'
  logical :: use_hdf5 = .false.
  ! ------------------------------------------------------------------------- !
  contains
  ! ------------------------------------------------------------------------- !
  ! Initialization
  procedure :: init_from_xctinfo => evecs_init_from_xctinfo
  procedure :: init => evecs_init
  ! Allocation / deallocation
  procedure :: alloc => evecs_alloc
  procedure :: alloc_global_arrays => evecs_alloc_global_arrays
  procedure :: alloc_local_arrays => evecs_alloc_local_arrays
  procedure :: free => evecs_free
  ! I/O interface
  procedure :: write_file => evecs_write
  procedure :: read_file => evecs_read
  procedure :: open_read => evecs_open_read
  procedure :: close_file => evecs_close
  procedure :: read_header => evecs_read_header
  procedure :: read_next_eigenvector => evecs_read_next_eigenvector
  ! Output
  procedure :: print_out_header_info => evecs_print_out_header_info
  procedure :: print_out_dimensions => evecs_print_out_dimensions
  ! Arrays manipulation
  procedure :: reshape_Avc => evecs_reshape_Avc
  ! Parallel setup
  procedure :: set_mpi_comm => evecs_set_mpi_comm
  procedure :: setup_paral_default => evecs_setup_paral_default
  procedure :: setup_paral_distribution => evecs_setup_paral_distribution
  procedure :: get_BSE_paral => evecs_get_BSE_paral
  ! MPI communication
  procedure :: broadcast_header => evecs_broadcast_header
  procedure :: broadcast_global_arrays => evecs_broadcast_global_arrays
  procedure :: broadcast_local_arrays => evecs_broadcast_local_arrays
  ! Testing
  procedure :: copy => evecs_copy
  procedure :: compare => evecs_compare
  ! Binary I/O
  procedure :: evecs_write_bin
  procedure :: evecs_read_bin
  procedure :: evecs_read_and_broadcast_bin
  procedure :: evecs_open_read_bin
  procedure :: evecs_close_bin
  procedure :: evecs_read_header_bin
  procedure :: evecs_read_dimensions_bin
  procedure :: evecs_count_neig_and_rewind_bin
  procedure :: evecs_read_single_eigenvector_bin
  procedure :: evecs_skip_eigenvectors_bin
  procedure :: evecs_read_and_distribute_bin
  ! HDF5 I/O
end type evecs_t
! =========================================================================== !
contains
! =========================================================================== !
!> Initialize manually, mostly for testing purposes.
subroutine evecs_init(this, ns, nk, nv, nc, neig, meig, tda, mf_header)
  class(evecs_t), intent(inout) :: this
  integer, intent(in) :: ns !< Number of spin
  integer, intent(in) :: nk !< Number of k-points
  integer, intent(in) :: nv !< Number of valence bands
  integer, intent(in) :: nc !< Number of conduction bands
  integer, intent(in), optional :: neig !< Number of eigenvalues
  integer, intent(in), optional :: meig !< Number of eigenvectors for array dimension
                                        !! Cannot be larger than neig
  logical, intent(in), optional :: tda !< Use Tamm-Dancoff approximation
  type (mf_header_t), intent(in), optional :: mf_header !< The mean field header
 
  this%ns = ns
  this%nk = nk
  this%nv = nv
  this%nc = nc
  if (present(neig)) then
    this%neig = neig
  else
    this%neig = this%ns * this%nk * this%nv * this%nc
  end if
  if (present(meig)) then
    this%meig = meig
  else
    this%meig = this%neig
  end if
  this%tda = .true.
  if (present(tda)) this%tda = tda
  this%nmat = this%ns * this%nk * this%nv * this%nc
  if (this%tda) then
    this%usize = this%nmat
  else
    this%usize = 2 * this%nmat
  end if
  allocate(this%kpts (3, this%nk))
  this%kpts = 0.0D0
  allocate(this%qpts (3, this%nq))
  this%qpts = 0.0D0
  this%ieig_read = 0
  this%nwrite = 0
  if (present(mf_header)) then
    this%mf_header = mf_header
    this%has_mf_header = .true.
  end if
  ! Get parallel info
  call this%setup_paral_default()
 
end subroutine evecs_init
! =========================================================================== !
!> Initialize the eigenvectors using the xctinfo object typically found
!! in a BSE calculation.
subroutine evecs_init_from_xctinfo(this, xct, kg, nmat, neig, meig, mf_header)
  class(evecs_t), intent(inout) :: this
  type (xctinfo), intent(in) :: xct !< Information on the BSE calculation
  type (grid), intent(in) :: kg !< Information on the k-point grid
  integer, intent(in) :: nmat
  integer, intent(in) :: neig
  integer, intent(in) :: meig
  type (mf_header_t), intent(in), optional :: mf_header !< The mean field header
  integer :: ii, ik
 
  this%ieig_read = 0
  this%nmat = nmat
  this%neig = neig
  this%meig = meig
  this%nwrite = neig
  this%ns = xct%nspin
  this%nv = xct%nvb_fi
  this%nc = xct%ncb_fi
  this%nk = xct%nkpt_fi
  this%tda = xct%tda
  this%krnl = xct%krnl
  this%use_hdf5 = (xct%use_hdf5 .and. xct%use_hdf5_output)
  allocate(this%kpts (3, this%nk))
  this%kpts = kg%f
  allocate(this%qpts (3, this%nq))
  this%qpts(:,:) = 0.0D0
  if (xct%qflag == 1) then
    continue
  else if (xct%qflag == 2 .or. xct%qflag == 0) then
    this%qpts(:,1) = xct%finiteq
  end if
  if (present(mf_header)) then
    this%mf_header = mf_header
    this%has_mf_header = .true.
  end if
  ! Get parallel info
  call this%get_BSE_paral()
 
end subroutine evecs_init_from_xctinfo
! =========================================================================== !
!> Copy object into a new instance.
!! For testing purposes.
subroutine evecs_copy(this, other)
  class(evecs_t), intent(inout) :: this
  class(evecs_t), intent(out) :: other
 
  other%ns = this%ns
  other%nk = this%nk
  other%nc = this%nc
  other%nv = this%nv
  other%neig = this%neig
  other%meig = this%meig
  other%nmat = this%nmat
  other%usize = this%usize
  other%tda = this%tda
  other%has_Avc = this%has_Avc
  other%has_flat_Avc = this%has_flat_Avc
  other%mf_header = this%mf_header
  other%has_mf_header = this%has_mf_header
  other%is_master = this%is_master
  other%is_distributed = this%is_distributed
  other%global_meig = this%global_meig
  other%is_open = this%is_open
  other%nwrite = this%nwrite
  other%ieig_read = this%ieig_read
  other%uur = this%uur
  other%uul = this%uul
  other%default_fname_bin = this%default_fname_bin
  other%default_fname_hdf5 = this%default_fname_hdf5
  other%use_hdf5 = this%use_hdf5
  call other%alloc(with_Avc=.true.)
  other%qpts = this%qpts
  if (allocated(this%kpts) .and. allocated(other%kpts)) other%kpts = this%kpts
  if (allocated(this%evals) .and. allocated(other%evals)) other%evals = this%evals
  if (allocated(this%u_r) .and. allocated(other%u_r)) other%u_r = this%u_r
  if (allocated(this%u_l) .and. allocated(other%u_l)) other%u_l = this%u_l
  if (allocated(this%Avc) .and. allocated(other%Avc)) other%Avc = this%Avc
  if (allocated(this%Avc_l) .and. allocated(other%Avc_l)) other%Avc_l = this%Avc_l
  if (allocated(this%Bvc_r) .and. allocated(other%Bvc_r)) other%Bvc_r = this%Bvc_r
  if (allocated(this%Bvc_l) .and. allocated(other%Bvc_l)) other%Bvc_l = this%Bvc_l
  if (allocated(this%global_ieigs) .and. allocated(other%global_ieigs)) other%global_ieigs = this%global_ieigs
  if (allocated(this%who_owns) .and. allocated(other%who_owns)) other%who_owns = this%who_owns
  if (allocated(this%local_ieigs) .and. allocated(other%local_ieigs)) other%local_ieigs = this%local_ieigs
 
end subroutine evecs_copy
! =========================================================================== !
!> Compare a evecs object against another, and report any discrepancy.
!! For testing purposes.
subroutine evecs_compare(this, other, ndiff, verbose)
  class(evecs_t), intent(in) :: this
  class(evecs_t), intent(in) :: other
  integer, intent(out) :: ndiff
  logical, intent(in), optional :: verbose
  logical :: verbose_ = .True.
  real(DP), parameter :: tol=1.0d-5
  real(DP) :: checksum
  integer :: unt=6
  if (present(verbose)) verbose_ = verbose
 
  ndiff = 0
  checksum = abs(sum(this%qpts) - sum(other%qpts))
  if (checksum .gt. tol) then
    if (verbose_) write(unt,*) 'Error in checksum of: qpt. ','delta=',checksum
    ndiff = ndiff + 1
  end if
  if (allocated(this%kpts) .and. allocated(other%kpts)) then
    checksum = abs(sum(this%kpts) - sum(other%kpts))
    if (checksum .gt. tol) then
      if (verbose_) write(unt,*) 'Error in checksum of: kpts. ','delta=',checksum
      ndiff = ndiff + 1
    end if
  end if
  if (allocated(this%Avc) .and. allocated(other%Avc)) then
    checksum = abs(sum(this%Avc) - sum(other%Avc))
    if (checksum .gt. tol) then
      if (verbose_) write(unt,*) 'Error in checksum of: Avc. ','delta=',checksum
      ndiff = ndiff + 1
    end if
  end if
  if (allocated(this%Avc_l) .and. allocated(other%Avc_l)) then
    checksum = abs(sum(this%Avc_l) - sum(other%Avc_l))
    if (checksum .gt. tol) then
      if (verbose_) write(unt,*) 'Error in checksum of: Avc_l. ','delta=',checksum
      ndiff = ndiff + 1
    end if
  end if
  if (allocated(this%Bvc_r) .and. allocated(other%Bvc_r)) then
    checksum = abs(sum(this%Bvc_r) - sum(other%Bvc_r))
    if (checksum .gt. tol) then
      if (verbose_) write(unt,*) 'Error in checksum of: Bvc_r. ','delta=',checksum
      ndiff = ndiff + 1
    end if
  end if
  if (allocated(this%Bvc_l) .and. allocated(other%Bvc_l)) then
    checksum = abs(sum(this%Bvc_l) - sum(other%Bvc_l))
    if (checksum .gt. tol) then
      if (verbose_) write(unt,*) 'Error in checksum of: Bvc_l. ','delta=',checksum
      ndiff = ndiff + 1
    end if
  end if
 
end subroutine evecs_compare
! =========================================================================== !
!> Pass the MPI communicator
subroutine evecs_set_mpi_comm(this, comm)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: comm !< MPI communicator
 
  this%is_master = (this%inode .eq. 0)
 
end subroutine evecs_set_mpi_comm
! =========================================================================== !
!> Setup the default parallel scheme, where each node owns a copy
!! of all eigenvectors. No information on the dimensions is needed.
subroutine evecs_setup_paral_default(this, comm, neig)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: comm !< MPI communicator
  integer, intent(in), optional :: neig !< Number of eigenvalues
  integer :: ieig
 
  this%is_distributed = .false.
  call this%set_mpi_comm(comm)
  this%global_meig = this%meig
  if (present(neig)) this%neig = neig
  ! No distribution: Master owns everything (but evecs might be broadcast)
  if (this%neig .gt. 0) then
    if (allocated(this%global_ieigs)) then
      if(allocated(this%global_ieigs))then;deallocate(this%global_ieigs);endif
    end if
    if (allocated(this%local_ieigs)) then
      if(allocated(this%local_ieigs))then;deallocate(this%local_ieigs);endif
    end if
    if (allocated(this%who_owns)) then
      if(allocated(this%who_owns))then;deallocate(this%who_owns);endif
    end if
    allocate(this%global_ieigs (this%neig))
    allocate(this%local_ieigs (this%neig))
    allocate(this%who_owns (this%neig))
    this%who_owns = 0 !< The index of the node holding each eigenvector
    this%local_ieigs = 0 !< The local index on the node of all eigenvector
    this%global_ieigs = 0 !< The global index of each local eigenvector
    do ieig=1,this%neig
      this%local_ieigs(ieig) = ieig
      this%global_ieigs(ieig) = ieig
    enddo
  end if
 
end subroutine evecs_setup_paral_default
! =========================================================================== !
!> Initialize parallel info by relying on peinfo, which has been set up for BSE
!! Assume that neig is known
subroutine evecs_get_BSE_paral(this)
  class(evecs_t), intent(inout) :: this
  integer :: ieig_local, ieig_global, peadd
 
  call this%setup_paral_default()
  if (allocated(peinf%neig)) then
    this%meig = peinf%neig(peinf%inode+1)
  end if
  this%global_meig = this%meig
  if (.not. allocated(this%global_ieigs)) then
    allocate(this%global_ieigs (this%global_meig))
  end if
  if (.not. allocated(this%local_ieigs)) then
    allocate(this%local_ieigs (this%neig))
  end if
  if (.not. allocated(this%who_owns)) then
    allocate(this%who_owns (this%neig))
  end if
  this%who_owns = 0 !< The index of the node holding each eigenvector
  this%local_ieigs = 0 !< The local index on the node of all eigenvector
  this%global_ieigs = 0 !< The global index of each local eigenvector
  if (allocated(peinf%peig)) then
    this%is_distributed = .true.
    do ieig_local=1,this%meig
      this%global_ieigs(ieig_local) = peinf%peig(peinf%inode+1, ieig_local)
    end do
    do peadd=1,peinf%npes
      do ieig_local=1,peinf%neig(peadd)
        ieig_global = peinf%peig(peadd, ieig_local)
        this%who_owns(ieig_global) = peadd - 1
        this%local_ieigs(ieig_global) = ieig_local
      end do
    end do
  else
    this%is_distributed = .false.
    do ieig_local=1,this%meig
      this%global_ieigs(ieig_local) = ieig_local
      this%local_ieigs(ieig_local) = ieig_local
    end do
  end if
 
end subroutine evecs_get_BSE_paral
! =========================================================================== !
!> Setup the distribution of eigenvectors among all nodes.
!! The number of eigenvalues must be known prior to calling this function
!! or be given explicitly.
!! If there are N nodes in the group numbered from 0 to N-1
!! and there are M eigenvalues numbered from 1 to M,
!! then eigenvalue i will be owned by node mod(i-1, N)
subroutine evecs_setup_paral_distribution(this, neig, comm)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: neig !< Number of eigenvectors
  integer, intent(in), optional :: comm !< MPI communicator
  integer :: ieig, ieig_local
  integer, allocatable :: tmp_local_ieigs(:)
 
  call this%set_mpi_comm(comm=comm)
  if (present(neig)) then
    this%neig = neig
  else if (this%neig .eq. 0) then
    call die('evecs_setup_paral_distribution was called without specifying neig, but neig is unknown.')
  end if
  if (allocated(this%global_ieigs)) then
    if(allocated(this%global_ieigs))then;deallocate(this%global_ieigs);endif
  end if
  if (allocated(this%local_ieigs)) then
    if(allocated(this%local_ieigs))then;deallocate(this%local_ieigs);endif
  end if
  if (allocated(this%who_owns)) then
    if(allocated(this%who_owns))then;deallocate(this%who_owns);endif
  end if
  ! This array only needs to be of length meig
  ! but we are using neig as an upper bound.
  allocate(this%global_ieigs (this%neig))
  allocate(this%local_ieigs (this%neig))
  allocate(this%who_owns (this%neig))
  this%who_owns = 0 !< The index of the node holding each eigenvector
  this%local_ieigs = 0 !< The local index on the node of all eigenvector
  this%global_ieigs = 0 !< The global index of each local eigenvector
  ! Manual distribution of the eigenvectors
  ieig_local = 0
  do ieig=1,this%neig
    this%who_owns(ieig) = mod(ieig-1, this%npes)
    if (this%who_owns(ieig) .eq. this%inode) then
      ieig_local = ieig_local + 1
      this%local_ieigs(ieig) = ieig_local
      this%global_ieigs(ieig_local) = ieig
    endif
  enddo
  this%meig = ieig_local
  ! Figure out maximum number of eigenvectors locally stored
  this%global_meig = this%meig
  ! Communicate all local indices
  this%is_distributed = .true.
 
end subroutine evecs_setup_paral_distribution
! =========================================================================== !
!> Allocate memory
subroutine evecs_alloc(this, with_Avc, with_flat_Avc)
  class(evecs_t), intent(inout) :: this
  logical, intent(in), optional :: with_Avc !< Allocate e-h coefficients
  logical, intent(in), optional :: with_flat_Avc !< Allocate u arrays
 
  call this%alloc_global_arrays()
  call this%alloc_local_arrays(with_Avc=with_Avc, with_flat_Avc=with_flat_Avc)
 
end subroutine evecs_alloc
!> Allocate memory
subroutine evecs_alloc_global_arrays(this)
  class(evecs_t), intent(inout) :: this
 
  ! Global arrays
  ! -------------
  ! k-points
  if (.not. allocated(this%kpts)) then
    allocate(this%kpts (3, this%nk))
  end if
  ! Eigenvalues
  if (.not. allocated(this%evals)) then
    allocate(this%evals (this%neig))
  end if
  if (.not. allocated(this%global_ieigs)) then
    allocate(this%global_ieigs (this%meig))
  end if
  ! q-points
  if (.not. allocated(this%qpts)) then
    allocate(this%qpts (3, this%nq))
    this%qpts = 0.0D0
  end if
 
end subroutine evecs_alloc_global_arrays
! =========================================================================== !
!> Allocate memory
subroutine evecs_alloc_local_arrays(this, with_Avc, with_flat_Avc)
  class(evecs_t), intent(inout) :: this
  logical, intent(in), optional :: with_Avc !< Whether Avc arrays should be allocated
  logical, intent(in), optional :: with_flat_Avc !< Whether u arrays should be allocated
  logical :: with_Avc_, with_flat_Avc_ !< Whether u arrays should be allocated
 
  with_Avc_ = .false.
  with_flat_Avc_ = .true.
  if (present(with_Avc)) then
    with_Avc_ = with_Avc
  end if
  if (present(with_flat_Avc)) then
    with_flat_Avc_ = with_flat_Avc
  end if
  ! Flat Avc
  if ( this%tda ) then
    this%usize = this%nmat
    if (with_flat_Avc_) then
      allocate(this%u_r (this%usize,this%meig))
    end if
  else
    this%usize = 2*this%nmat
    if (with_flat_Avc_) then
      allocate(this%u_r (this%usize, this%meig))
      allocate(this%u_l (this%usize, this%meig))
    end if
  end if
  ! Avc
  if (with_Avc_) then
    if ( this%tda ) then
      allocate(this%Avc (this%ns,this%nv,this%nc,this%nk,this%meig))
    else
      allocate(this%Avc (this%ns,this%nv,this%nc,this%nk, this%meig))
      allocate(this%Avc_l (this%ns,this%nv,this%nc,this%nk, this%meig))
      allocate(this%Bvc_r (this%ns,this%nv,this%nc,this%nk, this%meig))
      allocate(this%Bvc_l (this%ns,this%nv,this%nc,this%nk, this%meig))
    end if
  end if
 
end subroutine evecs_alloc_local_arrays
! =========================================================================== !
!> Free memory
subroutine evecs_free(this)
  class(evecs_t), intent(inout) :: this
 
  if(allocated(this%kpts))then;deallocate(this%kpts);endif
  if(allocated(this%qpts))then;deallocate(this%qpts);endif
  if(allocated(this%evals))then;deallocate(this%evals);endif
  if(allocated(this%u_r))then;deallocate(this%u_r);endif
  if(allocated(this%u_l))then;deallocate(this%u_l);endif
  if(allocated(this%Avc))then;deallocate(this%Avc);endif
  if(allocated(this%Avc_l))then;deallocate(this%Avc_l);endif
  if(allocated(this%Bvc_r))then;deallocate(this%Bvc_r);endif
  if(allocated(this%Bvc_l))then;deallocate(this%Bvc_l);endif
  if(allocated(this%global_ieigs))then;deallocate(this%global_ieigs);endif
  if(allocated(this%local_ieigs))then;deallocate(this%local_ieigs);endif
  if(allocated(this%who_owns))then;deallocate(this%who_owns);endif
  ! Reinitialize some variables
  this%has_Avc = .false.
  this%has_flat_Avc = .false.
  this%is_distributed = .false.
  this%ieig_read = 0
  this%nwrite = 0
 
end subroutine evecs_free
! =========================================================================== !
!> Transform the flat eigenvectors into the Avc components
subroutine evecs_reshape_Avc(this, force)
  class(evecs_t), intent(inout) :: this
  logical, intent(in), optional :: force !< Do it even if done before
  integer :: ieig,iflat,ik,ic,iv,is
  logical :: force_ !< Do it even if done before
 
  if (this%has_Avc) then
    force_ = .false.
    if (present(force)) then
      force_ = force
    end if
    if (.not. force_) return
  end if
  ! local worker may not even have the data
  if (.not. allocated(this%u_r)) then
    this%has_Avc = .true.
    return
  end if
  if (.not. allocated(this%Avc)) then
    call this%alloc_local_arrays(with_Avc=.true., with_flat_Avc=.false.)
  end if
  if ( this%tda ) then
    do ieig=1,this%meig
      iflat = 0
      do ik=1,this%nk
        do ic=1,this%nc
          do iv=1,this%nv
            do is=1,this%ns
              iflat = iflat + 1
              this%Avc(is,iv,ic,ik,ieig) = this%u_r(iflat,ieig)
            enddo
          enddo
        enddo
      enddo
    enddo
  else
    do ieig=1,this%meig
      iflat = 0
      do ik=1,this%nk
        do ic=1,this%nc
          do iv=1,this%nv
            do is=1,this%ns
              iflat = iflat + 1
              this%Avc(is,iv,ic,ik,ieig) = this%u_r(iflat,ieig)
              this%Bvc_r(is,iv,ic,ik,ieig) = this%u_r(this%nmat+iflat,ieig)
              this%Avc_l(is,iv,ic,ik,ieig) = this%u_l(iflat,ieig)
              this%Bvc_l(is,iv,ic,ik,ieig) = this%u_l(this%nmat+iflat,ieig)
            enddo
          enddo
        enddo
      enddo
    enddo
  end if !tda
  this%has_Avc = .true.
 
end subroutine evecs_reshape_Avc
! =========================================================================== !
!> Broadcast dimensions and k-points
subroutine evecs_broadcast_header(this, comm)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: comm !< MPI communicator
  integer :: dims(11)
  integer :: tda
 
  call this%set_mpi_comm(comm)
  if (this%is_master) then
    if (this%tda) then
      tda = 1
    else
      tda = 0
    end if
    dims(:) = [this%nmat, this%usize, this%neig, this%ns, this%nk, this%nc, this%nv, this%meig, this%nq, this%krnl, tda]
  end if
  ! Broadcast
  if (.not. this%is_master) then
    this%nmat = dims(1)
    this%usize = dims(2)
    this%neig = dims(3)
    this%ns = dims(4)
    this%nk = dims(5)
    this%nc = dims(6)
    this%nv = dims(7)
    this%meig = dims(8)
    this%nq = dims(9)
    this%krnl = dims(10)
    tda = dims(11)
    if (tda.eq.0) then
      this%tda = .false.
    else
      this%tda = .true.
    end if
  end if
 
end subroutine evecs_broadcast_header
! =========================================================================== !
!> Broadcast the arrays that should be owned by all processors
!! Assume that the arrays are already allocated
subroutine evecs_broadcast_global_arrays(this)
  class(evecs_t), intent(inout) :: this
 
 
end subroutine evecs_broadcast_global_arrays
! =========================================================================== !
!> Broadcast the eigenvectors
subroutine evecs_broadcast_local_arrays(this)
  class(evecs_t), intent(inout) :: this
  integer :: dims
 
  this%is_distributed = .false.
 
end subroutine evecs_broadcast_local_arrays
! =========================================================================== !
! Interface routines for binary / HDF5 selection
! =========================================================================== !
!> Main reading routine
subroutine evecs_read(this, fname, neig_read, ieig_offset, distribute, master_only, comm)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname
  integer, intent(in), optional :: neig_read !< Number of eigenvectors to be read
  integer, intent(in), optional :: ieig_offset !< Number of eigenvectors to skip befor reading
                                                !< The first eigenvector index is thus ieig_offset+1
  logical, intent(in), optional :: distribute !< Distribute eigenvectors among all processors
  logical, intent(in), optional :: master_only !< Only master reads eivenvectors
  integer, intent(in), optional :: comm !< MPI communicator
 
  if (this%use_hdf5) then
  else
    call this%evecs_read_bin(fname=fname, neig_read=neig_read, ieig_offset=ieig_offset, &
      with_Avc=.true., distribute=distribute, master_only=master_only, comm=comm)
  end if
 
end subroutine evecs_read
! =========================================================================== !
!> Open file for sequential reading
subroutine evecs_open_read(this, fname, file_id)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname !< File name
  integer, intent(out), optional :: file_id !< File ID for sequential read
 
  if (this%use_hdf5) then
  else
    call this%evecs_open_read_bin(fname)
  end if
 
end subroutine evecs_open_read
! =========================================================================== !
!> Interface function for reading header
subroutine evecs_read_header(this, fname, file_id, broadcast)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname !< File name
  integer, intent(inout), optional :: file_id !< File ID for HDF5 sequential read
  logical, intent(in), optional :: broadcast !< Broadcast after
 
  if (this%use_hdf5) then
  else
    call this%evecs_read_header_bin()
  end if
  if (present(broadcast)) then
    if (broadcast) then
      call this%broadcast_header()
      call this%broadcast_global_arrays()
    end if
  end if
 
end subroutine evecs_read_header
! =========================================================================== !
!> Close file after sequential reading
subroutine evecs_close(this, file_id)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: file_id
 
  if (this%use_hdf5) then
  else
    call this%evecs_close_bin()
  end if
 
end subroutine evecs_close
! =========================================================================== !
!> Read a single eigenvalue and eigenvector in serial
!! and store it in position 1 of the evals and evecs arrays
!!
!! If one wants needs this function to be compatible with both binary and hdf5,
!! then one should make sure the file is open beforehand and file_id is passed.
!! Else, for a call that only needs to be compatible with hdf5
!! the file_id can be omitted and the file will be closed after
!!
!! This is the sequential reading mode, which minimizes memory requirements.
!! It is used in summariza_eigenvectors.
subroutine evecs_read_next_eigenvector(this, fname, file_id)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname !< File name
  integer, intent(inout), optional :: file_id !< File ID for HDF5 sequential read
 
  if (this%use_hdf5) then
  else
    call this%evecs_read_single_eigenvector_bin()
  end if
 
end subroutine evecs_read_next_eigenvector
! =========================================================================== !
! Binary reading routines
! =========================================================================== !
!> Read eigenvectors in binary format
subroutine evecs_read_bin(this, fname, neig_read, ieig_offset, &
                          with_Avc, distribute, master_only, comm)
  class(evecs_t), intent(inout) :: this
  character(len=*), intent(in), optional :: fname
  integer, intent(in), optional :: neig_read !< Number of eigenvectors to be read
                                                !! It is optional only if HDF5 is used
  integer, intent(in), optional :: ieig_offset !< number of eigenvectors to skip befor reading
                                                !! the first eigenvector index is thus ieig_offset+1
                                                !! will not be used if eigenvectors are distributed
  logical, intent(in), optional :: with_Avc !< Reshape into Avc
  logical, intent(in), optional :: distribute !< spread the eigenvectors among processors
  logical, intent(in), optional :: master_only !< Only master reads and holds eigenvectors
  integer, intent(in), optional :: comm !< MPI communicator
  character(len=512) :: fname_
 
  ! Setup the distribution if requested
  if (present(distribute)) then
    if (distribute) then
      this%is_distributed = distribute
    endif
  endif
  if (this%is_distributed) then
    if (present(ieig_offset)) then
      call die('Cannot read and distribute binary file with ieig_offset')
    end if
    call this%evecs_read_and_distribute_bin(fname, neig_read=neig_read, &
                                            with_Avc=.true., comm=comm)
  else
    call this%evecs_read_and_broadcast_bin(fname, neig_read=neig_read, ieig_offset=ieig_offset, &
                    with_Avc=.true., master_only=master_only, comm=comm)
  end if
 
end subroutine evecs_read_bin
! =========================================================================== !
!> Open the binary eigenvectors file(s) for reading
subroutine evecs_open_read_bin(this, rootname)
  class(evecs_t), intent(inout) :: this
  character(*), intent(in), optional :: rootname !< File name, to be appended by _r or _l
                                                  !! for non-tda.
  character(60) :: fname
 
  if (present(rootname)) then
    fname = trim(rootname)
  else
    fname = this%default_fname_bin
  end if
  if (this%tda) then
    call open_file(unit=this%uur,file=trim(fname),form='unformatted',status='old')
  else
    call open_file(unit=this%uul,file=trim(trim(fname)//'_l'),form='unformatted',status='old')
    call open_file(unit=this%uur,file=trim(trim(fname)//'_r'),form='unformatted',status='old')
  end if
  this%is_open = .true.
 
end subroutine evecs_open_read_bin
! =========================================================================== !
!> Read the dimensions from binary file(s) then count number of eigenvectors
subroutine evecs_read_header_bin(this, neig)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: neig !< Provide the number of eigenvectors
                                         !! so it needs not be counted
 
  call this%evecs_read_dimensions_bin()
  if (present(neig)) then
    this%neig = neig
    this%meig = neig
  else
    call this%evecs_count_neig_and_rewind_bin()
  end if
  if (allocated(this%qpts)) then
    if(allocated(this%qpts))then;deallocate(this%qpts);endif
  end if
  allocate(this%qpts (3,this%nq))
  this%qpts = 0.0D0
 
end subroutine evecs_read_header_bin
! =========================================================================== !
!> Read the dimensions from binary file(s)
subroutine evecs_read_dimensions_bin(this)
  class(evecs_t), intent(inout) :: this
  real(DP), allocatable :: kpts_l(:,:)
  integer :: ns_, nv_, nc_, nk_, nmat_
  integer :: ik
 
  if (.not. this%is_open) then
    call this%evecs_open_read_bin()
  end if
  read(this%uur) this%ns
  read(this%uur) this%nv
  read(this%uur) this%nc
  read(this%uur) this%nk
  this%nmat = this%ns * this%nv * this%nc * this%nk
  this%neig = this%nmat !< Temporary value
  this%meig = this%neig
  if ( this%tda ) then
    this%usize = this%nmat
  else
    this%usize = 2*this%nmat
  end if
  if (.not. this%tda ) then
    read(this%uul) ns_
    read(this%uul) nv_
    read(this%uul) nc_
    read(this%uul) nk_
    nmat_ = ns_*nv_*nc_*nk_
    if( this%nmat .ne. nmat_) then
      call die('Error found: nmat should be equal for left and right eigenvalue matrices')
    end if
  end if
  if (allocated(this%kpts)) then
    if(allocated(this%kpts))then;deallocate(this%kpts);endif
  end if
  allocate(this%kpts (3,this%nk))
  read(this%uur) this%kpts(:,:)
  if( .not. this%tda ) then
    allocate(kpts_l (3,this%nk))
    read(this%uul) kpts_l(:,:)
! Check:
    do ik=1,this%nk
      if( abs(this%kpts(1,ik)-kpts_l(1,ik)) > 1.d06) then
        if( abs(this%kpts(2,ik)-kpts_l(2,ik)) > 1.d06) then
          if( abs(this%kpts(3,ik)-kpts_l(3,ik)) > 1.d06) then
            call die('Inconsistency in k-points found in left and right eigenvalues')
          end if
        end if
      end if
    end do
    if(allocated(kpts_l))then;deallocate(kpts_l);endif
  end if
 
end subroutine evecs_read_dimensions_bin
! =========================================================================== !
!> Read a single eigenvalue and eigenvector(s)
subroutine evecs_read_single_eigenvector_bin(this, ieig_evec, ieig_eval)
  class(evecs_t), intent(inout) :: this
  integer, intent(in), optional :: ieig_evec !< Index of the eigenvector for storage
  integer, intent(in), optional :: ieig_eval !< Index of the eigenvalue for storage
  integer :: ieig_evec_, ieig_eval_
  real(DP) :: eval_l
 
  ieig_eval_ = 1
  ieig_evec_ = 1
  if (present(ieig_eval)) ieig_eval_ = ieig_eval
  if (present(ieig_evec)) ieig_evec_ = ieig_evec
  if ( this%tda ) then
     read(this%uur) this%evals(ieig_eval_)
     read(this%uur) this%u_r(:,ieig_evec_)
  else
     read(this%uur) this%evals(ieig_eval_)
     read(this%uur) this%u_r(:,ieig_evec_)
     read(this%uul) eval_l
     read(this%uul) this%u_l(:,ieig_evec_)
     if( abs(eval_l - this%evals(ieig_eval_))>1.d-12) then
        call die('Inconsistency in energies in left and right eigenfuncion files')
     end if
  end if
  this%has_Avc = .false.
 
end subroutine evecs_read_single_eigenvector_bin
! =========================================================================== !
!> Skip eigenvectors in the file, but read and store the eigenvalues.
subroutine evecs_skip_eigenvectors_bin(this, nskip)
  class(evecs_t), intent(inout) :: this
  integer, intent(in) :: nskip !< Number of eigenvector to skip
  integer :: ii
  real(DP) :: eval_l
 
  if ( this%tda ) then
    do ii=1,nskip
      read(this%uur)
      read(this%uur)
    end do
  else
    do ii=1,nskip
      read(this%uur)
      read(this%uur)
      read(this%uul)
      read(this%uul)
    end do
  end if
 
end subroutine evecs_skip_eigenvectors_bin
! =========================================================================== !
!> Count number of eigenvectors in binary file then rewind
!! at the position after header.
subroutine evecs_count_neig_and_rewind_bin(this)
  class(evecs_t), intent(inout) :: this
  integer :: ii, neig
  integer :: iostatus
  real(dp) :: eval
 
  neig = 0
  if ( this%tda ) then
    do ii=1,this%nmat
      read(this%uur, iostat=iostatus) eval
      if (iostatus .eq. 0) then
        neig = neig + 1
      else if (iostatus .lt. 0) then
        exit
      else
        call die('Error counting eigenvectors')
      end if
      read(this%uur)
    end do
  else
    do ii=1,this%nmat
      read(this%uur, iostat=iostatus) eval
      if (iostatus .eq. 0) then
        neig = neig + 1
      else if (iostatus .lt. 0) then
        exit
      else
        call die('Error counting eigenvectors')
      end if
      read(this%uur)
      read(this%uul)
      read(this%uul)
    end do
  end if
  ! Now rewind the file
  rewind(this%uur)
  if (.not. this%tda) then
    rewind(this%uul)
  end if
  call this%evecs_read_dimensions_bin()
  this%neig = neig
  this%meig = neig
 
end subroutine evecs_count_neig_and_rewind_bin
! =========================================================================== !
!> Binary reading of all eigenvectors in serial and broadcast all eigenvectors
!! to every node.
subroutine evecs_read_and_broadcast_bin(this, rootname, neig_read, ieig_offset, with_Avc, master_only, comm)
  class(evecs_t), intent(inout) :: this
  character(*), intent(in), optional :: rootname
  integer, intent(in), optional :: neig_read !< Number of eigenvectors to be read
  integer, intent(in), optional :: ieig_offset !< Number of eigenvectors to skip befor reading
                                                !< The first eigenvector index is thus ieig_offset+1
  logical, intent(in), optional :: with_Avc
  logical, intent(in), optional :: master_only !< Only master reads eivenvectors
  integer, intent(in), optional :: comm !< MPI communicator
  logical :: master_only_
  logical :: with_Avc_
  integer :: ieig
  integer :: nread, nskip
 
  call this%set_mpi_comm(comm=comm)
  with_Avc_ = .true.
  if (present(with_Avc)) with_Avc_ = with_Avc
  if (present(master_only)) then
    master_only_ = master_only
  else
    master_only_ = .false.
  end if
  if (master_only_) this%is_distributed = .false.
  if (this%is_master) then
    call this%evecs_open_read_bin(rootname)
    call this%evecs_read_header_bin()
    nskip = 0
    if (present(ieig_offset)) nskip = ieig_offset
    if (present(neig_read)) then
      nread = neig_read
      this%neig = nread
      this%meig = nread
    else if (nskip > 0) then
      nread = this%neig - nskip
      this%neig = nread
      this%meig = nread
    else
      nread = this%neig
    end if
  end if
  call this%broadcast_header(comm=comm)
  call this%setup_paral_default(comm=comm)
  ! Allocate memory
  if (master_only_) then
    call this%alloc_global_arrays()
    if (this%is_master) then
      call this%alloc_local_arrays(with_Avc=with_Avc_)
    end if
  else
    call this%alloc(with_Avc=with_Avc_)
  end if
  if (this%is_master) then
    if (nskip > 0) then
      call this%evecs_skip_eigenvectors_bin(nskip)
    end if
    do ieig=1,nread
      call this%evecs_read_single_eigenvector_bin(ieig_evec=ieig, ieig_eval=ieig)
    end do
    call this%evecs_close_bin()
  end if
  this%has_flat_Avc = .true.
  if (with_Avc_) then
    call this%reshape_Avc()
  end if
 
end subroutine evecs_read_and_broadcast_bin
! =========================================================================== !
!> Binary reading of all eigenvectors in serial then distribute eigenvectors
!! among the nodes.
subroutine evecs_read_and_distribute_bin(this, rootname, neig_read, with_Avc, comm)
  class(evecs_t), intent(inout) :: this
  character(*), intent(in), optional :: rootname
  integer, intent(in), optional :: neig_read !< Number of eigenvectors to be read
  logical, intent(in), optional :: with_Avc !< Reshape into Avc
  integer, intent(in), optional :: comm !< MPI communicator
  logical :: with_Avc_
  integer :: ieig, local_ieig
  integer :: dest, tag
  real(DP), allocatable :: tmp_u_r(:,:), tmp_u_l(:,:)
 
  call this%set_mpi_comm(comm=comm)
  if (this%is_master) then
    call this%evecs_open_read_bin(rootname)
    call this%evecs_read_header_bin()
  end if
  call this%broadcast_header(comm=comm)
  call this%setup_paral_distribution(neig=neig_read, comm=comm)
  with_Avc_ = .true.
  if (present(with_Avc)) with_Avc_ = with_Avc
  call this%alloc(with_Avc=with_Avc_)
  ! Allocate temporary array for storage
  if (this%is_master) then
    allocate(tmp_u_r (this%usize, this%meig))
    if (.not. this%tda) then
      allocate(tmp_u_l (this%usize, this%meig))
    end if
  end if
  ! Read eigenvectors
  do ieig=1,this%neig
    if (this%is_master) then
      call this%evecs_read_single_eigenvector_bin(ieig_evec=1, ieig_eval=ieig)
    end if
    ! Communicate the eigenvector to the appropriate node
    dest = this%who_owns(ieig)
    local_ieig = this%local_ieigs(ieig)
    tag = 2000 + dest
    ! If destination is master, store in a temporary array
    if ((dest .eq. 0) .and. this%is_master) then
      tmp_u_r(:,local_ieig) = this%u_r(:,1)
    ! Otherwise, master sends the eigenvector
    else if (this%is_master) then
    ! node receives the eigenvector
    else if (this%inode .eq. dest) then
    end if
    if (.not. this%tda) then
      tag = 3000 + dest
      ! If destination is master, store in a temporary array
      if (dest .eq. 0) then
        tmp_u_l(:,local_ieig) = this%u_l(:,1)
      ! Otherwise, master sends the eigenvector
      else if (this%is_master) then
      ! node receives the eigenvector
      else if (this%inode .eq. dest) then
      end if
    end if
  end do
  if (this%is_master) then
    call this%evecs_close_bin()
    ! Copy eigenvector from temporary array
    do ieig=1,this%meig
      this%u_r(:,ieig) = tmp_u_r(:,ieig)
    end do
    if(allocated(tmp_u_r))then;deallocate(tmp_u_r);endif
    if (.not. this%tda) then
      do ieig=1,this%meig
        this%u_l(:,ieig) = tmp_u_l(:,ieig)
      end do
      if(allocated(tmp_u_l))then;deallocate(tmp_u_l);endif
    end if
  end if
  if (with_Avc_) then
    call this%reshape_Avc()
  end if
  call this%broadcast_global_arrays()
  this%has_flat_Avc = .true.
 
end subroutine evecs_read_and_distribute_bin
! =========================================================================== !
!> Close the binary eigenvectors file(s)
subroutine evecs_close_bin(this)
  class(evecs_t), intent(inout) :: this
 
  call close_file(this%uur)
  if ( .not. this%tda) call close_file(this%uul)
  this%is_open = .false.
 
end subroutine evecs_close_bin
! =========================================================================== !
! Writing routines
! =========================================================================== !
!> Generic routine to write BSE eigenvectors, whether binary or HDF5
subroutine evecs_write(this, nwrite, fname)
  class(evecs_t), intent(inout) :: this
  integer, intent(inout), optional :: nwrite
  character(len=*), intent(in), optional :: fname
 
  if (this%use_hdf5) then
  else
    call this%evecs_write_bin(nwrite, fname=fname)
  end if
 
end subroutine evecs_write
! =========================================================================== !
!> Write the binary eigenvectors file(s)
subroutine evecs_write_bin(this, nwrite, fname)
  class(evecs_t), intent(inout) :: this
  integer, intent(inout), optional :: nwrite !< Number of eigenvectors to write
  character(len=*), intent(in), optional :: fname !< File name
  integer :: nwrite_
  integer :: ii,jj,ik
  integer :: ieig, ieig_local
  integer :: who_owns
  integer :: rank_r, rank_l
  integer :: tag
  logical :: io_r, io_l
  character(len=512) :: fname_
  character(len=512) :: fname_r, fname_l
  real(DP), allocatable :: u_r(:), u_l(:)
 
  if (present(fname)) then
    fname_ = fname
  else
    fname_ = this%default_fname_bin
  end if
! Who can do io?
! FHJ: we do i/o for the right and left eigenvectors using different MPI ranks.
! The hope is to get better load balance on a lustre FS. The optimal solution,
! however, is to use HDF5.
  rank_r = 0
  io_r = (this%inode==rank_r)
  rank_l = this%npes-1
  io_l = ((this%inode==rank_l) .and. (.not. this%tda))
  nwrite_ = this%neig
  if (present(nwrite)) then
    if (nwrite.lt.0) then
      nwrite = this%neig
    else
      nwrite_ = nwrite
    end if
  end if
  this%nwrite = nwrite_
  ! Temporary arrays for communication
  allocate(u_r (this%usize))
  if (.not. this%tda) then
    allocate(u_l (this%usize))
  end if
  ! Open the file we will write to and write header information
  if (io_r) then
    write(6,*)
    if (.not. this%tda) then
      !fname_r = "eigenvectors_r"
      fname_r = trim(fname_)//"_r"
      write(6,'(1x,a,i0,a)') 'Writing ',nwrite_, ' right eigenvectors to file "'//trim(fname_r)//'"'
    else
      !fname_r = "eigenvectors"
      fname_r = trim(fname_)
      write(6,'(1x,a,i0,a)') 'Writing ',nwrite_, ' eigenvectors to file "'//trim(fname_r)//'"'
    endif
    write(6,'(1x,a,i0)') 'Length of each vector: ', this%usize
    write(6,*)
    call open_file(unit=this%uur,file=trim(fname_r),form='unformatted',status='replace')
    write(this%uur) this%ns
    write(this%uur) this%nv
    write(this%uur) this%nc
    write(this%uur) this%nk
    write(this%uur) ((this%kpts(jj,ik),jj=1,3),ik=1,this%nk)
  endif
  if (io_l) then
    !fname_l = "eigenvectors_l"
    fname_r = trim(fname_)//"_l"
    write(6,*)
    write(6,'(1x,a,i0,a)') 'Writing ',nwrite_, ' left eigenvectors to file "eigenvectors_l"'
    write(6,'(1x,a,i0)') 'Length of each vector: ', this%usize
    write(6,*)
    call open_file(unit=this%uul,file=trim(fname_l),form='unformatted',status='replace')
    write(this%uul) this%ns
    write(this%uul) this%nv
    write(this%uul) this%nc
    write(this%uul) this%nk
    write(this%uul) ((this%kpts(jj,ik),jj=1,3),ik=1,this%nk)
  endif
  ! Loop over states to be written
  do ieig=1,nwrite_
    ! Figure out which processor and column state ieig belongs to
    if (this%is_distributed) then
      who_owns = this%who_owns(ieig)
      ieig_local = this%local_ieigs(ieig)
    else
      who_owns = 0
      ieig_local = ieig
    end if
    ! Get the coeffs for state ieig into A (on all processors)
    if (this%inode==who_owns) then
      u_r(:) = this%u_r(:,ieig_local)
      if (.not. this%tda) then
        u_l(:) = this%u_l(:,ieig_local)
      endif
    endif
    ! Write to file
    if (io_r) then
      write(this%uur) this%evals(ieig)
      write(this%uur) (u_r(ii),ii=1,this%usize)
    endif
    if (io_l) then
      write(this%uul) this%evals(ieig)
      write(this%uul) (u_l(ii),ii=1,this%usize)
    endif
  end do
  if (io_r) call close_file(this%uur)
  if (io_l) call close_file(this%uul)
  if(allocated(u_r))then;deallocate(u_r);endif
  if(allocated(u_l))then;deallocate(u_l);endif
 
end subroutine evecs_write_bin
! =========================================================================== !
! HDF5 writing routines
! =========================================================================== !
! =========================================================================== !
! Output routines
! =========================================================================== !
!> Print out the information read in the header
subroutine evecs_print_out_header_info(this, unt)
  class(evecs_t), intent(inout) :: this
  integer, intent(in) :: unt !< The unit for writing
  integer :: ik
 
  if (this%is_master) then
    ! Print out information found in file
    write(unt,'(a)')
    write(unt,'(a,4i5)') ' ns, nv, nc, nk = ',this%ns,this%nv,this%nc,this%nk
    write(unt,'(a,i8)') ' nmat = ',this%nmat
    write(unt,'(a)')
    write(unt,'(a)') 'kpoints follow:'
    do ik=1, this%nk
      write(unt,'(i5,3f10.5)') ik, this%kpts(:,ik)
    enddo
  end if
 
end subroutine evecs_print_out_header_info
! =========================================================================== !
!> Print out current dimensions of the object
subroutine evecs_print_out_dimensions(this, unt)
  class(evecs_t), intent(inout) :: this
  integer, intent(in) :: unt !< The unit for writing
 
  if (this%is_master) then
    write(unt,'(a)') ' ---------------------------------------'
    write(unt,'(a)') ' BSE eigenvectors dimensions'
    write(unt,'(a,4i5)') ' ns, nv, nc, nk = ',this%ns,this%nv,this%nc,this%nk
    write(unt,'(a,i8)') ' nmat = ',this%nmat
    write(unt,'(a,i8)') ' neig = ',this%neig
    write(unt,'(a,i8)') ' meig = ',this%meig
    write(unt,'(a)') ' ---------------------------------------'
    write(unt,'(a)')
  end if
 
end subroutine evecs_print_out_dimensions
! =========================================================================== !
end module evecs_m
