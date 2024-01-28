!===============================================================================
!
! Routines:
!
!
!===============================================================================

module read_rho_vxc_m
  use global_m
  use gmap_m
  use wfn_rho_vxc_io_m
  implicit none
  private
  public :: read_vxc, read_rho
contains
  !---------------------------------------------------------------------------------------------------
  !> Read in the exchange-correlation potential and store in array sig%vxc
  subroutine read_vxc(sig, gvec, kp, syms, crys, isrti, isrt, vxc_type, tolerant_value)
    type(siginfo), intent(inout) :: sig
    type(gspace), intent(in) :: gvec
    type(kpoints), intent(in) :: kp
    type(symmetry), intent(in) :: syms
    type(crystal), intent(in) :: crys
    integer, intent(in) :: isrti(:) !< (gvec%ng)
    integer, intent(in) :: isrt(:) !< (gvec%ng)
    integer, intent(in) :: vxc_type
    logical, optional, intent(in) :: tolerant_value
    character*3 :: sheader
    integer :: iflavor, ig, ispin
    type(gspace) :: gvec_dummy
    type(crystal) :: crys_dummy
    type(symmetry) :: syms_dummy
    type(kpoints) :: kp_dummy
    real(DP) :: discrepancy
    logical :: tolerant_
   
    call logit('Reading VXC')
    if (vxc_type .eq. 1) then
      allocate(sig%vxc (gvec%ng,kp%nspin))
    else if (vxc_type .eq. 2) then
      allocate(sig%vxc2 (gvec%ng,kp%nspin))
    else if (vxc_type .eq. 3) then
      allocate(sig%dvxc (gvec%ng,kp%nspin))
    else
      call die("VXC type can only be 1 (VXC), 2 (VXC2), and 3 (dVXC).")
    endif
    if (vxc_type .eq. 1) then
      if(peinf%inode == 0) call open_file(96,file='VXC',form='unformatted',status='old')
    else if (vxc_type .eq. 2) then
      if(peinf%inode == 0) call open_file(96,file='VXC2',form='unformatted',status='old')
    else if (vxc_type .eq. 3) then
      if(peinf%inode == 0) call open_file(96,file='dVXC',form='unformatted',status='old')
    endif
    tolerant_ = .false.
    if(present(tolerant_value)) tolerant_ = tolerant_value
    sheader = 'VXC'
    iflavor = 0
    call read_binary_header_type(96, sheader, iflavor, kp_dummy, gvec_dummy, syms_dummy, &
         crys_dummy, warn = .false.)
    call check_header('WFN_inner', kp, gvec, syms, crys, 'VXC', kp_dummy, gvec_dummy, &
         syms_dummy, crys_dummy, is_wfn = .false., tolerant=tolerant_)
    allocate(gvec_dummy%components (3, gvec_dummy%ng))
    call read_binary_gvectors(96, gvec_dummy%ng, gvec_dummy%ng, gvec_dummy%components)
    do ig = 1, gvec%ng
      if(any(gvec_dummy%components(:,isrt(ig)) .ne. gvec%components(:,ig))) call die("gvec mismatch in VXC")
    enddo
    if(associated(gvec_dummy%components))then;deallocate(gvec_dummy%components);nullify(gvec_dummy%components);endif
    if (vxc_type .eq. 1) then
      call read_binary_data(96, gvec_dummy%ng, gvec_dummy%ng, kp%nspin, sig%vxc, gindex = isrti)
    else if (vxc_type .eq. 2) then
      call read_binary_data(96, gvec_dummy%ng, gvec_dummy%ng, kp%nspin, sig%vxc2, gindex = isrti)
    else if (vxc_type .eq. 3) then
      call read_binary_data(96, gvec_dummy%ng, gvec_dummy%ng, kp%nspin, sig%dvxc, gindex = isrti)
    endif
    if(peinf%inode == 0) then
      call close_file(96)
      if (vxc_type .eq. 1) then
        do ispin = 1, kp%nspin
          discrepancy = check_field_is_real(sig%vxc(:, ispin), gvec)
          if(discrepancy > TOL_Zero) then
            write(0,*) 'WARNING: VXC is not real in real space, with discrepancy ', discrepancy, ' for spin ', ispin
          endif
        enddo
      else if (vxc_type .eq. 2) then
        do ispin = 1, kp%nspin
          discrepancy = check_field_is_real(sig%vxc2(:, ispin), gvec)
          if(discrepancy > TOL_Zero) then
            write(0,*) 'WARNING: VXC2 is not real in real space, with discrepancy ', discrepancy, ' for spin ', ispin
          endif
        enddo
      else if (vxc_type .eq. 3) then
        write (0,*) 'NOTE: dVXC is complex in general.'
      endif
    endif
    call dealloc_header_type(sheader, crys_dummy, kp_dummy)
   
  end subroutine read_vxc
  !---------------------------------------------------------------------------------------------------
  !> Read in the charge density and store in array wpg%rho (formerly known as CD95)
  !! CD95 Ref: http://www.nature.com/nature/journal/v471/n7337/full/nature09897.html
  subroutine read_rho(wpg, gvec, kp, syms, crys, isrti, isrt, check_filename, tolerant_value)
    type(wpgen), intent(inout) :: wpg
    type(gspace), intent(in) :: gvec
    type(kpoints), intent(in) :: kp
    type(symmetry), intent(in) :: syms
    type(crystal), intent(in) :: crys
    integer, intent(in) :: isrti(:) !< (gvec%ng)
    integer, intent(in) :: isrt(:) !< (gvec%ng)
    character(len=*), intent(in) :: check_filename !< This is the file
                                            !< against which the header is
                                            !< checked
    logical, optional, intent(in) :: tolerant_value
    character*3 :: sheader
    integer :: iflavor, ig, ispin
    type(gspace) :: gvec_dummy
    type(crystal) :: crys_dummy
    type(symmetry) :: syms_dummy
    type(kpoints) :: kp_dummy
    real(DP) :: discrepancy
    logical :: tolerant_
   
    call logit('Reading RHO')
    allocate(wpg%rho (gvec%ng,kp%nspin))
    if(peinf%inode == 0) call open_file(95,file='RHO',form='unformatted',status='old')
    tolerant_ = .false.
    if(present(tolerant_value)) tolerant_ = tolerant_value
    sheader = 'RHO'
    iflavor = 0
    call read_binary_header_type(95, sheader, iflavor, kp_dummy, gvec_dummy, &
        syms_dummy, crys_dummy, warn = .false.)
    call check_header(trim(check_filename), kp, gvec, syms, crys, 'RHO', kp_dummy, &
        gvec_dummy, syms_dummy, crys_dummy, is_wfn = .false., tolerant=tolerant_)
    allocate(gvec_dummy%components (3, gvec_dummy%ng))
    call read_binary_gvectors(95, gvec_dummy%ng, gvec_dummy%ng, gvec_dummy%components)
    do ig = 1, gvec%ng
      if(any(gvec_dummy%components(:,isrt(ig)) .ne. gvec%components(:,ig))) call die("gvec mismatch in RHO")
    enddo
    if(associated(gvec_dummy%components))then;deallocate(gvec_dummy%components);nullify(gvec_dummy%components);endif
    call read_binary_data(95, gvec_dummy%ng, gvec_dummy%ng, kp%nspin, wpg%rho, gindex = isrti)
    if(peinf%inode == 0) call close_file(95)
    ! otherwise if nspin == 1, the 2 component may be uninitialized to NaN
    wpg%wpsq(1:2) = 0d0
    wpg%nelec(1:2) = 0d0
    ! since they are sorted, if G = 0 is present, it is the first one
    if(any(gvec%components(1:3, 1) /= 0)) call die("gvectors for RHO must include G = 0")
    ! otherwise, the code below will not do what we think it does
    do ispin=1,kp%nspin
      wpg%nelec(ispin)=dble(wpg%rho(1,ispin))
      wpg%wpsq(ispin)=ryd*ryd*16.0d0*PI_D*wpg%nelec(ispin)/crys%celvol
    enddo
    ! This is unacceptable because it means the number of electrons is negative,
    ! and the plasma frequency will be imaginary!
    if(any(wpg%nelec(1:kp%nspin) < TOL_Zero)) then
      write(0,*) wpg%nelec(:)
      call die("Charge density in RHO has negative part for G=0", only_root_writes = .true.)
    endif
    if(peinf%inode == 0) then
      do ispin = 1, kp%nspin
        discrepancy = check_field_is_real(wpg%rho(:, ispin), gvec)
        if(discrepancy > TOL_Zero) then
          write(0,*) 'WARNING: RHO is not real in real space, with discrepancy ', discrepancy, ' for spin ', ispin
        endif
      enddo
    endif
    call dealloc_header_type(sheader, crys_dummy, kp_dummy)
   
  end subroutine read_rho
  !---------------------------------------------------------------------------------------------------
  !> RHO and VXC must be real in real space, i.e. c(G) - c(-G)* = 0
  real(DP) function check_field_is_real(field, gvec)
    real(DP), intent(in) :: field(:)
    type(gspace), intent(in) :: gvec
    integer :: ig, umklapp(3)
    real(DP) :: diff
    type(symmetry) :: syms_inv
    integer, allocatable :: ind(:), identity(:)
    real(DP), allocatable :: phase(:)
   
    syms_inv%ntran = 1
    syms_inv%mtrx = 0
    syms_inv%tnp = 0d0
    syms_inv%mtrx(1:3, 1:3, 1) = reshape((/-1, 0, 0, 0, -1, 0, 0, 0, -1/), shape(syms_inv%mtrx(:,:,1)))
    umklapp(1:3) = 0
    allocate(phase (gvec%ng))
    allocate(ind (gvec%ng))
    allocate(identity (gvec%ng))
    do ig = 1, gvec%ng
      identity(ig) = ig
    enddo
    call gmap(gvec, syms_inv, gvec%ng, 1, umklapp, identity, identity, ind, phase, die_outside_sphere = .true.)
    if(any(abs(phase(1:gvec%ng) - 1d0) > TOL_Zero)) call die("non-unity phase in check_field_is_real")
    if(any(ind(1:gvec%ng) == 0)) call die("ind array from gmap has a zero")
    check_field_is_real = 0d0
    do ig = 1, gvec%ng
      diff = field(ig) - (field(ind(ig))*phase(ig))
      check_field_is_real = max(check_field_is_real, abs(diff))
    enddo
    if(allocated(phase))then;deallocate(phase);endif
    if(allocated(ind))then;deallocate(ind);endif
    if(allocated(identity))then;deallocate(identity);endif
   
  end function check_field_is_real
end module read_rho_vxc_m
