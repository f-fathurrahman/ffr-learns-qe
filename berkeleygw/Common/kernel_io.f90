!>=========================================================================
!!
!! Module:
!!
!! (1) kernel_io_m Originally by FHJ Last Modified 10/07/2013 (FHJ)
!!
!! Routines to read and write kernel files (bsedmat, ...)
!! Inspired on wfn_rho_vxc_io.F90.
!!
!!=========================================================================

module kernel_io_m
  use global_m
  use wfn_rho_vxc_io_m
  implicit none
  private
 !> For library usage, do not make global_m contents available
 !! to avoid namespace clashes.
  public :: &
    read_binary_kernel_header, &
    write_binary_kernel_header, &
    read_format_kernel_header, &
    write_format_kernel_header, &
    read_kernel_header, &
    write_kernel_header, &
    xctinfo_to_kernel_header, &
    check_xctinfo_kernel_header
contains
!===============================================================================
! FHJ: Below are the routines for kernel*mat files
!===============================================================================
!read_formatted_kernel_header
!=========================================================================
!
! Included from file kernel_io.F90.
! You might be expected to understand this. --FHJ
!
!=========================================================================
!> Defines a subroutine with the template {read,write}_{formatted,binary}_kernel_header
subroutine read_format_kernel_header(iunit, kernel)
  integer, intent(in) :: iunit !< unit number
  type(kernel_header_t), intent(out) :: kernel !< kernel_header_t type
  integer :: ik
 
  call read_mf_header(iunit, kernel%mf)
  if (kernel%mf%sheader/='KER') then
    write(0,*) 'ERROR: header mismatch (got "'//kernel%mf%sheader//'", expected "KER")'
    call die('Input file is not from a kernel calculation (header="'//kernel%mf%sheader//'")', &
      only_root_writes=.true.)
  endif
  if (peinf%inode==0) then
    ! General information
    read(iunit , *) kernel%iscreen, kernel%icutv, kernel%ecuts
    ! Variables specific to kernel files: kpts
    read(iunit , *) kernel%nk
    allocate(kernel%kpts (3,kernel%nk))
    do ik = 1, kernel%nk
      read(iunit , *) kernel%kpts(1:3, ik)
    enddo
    ! Variables specific to kernel files: everything else
    read(iunit , *) kernel%ns, kernel%nspinor, kernel%nvb, kernel%ncb, kernel%n1b, kernel%n2b
    read(iunit , *) kernel%theory, kernel%nmat, kernel%storage, kernel%nblocks
    ! Empty records: you can use these slots in the future to extend the file format
    read(iunit , *)
    read(iunit , *)
    read(iunit , *)
    read(iunit , *)
    read(iunit , *)
  endif
 
end subroutine read_format_kernel_header
! these undefs prevent lots of cpp warnings
!write_formatted_kernel_header
!=========================================================================
!
! Included from file kernel_io.F90.
! You might be expected to understand this. --FHJ
!
!=========================================================================
!> Defines a subroutine with the template {read,write}_{formatted,binary}_kernel_header
subroutine write_format_kernel_header(iunit, kernel)
  integer, intent(in) :: iunit !< unit number
  type(kernel_header_t), intent(in) :: kernel !< kernel_header_t type
  integer :: ik
 
  call write_mf_header(iunit, kernel%mf)
  if (peinf%inode==0) then
    ! General information
    write(iunit , *) kernel%iscreen, kernel%icutv, kernel%ecuts
    ! Variables specific to kernel files: kpts
    write(iunit , *) kernel%nk
    do ik = 1, kernel%nk
      write(iunit , *) kernel%kpts(1:3, ik)
    enddo
    ! Variables specific to kernel files: everything else
    write(iunit , *) kernel%ns, kernel%nspinor, kernel%nvb, kernel%ncb, kernel%n1b, kernel%n2b
    write(iunit , *) kernel%theory, kernel%nmat, kernel%storage, kernel%nblocks
    ! Empty records: you can use these slots in the future to extend the file format
    write(iunit , *)
    write(iunit , *)
    write(iunit , *)
    write(iunit , *)
    write(iunit , *)
  endif
 
end subroutine write_format_kernel_header
! these undefs prevent lots of cpp warnings
!read_binary_kernel_header
!=========================================================================
!
! Included from file kernel_io.F90.
! You might be expected to understand this. --FHJ
!
!=========================================================================
!> Defines a subroutine with the template {read,write}_{formatted,binary}_kernel_header
subroutine read_binary_kernel_header(iunit, kernel)
  integer, intent(in) :: iunit !< unit number
  type(kernel_header_t), intent(out) :: kernel !< kernel_header_t type
  integer :: ik
 
  call read_mf_header(iunit, kernel%mf)
  if (kernel%mf%sheader/='KER') then
    write(0,*) 'ERROR: header mismatch (got "'//kernel%mf%sheader//'", expected "KER")'
    call die('Input file is not from a kernel calculation (header="'//kernel%mf%sheader//'")', &
      only_root_writes=.true.)
  endif
  if (peinf%inode==0) then
    ! General information
    read(iunit ) kernel%iscreen, kernel%icutv, kernel%ecuts
    ! Variables specific to kernel files: kpts
    read(iunit ) kernel%nk
    allocate(kernel%kpts (3,kernel%nk))
    do ik = 1, kernel%nk
      read(iunit ) kernel%kpts(1:3, ik)
    enddo
    ! Variables specific to kernel files: everything else
    read(iunit ) kernel%ns, kernel%nspinor, kernel%nvb, kernel%ncb, kernel%n1b, kernel%n2b
    read(iunit ) kernel%theory, kernel%nmat, kernel%storage, kernel%nblocks
    ! Empty records: you can use these slots in the future to extend the file format
    read(iunit )
    read(iunit )
    read(iunit )
    read(iunit )
    read(iunit )
  endif
 
end subroutine read_binary_kernel_header
! these undefs prevent lots of cpp warnings
!write_binary_kernel_header
!=========================================================================
!
! Included from file kernel_io.F90.
! You might be expected to understand this. --FHJ
!
!=========================================================================
!> Defines a subroutine with the template {read,write}_{formatted,binary}_kernel_header
subroutine write_binary_kernel_header(iunit, kernel)
  integer, intent(in) :: iunit !< unit number
  type(kernel_header_t), intent(in) :: kernel !< kernel_header_t type
  integer :: ik
 
  call write_mf_header(iunit, kernel%mf)
  if (peinf%inode==0) then
    ! General information
    write(iunit ) kernel%iscreen, kernel%icutv, kernel%ecuts
    ! Variables specific to kernel files: kpts
    write(iunit ) kernel%nk
    do ik = 1, kernel%nk
      write(iunit ) kernel%kpts(1:3, ik)
    enddo
    ! Variables specific to kernel files: everything else
    write(iunit ) kernel%ns, kernel%nspinor, kernel%nvb, kernel%ncb, kernel%n1b, kernel%n2b
    write(iunit ) kernel%theory, kernel%nmat, kernel%storage, kernel%nblocks
    ! Empty records: you can use these slots in the future to extend the file format
    write(iunit )
    write(iunit )
    write(iunit )
    write(iunit )
    write(iunit )
  endif
 
end subroutine write_binary_kernel_header
! these undefs prevent lots of cpp warnings
!read_kernel_header
!=========================================================================
!
! Included from file kernel_io.F90.
! You might be expected to understand this. --FHJ
!
!=========================================================================
!> Automatically call {read,write}_{formatted,binary}_kernel_header
subroutine read_kernel_header(iunit, kernel)
  integer, intent(in) :: iunit !< unit number
  type(kernel_header_t), intent(out) :: kernel !< kernel_header_t type
  character(len=16) :: fmt_str
  logical :: is_fmt = .false.
 
  if (peinf%inode==0) then
    inquire(unit=iunit, form=fmt_str)
    if (TRUNC(fmt_str)=='FORMATTED') then
      is_fmt = .true.
    else if (TRUNC(fmt_str)/='UNFORMATTED') then
      call die('Unknown value for formatted string: '//TRUNC(fmt_str), &
        only_root_writes=.true.)
    endif
  endif
  ! FHJ: No need to bcast is_fmt because only inode==0 reads the file.
  if (is_fmt) then
    call read_format_kernel_header(iunit, kernel)
  else
    call read_binary_kernel_header(iunit, kernel)
  endif
 
end subroutine read_kernel_header
! these undefs prevent lots of cpp warnings
!write_kernel_header
!=========================================================================
!
! Included from file kernel_io.F90.
! You might be expected to understand this. --FHJ
!
!=========================================================================
!> Automatically call {read,write}_{formatted,binary}_kernel_header
subroutine write_kernel_header(iunit, kernel)
  integer, intent(in) :: iunit !< unit number
  type(kernel_header_t), intent(in) :: kernel !< kernel_header_t type
  character(len=16) :: fmt_str
  logical :: is_fmt = .false.
 
  if (peinf%inode==0) then
    inquire(unit=iunit, form=fmt_str)
    if (TRUNC(fmt_str)=='FORMATTED') then
      is_fmt = .true.
    else if (TRUNC(fmt_str)/='UNFORMATTED') then
      call die('Unknown value for formatted string: '//TRUNC(fmt_str), &
        only_root_writes=.true.)
    endif
  endif
  ! FHJ: No need to bcast is_fmt because only inode==0 reads the file.
  if (is_fmt) then
    call write_format_kernel_header(iunit, kernel)
  else
    call write_binary_kernel_header(iunit, kernel)
  endif
 
end subroutine write_kernel_header
! these undefs prevent lots of cpp warnings
  !> Populate the non-MF part of a kernel_header_t type. We assume that the MF
  !! part, i.e., kernel%mf, is already set (although we do overwrite
  !! kernel%mf%sheader='BSE' to play safe).
  !! Unlike WFN files, you can`t specify the flavor manually, it always matches 1 .
  subroutine xctinfo_to_kernel_header(xct, kpts, kernel, nmat)
    type(xctinfo), intent(in) :: xct
    real(DP), intent(in) :: kpts(:,:) !< (3, nk)
    type(kernel_header_t), intent(inout) :: kernel
    integer, intent(in) :: nmat !< 1 for kernelxmat, 3 for kerneldmat
   
    ! Generic header
    kernel%mf%sheader = 'KER'
    kernel%version = VER_BSE_HDF5
    kernel%iflavor = 1
    ! General paramters
    kernel%iscreen = xct%iscreen
    kernel%icutv = xct%icutv
    kernel%ecuts = xct%ecute
    kernel%ecutg = xct%ecutg
    kernel%efermi = xct%efermi
    kernel%theory = xct%theory
    kernel%nblocks = 1
    if (xct%extended_kernel) kernel%nblocks = 4
    kernel%storage = 0 ! Hard coded for now
    kernel%nmat = nmat
    kernel%energy_loss = xct%energy_loss
    ! K-point stuff
    kernel%nk = xct%nkpt_co
    allocate(kernel%kpts (3,kernel%nk))
    kernel%kpts = kpts(1:3, 1:kernel%nk)
    kernel%kgrid = kernel%mf%kp%kgrid
    kernel%qflag = xct%qflag
    kernel%exciton_Q_shift = xct%finiteq
    kernel%patched_sampling = xct%patched_sampling_co
    ! Bands stuff
    kernel%nvb = xct%nvb_co
    kernel%ncb = xct%ncb_co
    kernel%n1b = xct%n1b_co
    kernel%n2b = xct%n2b_co
    kernel%ns = xct%nspin
    kernel%nspinor = xct%nspinor
   
  end subroutine xctinfo_to_kernel_header
  subroutine check_xctinfo_kernel_header(fname, xct, kernel)
    character(len=*), intent(in) :: fname !< file name
    type(xctinfo), intent(in) :: xct
    type(kernel_header_t), intent(in) :: kernel
    integer :: nblocks
   
    ! Generic variables
    ! FHJ: absorption doesn`t care about WFN cutoff (for now!)
    !call check_R('WFN cutoff', xct%ecutg, kernel%mf%kp%ecutwfc)
    call check_I('screening flag', xct%iscreen, kernel%iscreen)
    call check_I('truncation flag', xct%icutv, kernel%icutv)
    ! FHJ: absorption doesn`t care about epsilon cutoff
    !call check_R('epsilon cutoff', xct%ecute, kernel%ecuts)
    ! Specific to bsemat
    call check_I('# of k-point', xct%nkpt_co, kernel%nk)
    call check_I('# of spins', xct%nspin, kernel%ns)
    call check_I('# of spinor components', xct%nspinor, kernel%nspinor)
    call check_I('# of val. bands', xct%nvb_co, kernel%nvb)
    call check_I('# of cond. bands', xct%ncb_co, kernel%ncb)
    call check_I('# of bands 1', xct%n1b_co, kernel%n1b)
    call check_I('# of bands 2', xct%n2b_co, kernel%n2b)
    call check_I('theory level', 0, kernel%theory) ! Hard coded for now
    nblocks = 1
    if (xct%extended_kernel) nblocks = 4
    call check_I('number of transition blocks', nblocks, kernel%nblocks)
    call check_I('storage format', 0, kernel%storage) ! Hard coded for now
   
    contains
      subroutine check_I(label, ref, got)
        character(len=*), intent(in) :: label
        integer, intent(in) :: ref
        integer, intent(in) :: got
       
        if (ref/=got) then
          if (peinf%inode==0) then
            write(0,'(1x,3a)') 'ERROR: incompatible values found in file "', fname, '".'
            write(0,'(1x,2a)') 'Quantity: ', label
            write(0,'(1x,a,i0)') 'Expected: ', ref
            write(0,'(1x,a,i0)') 'Got: ', got
          endif
          call die('Incompatible values in file "'+fname+'".', &
            only_root_writes=.true.)
        endif
       
      end subroutine check_I
      subroutine check_R(label, ref, got)
        character(len=*), intent(in) :: label
        real(DP), intent(in) :: ref
        real(DP), intent(in) :: got
        real(DP) :: tol=TOL_Small
       
        if (dabs(ref-got)>tol) then
          if (peinf%inode==0) then
            write(0,'(1x,3a)') 'ERROR: incompatible values found in file "', fname, '".'
            write(0,'(1x,2a)') 'Quantity: ', label
            write(0,'(1x,a,f0.8)') 'Expected: ', ref
            write(0,'(1x,a,f0.8)') 'Got: ', got
          endif
          call die('Incompatible values in file "'+fname+'".', &
            only_root_writes=.true.)
        endif
       
      end subroutine check_R
  end subroutine check_xctinfo_kernel_header
end module kernel_io_m
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
