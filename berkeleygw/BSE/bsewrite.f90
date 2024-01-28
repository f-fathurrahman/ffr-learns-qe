!==============================================================================
!
! Routines:
!
! (1) bsewrite() Originally By MLT Last Modified 7/1/2008 (JRD)
!
! input: xct types
! bsedbody,bsedhead,bsedwing,bsex
! output: binary files "bsedmat", "bsexmat"
!
! Write out all interaction matrices elements in "bsedmat", "bsexmat"
!
! Since all PEs write into the same units, the access must be organized:
! one PE at a time
!
! It appears the imatrix code is: 1 dhead, 2 dwing, 3 dbody, 4 x
!
! (2) bse_hdf5_write() Orgiginally By JRD Last Modified 4/27/2012 (JRD)
!
! hdf5 writer routine for bsemat.h5 files
!
!=============================================================================

module bsewrite_m
  use global_m
  use io_utils_m
  use kernel_io_m
  use wfn_rho_vxc_io_m
  implicit none
  public :: bsewrite
  private
contains
subroutine bsewrite(xct,iownsize,bsedbody,bsedhead,bsex,kg,kp,gvec,syms,crys,&
    bsedwing,bset)
  type (xctinfo), intent(in) :: xct
  type (grid), intent(in) :: kg
  integer, intent(in) :: iownsize
  real(DP), intent(in) :: bsedbody(:,:,:), bsedhead(:,:,:), &
    bsex(:,:,:) !< (iownsize*peinf%myown,xct%nspin,xct%nspin)
  type (kpoints), intent(in) :: kp
  type (gspace), intent(in) :: gvec
  type (symmetry), intent(in) :: syms
  type (crystal), intent(in) :: crys
  real(DP), intent(in), optional :: &
    bsedwing(:,:,:) !< (iownsize*peinf%myown,xct%nspin,xct%nspin)
  real(DP), intent(in), optional :: &
    bset(:,:,:) !< (iownsize*peinf%myown,xct%nspin,xct%nspin)
  real(DP), allocatable :: bsemt(:,:,:,:,:), bsemtt(:,:,:,:,:)
  real(DP) :: bsem
  integer :: nmatrices,imatrix,iunit,error,version
  integer :: ic,icp,ik,ikp,is1,is2,iv,ivp,it
  real(DP) :: bsedh_fro2, bsedw_fro2, bsedb_fro2, bsex_fro2, bset_fro2
  real(DP) :: bsedh_fro, bsedw_fro, bsedb_fro, bsex_fro, bset_fro
  real(DP) :: approx_sz
  type(progress_info) :: prog_info
  type(kernel_header_t) :: kern
 
  approx_sz = (xct%n1b_co*xct%n2b_co*xct%nkpt_co*xct%nspin + 8d0)**2 * 8d0
  approx_sz = 1 * approx_sz / (1024d0 * 1024d0)
  if (peinf%inode .eq. 0) then
    write(6,'(1x,a)') 'Expected size of the BSE matrices:'
    if (xct%use_hdf5) then
      write(6,'(1x,a,f0.3,a)') '  bsemat.h5: ', approx_sz*4, ' MB'
    else
      write(6,'(1x,a,f0.3,a)') '  bsexmat: ', approx_sz, ' MB'
      write(6,'(1x,a,f0.3,a)') '  bsedmat: ', approx_sz*3, ' MB'
    endif
    write(6,*)
  endif
  ! FHJ: Min. number of matrices we save is 3 => X + W_body + W_head
  nmatrices = 3
  if (present(bsedwing)) nmatrices = nmatrices + 1
  if (present(bset)) nmatrices = nmatrices + 1
  ! FHJ: Compue local contribution to the Frobenius norm of the matrices
  bsedh_fro2 = sum(((bsedhead(:,:,:))**2))
  if (present(bsedwing)) bsedw_fro2 = sum(((bsedwing(:,:,:))**2))
  bsedb_fro2 = sum(((bsedbody(:,:,:))**2))
  bsex_fro2 = sum(((bsex(:,:,:))**2))
  if (present(bset)) bset_fro2 = sum(((bset(:,:,:))**2))
  version = VER_BSE_FORT
  call init_mf_header_from_types(kern%mf, 'KER', 1, version, kp, gvec, syms, crys)
  call xctinfo_to_kernel_header(xct, kg%f, kern, 3)
  if (present(bsedwing)) kern%nmat = kern%nmat + 1
  if (present(bset)) kern%nmat = kern%nmat + 1
    allocate(bsemt (xct%nspin,xct%nspin,xct%n1b_co,xct%n2b_co,xct%nkpt_co))
    allocate(bsemtt (xct%nspin,xct%nspin,xct%n1b_co,xct%n2b_co,xct%nkpt_co))
    if(peinf%inode .eq. 0 ) then
      call open_file(unit=11,file='bsedmat',form='unformatted',status='replace')
      call open_file(unit=12,file='bsexmat',form='unformatted',status='replace')
      if (present(bset)) &
        call open_file(unit=13,file='bsetmat',form='unformatted',status='replace')
    endif
    if (present(bsedwing)) then
      kern%nmat = 3
    else
      kern%nmat = 2
    endif
    call write_binary_kernel_header(11, kern)
    kern%nmat = 1
    call write_binary_kernel_header(12, kern)
    if (present(bset)) then
      kern%nmat = 1
      call write_binary_kernel_header(13, kern)
    endif
    ! FHJ: this is to generate nice output / time estimate
    call progress_init(prog_info, 'writing out BSE matrices', 'records', &
      xct%nkpt_co*nmatrices*xct%n2b_co*xct%n1b_co)
    do imatrix = 1, nmatrices
      do ikp=1,xct%nkpt_co
        do icp=1,xct%n2b_co
          do ivp=1,xct%n1b_co
            call progress_step(prog_info)
            bsemt=0.d0
            do ik=1,xct%nkpt_co
              do ic=1,xct%n2b_co
                do iv=1,xct%n1b_co
                  if (xct%icpar .eq. 0) then
                    it = peinf%wown(1,1,ikp,1,1,ik)
                    if (it .ne. 0) then
                      it = peinf%wown(1,1,ikp,1,1,ik) + xct%n1b_co*xct%n1b_co*xct%n2b_co*(icp-1) &
                        + xct%n1b_co*xct%n1b_co*(ic-1) + xct%n1b_co*(ivp-1) + iv -1
                    endif
                  else if (xct%ivpar .eq. 0) then
                    it = peinf%wown(1,icp,ikp,1,ic,ik)
                    if (it .ne. 0) then
                      it = it + xct%n1b_co*(ivp-1) + iv -1
                    endif
                  else
                    it = peinf%wown(ivp,icp,ikp,iv,ic,ik)
                  endif
                  if (it .ne. 0) then
                    do is1=1,xct%nspin
                      do is2=1,xct%nspin
                        if (present(bsedwing)) then
                          if ( imatrix .eq. 1) then
                            bsem=bsedhead(it,is1,is2)
                          else if ( imatrix .eq. 2) then
                            bsem=bsedwing(it,is1,is2)
                          else if ( imatrix .eq. 3) then
                            bsem=bsedbody(it,is1,is2)
                          else if ( imatrix .eq. 4) then
                            bsem=bsex(it,is1,is2)
                          endif
                        else if (present(bset)) then
                          if ( imatrix .eq. 1) then
                            bsem=bsedhead(it,is1,is2)
                          else if ( imatrix .eq. 2) then
                            bsem=bsedbody(it,is1,is2)
                          else if ( imatrix .eq. 3) then
                            bsem=bsex(it,is1,is2)
                          else if ( imatrix .eq. 4) then
                            bsem=bset(it,is1,is2)
                          endif
                        endif
                        bsemt(is2,is1,iv,ic,ik) = bsem
                      enddo !is2
                    enddo !is1
                  endif
                enddo !iv
              enddo !ic
            enddo !ik
            bsemtt=0D0
            call timacc(75,1)
            call timacc(75,2)
            call timacc(76,1)
            if (peinf%inode .eq. 0) then
              if (present(bsedwing)) then
                if ( imatrix .lt. 4) then
                  iunit = 11
                else if (imatrix .eq. 4) then
                  iunit = 12
                endif
              else
                if ( imatrix .lt. 3) then
                  iunit = 11
                else if (imatrix .eq. 3) then
                  iunit = 12
                else if (imatrix .eq. 4) then
                  iunit = 13
                endif
              endif
              write(iunit) ikp,icp,ivp, bsemt(:,:,:,:,:)
            endif
            call timacc(76,2)
          enddo !ivp
        enddo !icp
      enddo !ikp
    enddo !imatrix
    call progress_free(prog_info)
    if(peinf%inode .eq. 0 ) then
      call close_file(11)
      call close_file(12)
      if (present(bset)) call close_file(13)
    endif
! FHJ: Reduce and print the Frobenius norm from all matrices. We use Frobenius
! instead of the max norm because the Frobenius norm is invariant under unitary transf.
  bsedh_fro = bsedh_fro2
  if (present(bsedwing)) &
    bsedw_fro = bsedw_fro2
  bsedb_fro = bsedb_fro2
  bsex_fro = bsex_fro2
  if (present(bset)) &
    bset_fro = bset_fro2
  if (peinf%inode.eq.0) then
    write(6,*)
    write(6,'(1x,a)') 'Frobenius norm of the matrices per spin:'
    write(6,'(1x,a,es21.12e4)') '- Head : ', sqrt(bsedh_fro)/xct%nspin
    if (present(bsedwing)) &
      write(6,'(1x,a,es21.12e4)') '- Wing : ', sqrt(bsedw_fro)/xct%nspin
    write(6,'(1x,a,es21.12e4)') '- Body : ', sqrt(bsedb_fro)/xct%nspin
    write(6,'(1x,a,es21.12e4)') '- X    : ', sqrt(bsex_fro)/xct%nspin
    if (present(bset)) &
      write(6,'(1x,a,es21.12e4)') '- T    : ', sqrt(bset_fro)/xct%nspin
    write(6,*)
  endif
  if(.not. xct%use_hdf5) then
    if(allocated(bsemt))then;deallocate(bsemt);endif
    if(allocated(bsemtt))then;deallocate(bsemtt);endif
  endif
 
  return
end subroutine bsewrite
!============================================================================
!
! (2) bse_hdf5_write
!
!============================================================================
end module bsewrite_m
