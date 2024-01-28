!=================================================================================
!
! Module write_matrix_m
!
! (1) write_matrix_d() Originally by JRD Last Modified 5/1/2008 (JRD)
!
! This program writes a distributed matrix like chimat or epsmat to file.
!
! (2) write_matrix_f() Originally by JRD Last Modified 2/5/2009 (CHP)
!
! Modification of write_matrix_d for full-frequency.
!
!=================================================================================

module write_matrix_m
  use, intrinsic :: iso_c_binding
  use global_m
  use scalapack_m
  use io_utils_m
  use timing_m, only: timing => common_timing
  implicit none
  private
  public :: &
    write_matrix_d, &
    write_matrix_d_sub, &
    write_matrix_f
contains
!===================================================================================
subroutine write_matrix_d(scal,matrix,nmtx,iunit)
  type(scalapack), intent(in) :: scal
  real(DP), intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iunit
  integer :: ii, jj
  type(progress_info) :: prog_info !< a user-friendly progress report
 
  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nmtx, iunit
    write(6,*)
  endif
  if (peinf%inode .eq. 0) then
    call progress_init(prog_info, 'writing matrix', 'column', nmtx)
    do jj = 1, nmtx
      call progress_step(prog_info, jj)
      write(iunit) (matrix(ii, jj), ii = 1, nmtx)
    enddo
    call progress_free(prog_info)
  endif
 
  return
end subroutine write_matrix_d
subroutine write_matrix_d_sub(scal,matrix,nmtx,iunit,neig_sub)
  type(scalapack), intent(in) :: scal
  complex(DPC), intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iunit
  integer, intent(in), optional :: neig_sub
  integer :: ii, jj, nmtx_col
  type(progress_info) :: prog_info !< a user-friendly progress report
 
  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nmtx, iunit
    write(6,*)
  endif
  ! neig_sub allows to write only neig_sub columns of the actual matrix
  nmtx_col = nmtx
  IF(PRESENT(neig_sub)) nmtx_col = neig_sub
  if (peinf%inode .eq. 0) then
    call progress_init(prog_info, 'writing matrix', 'column', nmtx_col)
    do jj = 1, nmtx_col
      call progress_step(prog_info, jj)
      write(iunit) (matrix(ii, jj), ii = 1, nmtx)
    enddo
    call progress_free(prog_info)
    !XXXX
    do jj = nmtx_col + 1, nmtx
      write(iunit)
    end do
    !XXXX
  endif
 
  return
end subroutine write_matrix_d_sub
!=================================================================================
subroutine write_matrix_f(scal,nfreq,retarded,nmtx,iunit,nfreq_group,advanced)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq
  complex(DPC), intent(in) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iunit
  integer, intent(in) :: nfreq_group
  complex(DPC), optional, intent(in) :: advanced(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  integer :: ii, jj, ifreq
  type(progress_info) :: prog_info !< a user-friendly progress report
  logical :: has_advanced
 
  if (peinf%verb_debug .and. peinf%inode==0) then
    write(6,*) 'Writing matrix: ', nfreq, nmtx, iunit
    write(6,*)
  endif
  has_advanced = present(advanced)
  if(peinf%inode .eq. 0) then
    call progress_init(prog_info, 'writing matrix', 'column', nmtx)
    do jj = 1, nmtx
      call progress_step(prog_info, jj)
      do ii = 1, nmtx
        write(iunit) (retarded(ii, jj, ifreq), ifreq= 1, nfreq)
      enddo
    enddo
    call progress_free(prog_info)
  endif
 
  return
end subroutine write_matrix_f
end module write_matrix_m
