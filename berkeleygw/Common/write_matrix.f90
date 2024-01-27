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
module write_matrix_m
  use, intrinsic :: iso_c_binding
  use global_m
  use hdf5_io_m
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
