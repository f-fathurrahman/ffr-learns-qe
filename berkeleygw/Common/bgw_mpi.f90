! File bgw_mpi.f90 automatically created from bgw_mpi.f90p by mako_preprocess.py.
! Do not edit the resulting file (bgw_mpi.f90) directly!
!>=========================================================================
!!
!! Module:
!!
!! bgw_mpi_m Originally by FHJ Last Modified 03/2018 (FHJ)
!!
!! Wrappers around MPI functions. All functions work both with
!! MPI and serial runs. Right now, only MPI_Bcast is supported.
!! TODO: implement wrappers for MPI_Send, MPI_Recv, gather, etc.
!!
!! This module uses python Mako template to generate Fortran code.
!! Make sure to run the `mako_preprocess.py` if you edit the source
!! .f90p file, and don`t edit the generated .f90 directly.
!!
!!=========================================================================
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
module bgw_mpi_m
  use, intrinsic :: iso_c_binding
  use global_m
  implicit none
  interface bgw_bcast
     module procedure &
       bgw_bcast_str, &
       bgw_bcast_int_0, &
       bgw_bcast_int_1, &
       bgw_bcast_int_2, &
       bgw_bcast_int_3, &
       bgw_bcast_int_4, &
       bgw_bcast_int_5, &
       bgw_bcast_int_6, &
       bgw_bcast_int_7, &
       bgw_bcast_real_0, &
       bgw_bcast_real_1, &
       bgw_bcast_real_2, &
       bgw_bcast_real_3, &
       bgw_bcast_real_4, &
       bgw_bcast_real_5, &
       bgw_bcast_real_6, &
       bgw_bcast_real_7, &
       bgw_bcast_complex_0, &
       bgw_bcast_complex_1, &
       bgw_bcast_complex_2, &
       bgw_bcast_complex_3, &
       bgw_bcast_complex_4, &
       bgw_bcast_complex_5, &
       bgw_bcast_complex_6, &
       bgw_bcast_complex_7, &
       bgw_bcast_logical_0, &
       bgw_bcast_logical_1, &
       bgw_bcast_logical_2, &
       bgw_bcast_logical_3, &
       bgw_bcast_logical_4, &
       bgw_bcast_logical_5, &
       bgw_bcast_logical_6, &
       bgw_bcast_logical_7
  end interface
  interface bgw_gather
     module procedure &
       bgw_gather_int_0, &
       bgw_gather_int_1, &
       bgw_gather_int_2, &
       bgw_gather_int_3, &
       bgw_gather_int_4, &
       bgw_gather_int_5, &
       bgw_gather_int_6, &
       bgw_gather_int_7, &
       bgw_gather_real_0, &
       bgw_gather_real_1, &
       bgw_gather_real_2, &
       bgw_gather_real_3, &
       bgw_gather_real_4, &
       bgw_gather_real_5, &
       bgw_gather_real_6, &
       bgw_gather_real_7, &
       bgw_gather_complex_0, &
       bgw_gather_complex_1, &
       bgw_gather_complex_2, &
       bgw_gather_complex_3, &
       bgw_gather_complex_4, &
       bgw_gather_complex_5, &
       bgw_gather_complex_6, &
       bgw_gather_complex_7, &
       bgw_gather_logical_0, &
       bgw_gather_logical_1, &
       bgw_gather_logical_2, &
       bgw_gather_logical_3, &
       bgw_gather_logical_4, &
       bgw_gather_logical_5, &
       bgw_gather_logical_6, &
       bgw_gather_logical_7
  end interface
  interface bgw_allgather
     module procedure &
       bgw_allgather_int_0, &
       bgw_allgather_int_1, &
       bgw_allgather_int_2, &
       bgw_allgather_int_3, &
       bgw_allgather_int_4, &
       bgw_allgather_int_5, &
       bgw_allgather_int_6, &
       bgw_allgather_int_7, &
       bgw_allgather_real_0, &
       bgw_allgather_real_1, &
       bgw_allgather_real_2, &
       bgw_allgather_real_3, &
       bgw_allgather_real_4, &
       bgw_allgather_real_5, &
       bgw_allgather_real_6, &
       bgw_allgather_real_7, &
       bgw_allgather_complex_0, &
       bgw_allgather_complex_1, &
       bgw_allgather_complex_2, &
       bgw_allgather_complex_3, &
       bgw_allgather_complex_4, &
       bgw_allgather_complex_5, &
       bgw_allgather_complex_6, &
       bgw_allgather_complex_7, &
       bgw_allgather_logical_0, &
       bgw_allgather_logical_1, &
       bgw_allgather_logical_2, &
       bgw_allgather_logical_3, &
       bgw_allgather_logical_4, &
       bgw_allgather_logical_5, &
       bgw_allgather_logical_6, &
       bgw_allgather_logical_7
  end interface
  interface bgw_allgatherv
     module procedure &
       bgw_allgatherv_int_0, &
       bgw_allgatherv_int_1, &
       bgw_allgatherv_int_2, &
       bgw_allgatherv_int_3, &
       bgw_allgatherv_int_4, &
       bgw_allgatherv_int_5, &
       bgw_allgatherv_int_6, &
       bgw_allgatherv_int_7, &
       bgw_allgatherv_real_0, &
       bgw_allgatherv_real_1, &
       bgw_allgatherv_real_2, &
       bgw_allgatherv_real_3, &
       bgw_allgatherv_real_4, &
       bgw_allgatherv_real_5, &
       bgw_allgatherv_real_6, &
       bgw_allgatherv_real_7, &
       bgw_allgatherv_complex_0, &
       bgw_allgatherv_complex_1, &
       bgw_allgatherv_complex_2, &
       bgw_allgatherv_complex_3, &
       bgw_allgatherv_complex_4, &
       bgw_allgatherv_complex_5, &
       bgw_allgatherv_complex_6, &
       bgw_allgatherv_complex_7, &
       bgw_allgatherv_logical_0, &
       bgw_allgatherv_logical_1, &
       bgw_allgatherv_logical_2, &
       bgw_allgatherv_logical_3, &
       bgw_allgatherv_logical_4, &
       bgw_allgatherv_logical_5, &
       bgw_allgatherv_logical_6, &
       bgw_allgatherv_logical_7
  end interface
  public :: bgw_bcast, bgw_gather, bgw_allgather, bgw_allgatherv, &
    bgw_barrier, bgw_comm_size
contains
!==============================================================================
! Auxiliary routines
!==============================================================================
subroutine bgw_mpi_error(bgw_routine, mpi_routine, ierr, comm, root)
  character(len=*), intent(in) :: bgw_routine
  character(len=*), intent(in) :: mpi_routine
  integer, intent(in) :: ierr
  integer, intent(in), optional :: comm
  integer, intent(in), optional :: root
  character(len=256) :: error_str
  integer :: comm_, root_
 
 
end subroutine bgw_mpi_error
!==============================================================================
! High-level routines
!==============================================================================
!------------------------------------------------------------------------------
! Comm_size
!------------------------------------------------------------------------------
integer function bgw_comm_size(comm)
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  bgw_comm_size = 1
 
end function bgw_comm_size
!------------------------------------------------------------------------------
! Comm_rank
!------------------------------------------------------------------------------
integer function bgw_comm_rank(comm)
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  bgw_comm_rank = 0
 
end function bgw_comm_rank
!------------------------------------------------------------------------------
! Barrier
!------------------------------------------------------------------------------
subroutine bgw_barrier(comm)
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
 
end subroutine bgw_barrier
!------------------------------------------------------------------------------
! Bcast
!------------------------------------------------------------------------------
subroutine bgw_bcast_int_0(x, root, comm)
  integer, intent(inout), target :: x
  integer, intent(in), optional :: root, comm
  integer, pointer :: x_1d(:)
 
  call c_f_pointer(c_loc(x), x_1d, [1])
  call bgw_bcast_int_lowlevel(x_1d, 1, root=root, comm=comm)
 
end subroutine bgw_bcast_int_0
subroutine bgw_bcast_int_1(x, root, comm)
  integer, intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_int_1
subroutine bgw_bcast_int_2(x, root, comm)
  integer, intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_int_2
subroutine bgw_bcast_int_3(x, root, comm)
  integer, intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_int_3
subroutine bgw_bcast_int_4(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_int_4
subroutine bgw_bcast_int_5(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_int_5
subroutine bgw_bcast_int_6(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_int_6
subroutine bgw_bcast_int_7(x, root, comm)
  integer, intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_int_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_int_7
subroutine bgw_bcast_real_0(x, root, comm)
  real(DP), intent(inout), target :: x
  integer, intent(in), optional :: root, comm
  real(DP), pointer :: x_1d(:)
 
  call c_f_pointer(c_loc(x), x_1d, [1])
  call bgw_bcast_real_lowlevel(x_1d, 1, root=root, comm=comm)
 
end subroutine bgw_bcast_real_0
subroutine bgw_bcast_real_1(x, root, comm)
  real(DP), intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_real_1
subroutine bgw_bcast_real_2(x, root, comm)
  real(DP), intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_real_2
subroutine bgw_bcast_real_3(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_real_3
subroutine bgw_bcast_real_4(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_real_4
subroutine bgw_bcast_real_5(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_real_5
subroutine bgw_bcast_real_6(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_real_6
subroutine bgw_bcast_real_7(x, root, comm)
  real(DP), intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_real_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_real_7
subroutine bgw_bcast_complex_0(x, root, comm)
  complex(DPC), intent(inout), target :: x
  integer, intent(in), optional :: root, comm
  complex(DPC), pointer :: x_1d(:)
 
  call c_f_pointer(c_loc(x), x_1d, [1])
  call bgw_bcast_complex_lowlevel(x_1d, 1, root=root, comm=comm)
 
end subroutine bgw_bcast_complex_0
subroutine bgw_bcast_complex_1(x, root, comm)
  complex(DPC), intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_complex_1
subroutine bgw_bcast_complex_2(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_complex_2
subroutine bgw_bcast_complex_3(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_complex_3
subroutine bgw_bcast_complex_4(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_complex_4
subroutine bgw_bcast_complex_5(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_complex_5
subroutine bgw_bcast_complex_6(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_complex_6
subroutine bgw_bcast_complex_7(x, root, comm)
  complex(DPC), intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_complex_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_complex_7
subroutine bgw_bcast_logical_0(x, root, comm)
  logical, intent(inout), target :: x
  integer, intent(in), optional :: root, comm
  logical, pointer :: x_1d(:)
 
  call c_f_pointer(c_loc(x), x_1d, [1])
  call bgw_bcast_logical_lowlevel(x_1d, 1, root=root, comm=comm)
 
end subroutine bgw_bcast_logical_0
subroutine bgw_bcast_logical_1(x, root, comm)
  logical, intent(inout) :: x(:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_logical_1
subroutine bgw_bcast_logical_2(x, root, comm)
  logical, intent(inout) :: x(:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_logical_2
subroutine bgw_bcast_logical_3(x, root, comm)
  logical, intent(inout) :: x(:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_logical_3
subroutine bgw_bcast_logical_4(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_logical_4
subroutine bgw_bcast_logical_5(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_logical_5
subroutine bgw_bcast_logical_6(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_logical_6
subroutine bgw_bcast_logical_7(x, root, comm)
  logical, intent(inout) :: x(:,:,:,:,:,:,:)
  integer, intent(in), optional :: root, comm
 
  call bgw_bcast_logical_lowlevel(x, size(x), root=root, comm=comm)
 
end subroutine bgw_bcast_logical_7
!------------------------------------------------------------------------------
! Gather
!------------------------------------------------------------------------------
subroutine bgw_gather_int_0(x_in, x_out, root, comm)
  integer, intent(in), target :: x_in
  integer, intent(out), target :: x_out
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
  integer, pointer :: x_in_1d(:)
  integer, pointer :: x_out_1d(:)
 
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_gather_int_lowlevel(x_in_1d, x_out_1d, 1, root, comm=comm)
 
end subroutine bgw_gather_int_0
subroutine bgw_gather_int_1(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:)
  integer, intent(out) :: x_out(:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_int_1
subroutine bgw_gather_int_2(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:)
  integer, intent(out) :: x_out(:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_int_2
subroutine bgw_gather_int_3(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:,:)
  integer, intent(out) :: x_out(:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_int_3
subroutine bgw_gather_int_4(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_int_4
subroutine bgw_gather_int_5(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_int_5
subroutine bgw_gather_int_6(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_int_6
subroutine bgw_gather_int_7(x_in, x_out, root, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_int_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_int_7
subroutine bgw_gather_real_0(x_in, x_out, root, comm)
  real(DP), intent(in), target :: x_in
  real(DP), intent(out), target :: x_out
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
  real(DP), pointer :: x_in_1d(:)
  real(DP), pointer :: x_out_1d(:)
 
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_gather_real_lowlevel(x_in_1d, x_out_1d, 1, root, comm=comm)
 
end subroutine bgw_gather_real_0
subroutine bgw_gather_real_1(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:)
  real(DP), intent(out) :: x_out(:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_real_1
subroutine bgw_gather_real_2(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:)
  real(DP), intent(out) :: x_out(:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_real_2
subroutine bgw_gather_real_3(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:,:)
  real(DP), intent(out) :: x_out(:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_real_3
subroutine bgw_gather_real_4(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_real_4
subroutine bgw_gather_real_5(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_real_5
subroutine bgw_gather_real_6(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_real_6
subroutine bgw_gather_real_7(x_in, x_out, root, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_real_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_real_7
subroutine bgw_gather_complex_0(x_in, x_out, root, comm)
  complex(DPC), intent(in), target :: x_in
  complex(DPC), intent(out), target :: x_out
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
  complex(DPC), pointer :: x_in_1d(:)
  complex(DPC), pointer :: x_out_1d(:)
 
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_gather_complex_lowlevel(x_in_1d, x_out_1d, 1, root, comm=comm)
 
end subroutine bgw_gather_complex_0
subroutine bgw_gather_complex_1(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:)
  complex(DPC), intent(out) :: x_out(:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_complex_1
subroutine bgw_gather_complex_2(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:)
  complex(DPC), intent(out) :: x_out(:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_complex_2
subroutine bgw_gather_complex_3(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_complex_3
subroutine bgw_gather_complex_4(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_complex_4
subroutine bgw_gather_complex_5(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_complex_5
subroutine bgw_gather_complex_6(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_complex_6
subroutine bgw_gather_complex_7(x_in, x_out, root, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_complex_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_complex_7
subroutine bgw_gather_logical_0(x_in, x_out, root, comm)
  logical, intent(in), target :: x_in
  logical, intent(out), target :: x_out
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
  logical, pointer :: x_in_1d(:)
  logical, pointer :: x_out_1d(:)
 
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_gather_logical_lowlevel(x_in_1d, x_out_1d, 1, root, comm=comm)
 
end subroutine bgw_gather_logical_0
subroutine bgw_gather_logical_1(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:)
  logical, intent(out) :: x_out(:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_logical_1
subroutine bgw_gather_logical_2(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:)
  logical, intent(out) :: x_out(:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_logical_2
subroutine bgw_gather_logical_3(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:,:)
  logical, intent(out) :: x_out(:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_logical_3
subroutine bgw_gather_logical_4(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_logical_4
subroutine bgw_gather_logical_5(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_logical_5
subroutine bgw_gather_logical_6(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_logical_6
subroutine bgw_gather_logical_7(x_in, x_out, root, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
 
  call bgw_gather_logical_lowlevel(x_in, x_out, size(x_in), root, comm=comm)
 
end subroutine bgw_gather_logical_7
!------------------------------------------------------------------------------
! Allgather
!------------------------------------------------------------------------------
subroutine bgw_allgather_int_0(x_in, x_out, comm)
  integer, intent(in), target :: x_in
  integer, intent(out), target :: x_out
  integer, intent(in), optional :: comm
  integer, pointer :: x_in_1d(:)
  integer, pointer :: x_out_1d(:)
 
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgather_int_lowlevel(x_in_1d, x_out_1d, 1, comm=comm)
 
end subroutine bgw_allgather_int_0
subroutine bgw_allgather_int_1(x_in, x_out, comm)
  integer, intent(in) :: x_in(:)
  integer, intent(out) :: x_out(:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_int_1
subroutine bgw_allgather_int_2(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:)
  integer, intent(out) :: x_out(:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_int_2
subroutine bgw_allgather_int_3(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:,:)
  integer, intent(out) :: x_out(:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_int_3
subroutine bgw_allgather_int_4(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_int_4
subroutine bgw_allgather_int_5(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_int_5
subroutine bgw_allgather_int_6(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_int_6
subroutine bgw_allgather_int_7(x_in, x_out, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_int_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_int_7
subroutine bgw_allgather_real_0(x_in, x_out, comm)
  real(DP), intent(in), target :: x_in
  real(DP), intent(out), target :: x_out
  integer, intent(in), optional :: comm
  real(DP), pointer :: x_in_1d(:)
  real(DP), pointer :: x_out_1d(:)
 
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgather_real_lowlevel(x_in_1d, x_out_1d, 1, comm=comm)
 
end subroutine bgw_allgather_real_0
subroutine bgw_allgather_real_1(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:)
  real(DP), intent(out) :: x_out(:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_real_1
subroutine bgw_allgather_real_2(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:)
  real(DP), intent(out) :: x_out(:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_real_2
subroutine bgw_allgather_real_3(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:,:)
  real(DP), intent(out) :: x_out(:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_real_3
subroutine bgw_allgather_real_4(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_real_4
subroutine bgw_allgather_real_5(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_real_5
subroutine bgw_allgather_real_6(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_real_6
subroutine bgw_allgather_real_7(x_in, x_out, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_real_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_real_7
subroutine bgw_allgather_complex_0(x_in, x_out, comm)
  complex(DPC), intent(in), target :: x_in
  complex(DPC), intent(out), target :: x_out
  integer, intent(in), optional :: comm
  complex(DPC), pointer :: x_in_1d(:)
  complex(DPC), pointer :: x_out_1d(:)
 
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgather_complex_lowlevel(x_in_1d, x_out_1d, 1, comm=comm)
 
end subroutine bgw_allgather_complex_0
subroutine bgw_allgather_complex_1(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:)
  complex(DPC), intent(out) :: x_out(:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_complex_1
subroutine bgw_allgather_complex_2(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:)
  complex(DPC), intent(out) :: x_out(:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_complex_2
subroutine bgw_allgather_complex_3(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_complex_3
subroutine bgw_allgather_complex_4(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_complex_4
subroutine bgw_allgather_complex_5(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_complex_5
subroutine bgw_allgather_complex_6(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_complex_6
subroutine bgw_allgather_complex_7(x_in, x_out, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_complex_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_complex_7
subroutine bgw_allgather_logical_0(x_in, x_out, comm)
  logical, intent(in), target :: x_in
  logical, intent(out), target :: x_out
  integer, intent(in), optional :: comm
  logical, pointer :: x_in_1d(:)
  logical, pointer :: x_out_1d(:)
 
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgather_logical_lowlevel(x_in_1d, x_out_1d, 1, comm=comm)
 
end subroutine bgw_allgather_logical_0
subroutine bgw_allgather_logical_1(x_in, x_out, comm)
  logical, intent(in) :: x_in(:)
  logical, intent(out) :: x_out(:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_logical_1
subroutine bgw_allgather_logical_2(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:)
  logical, intent(out) :: x_out(:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_logical_2
subroutine bgw_allgather_logical_3(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:,:)
  logical, intent(out) :: x_out(:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_logical_3
subroutine bgw_allgather_logical_4(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_logical_4
subroutine bgw_allgather_logical_5(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_logical_5
subroutine bgw_allgather_logical_6(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_logical_6
subroutine bgw_allgather_logical_7(x_in, x_out, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, intent(in), optional :: comm
 
  call bgw_allgather_logical_lowlevel(x_in, x_out, size(x_in), comm=comm)
 
end subroutine bgw_allgather_logical_7
!------------------------------------------------------------------------------
! Allgatherv
!------------------------------------------------------------------------------
subroutine bgw_allgatherv_int_0(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in), target :: x_in
  integer, intent(out), target :: x_out
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: x_in_1d(:)
  integer, pointer :: x_out_1d(:)
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    siz = 1
  else
    siz = 1
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgatherv_int_lowlevel(x_in_1d, x_out_1d, 1, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_int_0
subroutine bgw_allgatherv_int_1(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:)
  integer, intent(out) :: x_out(:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_int_1
subroutine bgw_allgatherv_int_2(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:)
  integer, intent(out) :: x_out(:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_int_2
subroutine bgw_allgatherv_int_3(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:,:)
  integer, intent(out) :: x_out(:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_int_3
subroutine bgw_allgatherv_int_4(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_int_4
subroutine bgw_allgatherv_int_5(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_int_5
subroutine bgw_allgatherv_int_6(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_int_6
subroutine bgw_allgatherv_int_7(x_in, x_out, recvcounts, displs, comm)
  integer, intent(in) :: x_in(:,:,:,:,:,:,:)
  integer, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_int_7
subroutine bgw_allgatherv_real_0(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in), target :: x_in
  real(DP), intent(out), target :: x_out
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  real(DP), pointer :: x_in_1d(:)
  real(DP), pointer :: x_out_1d(:)
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    siz = 1
  else
    siz = 1
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgatherv_real_lowlevel(x_in_1d, x_out_1d, 1, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_real_0
subroutine bgw_allgatherv_real_1(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:)
  real(DP), intent(out) :: x_out(:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_real_1
subroutine bgw_allgatherv_real_2(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:)
  real(DP), intent(out) :: x_out(:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_real_2
subroutine bgw_allgatherv_real_3(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:,:)
  real(DP), intent(out) :: x_out(:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_real_3
subroutine bgw_allgatherv_real_4(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_real_4
subroutine bgw_allgatherv_real_5(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_real_5
subroutine bgw_allgatherv_real_6(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_real_6
subroutine bgw_allgatherv_real_7(x_in, x_out, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(:,:,:,:,:,:,:)
  real(DP), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_real_7
subroutine bgw_allgatherv_complex_0(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in), target :: x_in
  complex(DPC), intent(out), target :: x_out
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  complex(DPC), pointer :: x_in_1d(:)
  complex(DPC), pointer :: x_out_1d(:)
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    siz = 1
  else
    siz = 1
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgatherv_complex_lowlevel(x_in_1d, x_out_1d, 1, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_complex_0
subroutine bgw_allgatherv_complex_1(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:)
  complex(DPC), intent(out) :: x_out(:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_complex_1
subroutine bgw_allgatherv_complex_2(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:)
  complex(DPC), intent(out) :: x_out(:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_complex_2
subroutine bgw_allgatherv_complex_3(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_complex_3
subroutine bgw_allgatherv_complex_4(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_complex_4
subroutine bgw_allgatherv_complex_5(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_complex_5
subroutine bgw_allgatherv_complex_6(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_complex_6
subroutine bgw_allgatherv_complex_7(x_in, x_out, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(:,:,:,:,:,:,:)
  complex(DPC), intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_complex_7
subroutine bgw_allgatherv_logical_0(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in), target :: x_in
  logical, intent(out), target :: x_out
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  logical, pointer :: x_in_1d(:)
  logical, pointer :: x_out_1d(:)
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    siz = 1
  else
    siz = 1
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call c_f_pointer(c_loc(x_in), x_in_1d, [1])
  call c_f_pointer(c_loc(x_out), x_out_1d, [1])
  call bgw_allgatherv_logical_lowlevel(x_in_1d, x_out_1d, 1, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_logical_0
subroutine bgw_allgatherv_logical_1(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:)
  logical, intent(out) :: x_out(:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_logical_1
subroutine bgw_allgatherv_logical_2(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:)
  logical, intent(out) :: x_out(:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_logical_2
subroutine bgw_allgatherv_logical_3(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:,:)
  logical, intent(out) :: x_out(:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_logical_3
subroutine bgw_allgatherv_logical_4(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_logical_4
subroutine bgw_allgatherv_logical_5(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_logical_5
subroutine bgw_allgatherv_logical_6(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_logical_6
subroutine bgw_allgatherv_logical_7(x_in, x_out, recvcounts, displs, comm)
  logical, intent(in) :: x_in(:,:,:,:,:,:,:)
  logical, intent(out) :: x_out(:,:,:,:,:,:,:)
  integer, optional, target, intent(in) :: recvcounts(:)
  integer, optional, target, intent(in) :: displs(:)
  integer, intent(in), optional :: comm
  integer, pointer :: recvcounts_(:)
  integer, pointer :: displs_(:)
  integer :: npes, siz, ipe
 
  if ((.not.present(recvcounts)) .or. (.not.present(displs))) then
    npes = peinf%npes
    if (present(comm)) npes = bgw_comm_size(comm)
  endif
  if (present(recvcounts)) then
    recvcounts_ => recvcounts
    if (present(comm)) then
      ipe = bgw_comm_rank(comm)
    else
      ipe = peinf%inode
    endif
    siz = recvcounts(ipe+1)
  else
    siz = size(x_in)
    allocate(recvcounts_ (npes))
    recvcounts_(:) = 0
    call bgw_allgather([siz], recvcounts_)
  endif
  if (present(displs)) then
    displs_ => displs
  else
    allocate(displs_ (npes))
    displs_(1) = 0
    do ipe = 2, npes
      displs_(ipe) = displs_(ipe-1) + recvcounts_(ipe-1)
    enddo
  endif
  call bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts_, displs_, comm)
  if (.not.present(recvcounts)) then
    if(associated(recvcounts_))then;deallocate(recvcounts_);nullify(recvcounts_);endif
  endif
  if (.not.present(displs)) then
    if(associated(displs_))then;deallocate(displs_);nullify(displs_);endif
  endif
 
end subroutine bgw_allgatherv_logical_7
!==============================================================================
! Low-level implementations
!==============================================================================
!------------------------------------------------------------------------------
! Bcast: str
!------------------------------------------------------------------------------
subroutine bgw_bcast_str(x, root, comm)
  character(len=*), intent(inout) :: x
  integer, intent(in), optional :: root, comm
  integer :: root_, comm_, ierr
 
 
end subroutine bgw_bcast_str
!------------------------------------------------------------------------------
! Bcast: int, real, complex, logical
!------------------------------------------------------------------------------
subroutine bgw_bcast_int_lowlevel(x, siz, root, comm)
  integer, intent(inout) :: x(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: root, comm
  integer :: root_, comm_, ierr
 
 
end subroutine bgw_bcast_int_lowlevel
subroutine bgw_bcast_real_lowlevel(x, siz, root, comm)
  real(DP), intent(inout) :: x(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: root, comm
  integer :: root_, comm_, ierr
 
 
end subroutine bgw_bcast_real_lowlevel
subroutine bgw_bcast_complex_lowlevel(x, siz, root, comm)
  complex(DPC), intent(inout) :: x(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: root, comm
  integer :: root_, comm_, ierr
 
 
end subroutine bgw_bcast_complex_lowlevel
subroutine bgw_bcast_logical_lowlevel(x, siz, root, comm)
  logical, intent(inout) :: x(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: root, comm
  integer :: root_, comm_, ierr
 
 
end subroutine bgw_bcast_logical_lowlevel
!------------------------------------------------------------------------------
! Gather: int, real, complex, logical
!------------------------------------------------------------------------------
subroutine bgw_gather_int_lowlevel(x_in, x_out, siz, root, comm)
  integer, intent(in) :: x_in(*)
  integer, intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_gather_int_lowlevel
subroutine bgw_gather_real_lowlevel(x_in, x_out, siz, root, comm)
  real(DP), intent(in) :: x_in(*)
  real(DP), intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_gather_real_lowlevel
subroutine bgw_gather_complex_lowlevel(x_in, x_out, siz, root, comm)
  complex(DPC), intent(in) :: x_in(*)
  complex(DPC), intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_gather_complex_lowlevel
subroutine bgw_gather_logical_lowlevel(x_in, x_out, siz, root, comm)
  logical, intent(in) :: x_in(*)
  logical, intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in) :: root
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_gather_logical_lowlevel
!------------------------------------------------------------------------------
! Allgather: int, real, complex, logical
!------------------------------------------------------------------------------
subroutine bgw_allgather_int_lowlevel(x_in, x_out, siz, comm)
  integer, intent(in) :: x_in(*)
  integer, intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_allgather_int_lowlevel
subroutine bgw_allgather_real_lowlevel(x_in, x_out, siz, comm)
  real(DP), intent(in) :: x_in(*)
  real(DP), intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_allgather_real_lowlevel
subroutine bgw_allgather_complex_lowlevel(x_in, x_out, siz, comm)
  complex(DPC), intent(in) :: x_in(*)
  complex(DPC), intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_allgather_complex_lowlevel
subroutine bgw_allgather_logical_lowlevel(x_in, x_out, siz, comm)
  logical, intent(in) :: x_in(*)
  logical, intent(out) :: x_out(*)
  integer, intent(in) :: siz
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_allgather_logical_lowlevel
!------------------------------------------------------------------------------
! Allgatherv: int, real, complex, logical
!------------------------------------------------------------------------------
subroutine bgw_allgatherv_int_lowlevel(x_in, x_out, siz, recvcounts, displs, comm)
  integer, intent(in) :: x_in(*)
  integer, intent(out) :: x_out(*)
  integer, intent(in) :: siz, recvcounts(*), displs(*)
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_allgatherv_int_lowlevel
subroutine bgw_allgatherv_real_lowlevel(x_in, x_out, siz, recvcounts, displs, comm)
  real(DP), intent(in) :: x_in(*)
  real(DP), intent(out) :: x_out(*)
  integer, intent(in) :: siz, recvcounts(*), displs(*)
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_allgatherv_real_lowlevel
subroutine bgw_allgatherv_complex_lowlevel(x_in, x_out, siz, recvcounts, displs, comm)
  complex(DPC), intent(in) :: x_in(*)
  complex(DPC), intent(out) :: x_out(*)
  integer, intent(in) :: siz, recvcounts(*), displs(*)
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_allgatherv_complex_lowlevel
subroutine bgw_allgatherv_logical_lowlevel(x_in, x_out, siz, recvcounts, displs, comm)
  logical, intent(in) :: x_in(*)
  logical, intent(out) :: x_out(*)
  integer, intent(in) :: siz, recvcounts(*), displs(*)
  integer, intent(in), optional :: comm
  integer :: comm_, ierr
 
  if (peinf%npes>1) then
  else
    x_out(1:siz) = x_in(1:siz)
  endif
 
end subroutine bgw_allgatherv_logical_lowlevel
end module bgw_mpi_m
