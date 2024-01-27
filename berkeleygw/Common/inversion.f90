!==============================================================================
!
! Routines:
!
! (1) Xinvert_with_scalapack_d() Originally by JRD Last Modified 02/2015 (FHJ)
!
! This routine inverts a matrix which is already distributed in block
! cyclic form with ScaLAPACK.
!
! (2) Xinvert_serial() Originally by JRD Last Modified 02/2015 (FHJ)
!
! Inverts a matrix using LAPACK.
!
!==============================================================================
module inversion_m
  use global_m
  use lapack_m
  use scalapack_m
  implicit none
  private
  public :: &
    dinvert_serial, &
    zinvert_serial
contains
!overrules flavor.mk
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
!===========================================================================
!
! Included from inversion.F90
!
!============================================================================
!---------------------- Use scaLAPACK For Inversion -----------------------------------
!------------------------------------------------------------
subroutine dinvert_serial(nmtx, matrix)
  integer, intent(in) :: nmtx
  real(DP), intent(inout) :: matrix(nmtx,nmtx)
  integer :: ii, info, lwork, ipiv(nmtx)
  real(DP), allocatable :: work(:)
 
  ! FHJ: LU factorization of the matrix
  call dgetrf(nmtx, nmtx, matrix, nmtx, ipiv, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getrf'
    call die('?getrf failed')
  endif
  ! FHJ: tringular inversion of LU decomposition
  allocate(work (10))
  call dgetri(nmtx, matrix, nmtx, ipiv, work, -1, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
    call die('?getri failed for query mode')
  endif
  lwork = max(1,int(work(1)))
  if(allocated(work))then;deallocate(work);endif
  allocate(work (lwork))
  call dgetri(nmtx, matrix, nmtx, ipiv, work, lwork, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
    call die('?getri failed')
  endif
  if(allocated(work))then;deallocate(work);endif
 
  return
end subroutine dinvert_serial
! use between inclusions of f_defs.h in template modules
! list here everything defined differently by flavor in f_defs.h
! these undefs prevent lots of warnings from cpp
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
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
!===========================================================================
!
! Included from inversion.F90
!
!============================================================================
!---------------------- Use scaLAPACK For Inversion -----------------------------------
!------------------------------------------------------------
subroutine zinvert_serial(nmtx, matrix)
  integer, intent(in) :: nmtx
  complex(DPC), intent(inout) :: matrix(nmtx,nmtx)
  integer :: ii, info, lwork, ipiv(nmtx)
  complex(DPC), allocatable :: work(:)
 
  ! FHJ: LU factorization of the matrix
  call zgetrf(nmtx, nmtx, matrix, nmtx, ipiv, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getrf'
    call die('?getrf failed')
  endif
  ! FHJ: tringular inversion of LU decomposition
  allocate(work (10))
  call zgetri(nmtx, matrix, nmtx, ipiv, work, -1, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
    call die('?getri failed for query mode')
  endif
  lwork = max(1,int(work(1)))
  if(allocated(work))then;deallocate(work);endif
  allocate(work (lwork))
  call zgetri(nmtx, matrix, nmtx, ipiv, work, lwork, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
    call die('?getri failed')
  endif
  if(allocated(work))then;deallocate(work);endif
 
  return
end subroutine zinvert_serial
end module inversion_m
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
