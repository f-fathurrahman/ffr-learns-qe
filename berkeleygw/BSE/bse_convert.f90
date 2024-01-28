!============================================================================
!
! Module bse_convert_m
!
!============================================================================
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
module bse_convert_m
  use global_m
  use kernel_io_m
  implicit none
  private
  public :: &
    bsemat_binasc, &
    dtmat_binasc, &
    vmtxel_binasc, &
    eps2_moments_binasc, &
    eigenvectors_binasc, &
    bsemat_ascbin, &
    dtmat_ascbin, &
    vmtxel_ascbin, &
    eps2_moments_ascbin, &
    eigenvectors_ascbin
contains
! this is a variadic macro
!============================================================================
!
! included from bse_convert.F90
!
!============================================================================
subroutine bsemat_binasc(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  type(kernel_header_t) :: kernel
  real(DP), allocatable :: bsemat(:,:,:,:,:)
  integer :: ik, ikp, imat, ibb, ic, iv
 
  call read_kernel_header(iunit, kernel)
  call write_kernel_header(ounit, kernel)
  allocate(bsemat (kernel%nk, kernel%n2b, kernel%n1b, kernel%ns, kernel%ns))
  do ik = 1, kernel%nk
    do imat = 1, kernel%nmat
      do ibb = 1, kernel%n1b*kernel%n2b
        read(iunit) ikp,ic,iv,bsemat(:,:,:,:,:); write(ounit, *) ikp,ic,iv,bsemat(:,:,:,:,:)
      enddo
    end do
  end do
 
end subroutine bsemat_binasc
!============================================================================
! compare read and write in intwfn.f90
subroutine dtmat_binasc(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  integer :: nk, nc, nv, nkf, ncf, nvf, ns, nmat, &
    ik, ic, jk, jc, is, ii, iv, jv, npts, ndims
  logical :: per(3)
  real(DP) :: kk(1:3)
  real(DP) :: dcc, dvv
 
  read(iunit) ndims, per(1:3), npts, nk; write(ounit, *) ndims, per(1:3), npts, nk
  read(iunit) nk,nc,nv,nkf,ncf,nvf,ns; write(ounit, *) nk,nc,nv,nkf,ncf,nvf,ns
  do ik=1, nk
    read(iunit) kk(1:3); write(ounit, *) kk(1:3)
  enddo
  ! cc
  nmat=nkf*ncf*nc*ns
  do ii=1,nmat
    read(iunit) ik,ic,jk,jc,is,dcc; write(ounit, *) ik,ic,jk,jc,is,dcc
  enddo
  ! vv
  nmat=nkf*nvf*nv*ns
  do ii=1,nmat
    read(iunit) ik,iv,jk,jv,is,dvv; write(ounit, *) ik,iv,jk,jv,is,dvv
  enddo
 
  return
end subroutine dtmat_binasc
!============================================================================
subroutine vmtxel_binasc(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  integer :: nkf, ncf, nvf, ns, ic, nmat, ii
  real(DP), allocatable :: s1(:)
 
  read(iunit) nkf,ncf,nvf,ns,ic; write(ounit, *) nkf,ncf,nvf,ns,ic
  nmat=nkf*ncf*nvf*ns
  allocate(s1 (nmat))
  read(iunit) (s1(ii),ii=1,nmat); write(ounit, *) (s1(ii),ii=1,nmat)
  if(allocated(s1))then;deallocate(s1);endif
 
  return
end subroutine vmtxel_binasc
!============================================================================
subroutine eps2_moments_binasc(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  integer :: nn, nmat, ii
  real(DP) :: tmp(2)
  real(DP), allocatable :: s1(:)
  real(DP), allocatable :: array(:)
 
  read(iunit) nn,tmp(1:2),nmat,ii; write(ounit, *) nn,tmp(1:2),nmat,ii
  allocate(array (nn))
  ! an
  read(iunit) (array(ii),ii=1,nn); write(ounit, *) (array(ii),ii=1,nn)
  ! bn
  read(iunit) (array(ii),ii=1,nn); write(ounit, *) (array(ii),ii=1,nn)
  if(allocated(array))then;deallocate(array);endif
  allocate(s1 (nmat))
  read(iunit) (s1(ii),ii=1,nmat); write(ounit, *) (s1(ii),ii=1,nmat)
  read(iunit) (s1(ii),ii=1,nmat); write(ounit, *) (s1(ii),ii=1,nmat)
  if(allocated(s1))then;deallocate(s1);endif
 
  return
end subroutine eps2_moments_binasc
!============================================================================
subroutine eigenvectors_binasc(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  integer :: ns, nv, nc, nk, nmat, ii, ik, isvck
  real(DP), allocatable :: kg(:,:)
  real(DP) :: energy
  real(DP), allocatable :: Asvck(:)
 
  read(iunit) ns; write(ounit, *) ns
  read(iunit) nv; write(ounit, *) nv
  read(iunit) nc; write(ounit, *) nc
  read(iunit) nk; write(ounit, *) nk
  allocate(kg (3,nk))
  read(iunit) ((kg(ii,ik),ii=1,3),ik=1,nk); write(ounit, *) ((kg(ii,ik),ii=1,3),ik=1,nk)
  if(allocated(kg))then;deallocate(kg);endif
  nmat = ns*nv*nc*nk
  allocate(Asvck (nmat))
  ! FIXME: there can be between 1 and nmat eigenvectors here, actually
  ! The code will crash if it is less than nmat, though what is written
  ! will be correct and complete nonetheless.
  do isvck = 1, nmat
    read(iunit) energy; write(ounit, *) energy
    read(iunit) (Asvck(ik),ik=1,nmat); write(ounit, *) (Asvck(ik),ik=1,nmat)
  enddo
  if(allocated(Asvck))then;deallocate(Asvck);endif
 
  return
end subroutine eigenvectors_binasc
!============================================================================
!
! included from bse_convert.F90
!
!============================================================================
subroutine bsemat_ascbin(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  type(kernel_header_t) :: kernel
  real(DP), allocatable :: bsemat(:,:,:,:,:)
  integer :: ik, ikp, imat, ibb, ic, iv
 
  call read_kernel_header(iunit, kernel)
  call write_kernel_header(ounit, kernel)
  allocate(bsemat (kernel%nk, kernel%n2b, kernel%n1b, kernel%ns, kernel%ns))
  do ik = 1, kernel%nk
    do imat = 1, kernel%nmat
      do ibb = 1, kernel%n1b*kernel%n2b
        read(iunit, *) ikp,ic,iv,bsemat(:,:,:,:,:); write(ounit) ikp,ic,iv,bsemat(:,:,:,:,:)
      enddo
    end do
  end do
 
end subroutine bsemat_ascbin
!============================================================================
! compare read and write in intwfn.f90
subroutine dtmat_ascbin(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  integer :: nk, nc, nv, nkf, ncf, nvf, ns, nmat, &
    ik, ic, jk, jc, is, ii, iv, jv, npts, ndims
  logical :: per(3)
  real(DP) :: kk(1:3)
  real(DP) :: dcc, dvv
 
  read(iunit, *) ndims, per(1:3), npts, nk; write(ounit) ndims, per(1:3), npts, nk
  read(iunit, *) nk,nc,nv,nkf,ncf,nvf,ns; write(ounit) nk,nc,nv,nkf,ncf,nvf,ns
  do ik=1, nk
    read(iunit, *) kk(1:3); write(ounit) kk(1:3)
  enddo
  ! cc
  nmat=nkf*ncf*nc*ns
  do ii=1,nmat
    read(iunit, *) ik,ic,jk,jc,is,dcc; write(ounit) ik,ic,jk,jc,is,dcc
  enddo
  ! vv
  nmat=nkf*nvf*nv*ns
  do ii=1,nmat
    read(iunit, *) ik,iv,jk,jv,is,dvv; write(ounit) ik,iv,jk,jv,is,dvv
  enddo
 
  return
end subroutine dtmat_ascbin
!============================================================================
subroutine vmtxel_ascbin(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  integer :: nkf, ncf, nvf, ns, ic, nmat, ii
  real(DP), allocatable :: s1(:)
 
  read(iunit, *) nkf,ncf,nvf,ns,ic; write(ounit) nkf,ncf,nvf,ns,ic
  nmat=nkf*ncf*nvf*ns
  allocate(s1 (nmat))
  read(iunit, *) (s1(ii),ii=1,nmat); write(ounit) (s1(ii),ii=1,nmat)
  if(allocated(s1))then;deallocate(s1);endif
 
  return
end subroutine vmtxel_ascbin
!============================================================================
subroutine eps2_moments_ascbin(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  integer :: nn, nmat, ii
  real(DP) :: tmp(2)
  real(DP), allocatable :: s1(:)
  real(DP), allocatable :: array(:)
 
  read(iunit, *) nn,tmp(1:2),nmat,ii; write(ounit) nn,tmp(1:2),nmat,ii
  allocate(array (nn))
  ! an
  read(iunit, *) (array(ii),ii=1,nn); write(ounit) (array(ii),ii=1,nn)
  ! bn
  read(iunit, *) (array(ii),ii=1,nn); write(ounit) (array(ii),ii=1,nn)
  if(allocated(array))then;deallocate(array);endif
  allocate(s1 (nmat))
  read(iunit, *) (s1(ii),ii=1,nmat); write(ounit) (s1(ii),ii=1,nmat)
  read(iunit, *) (s1(ii),ii=1,nmat); write(ounit) (s1(ii),ii=1,nmat)
  if(allocated(s1))then;deallocate(s1);endif
 
  return
end subroutine eps2_moments_ascbin
!============================================================================
subroutine eigenvectors_ascbin(iunit, ounit)
  integer, intent(in) :: iunit, ounit
  integer :: ns, nv, nc, nk, nmat, ii, ik, isvck
  real(DP), allocatable :: kg(:,:)
  real(DP) :: energy
  real(DP), allocatable :: Asvck(:)
 
  read(iunit, *) ns; write(ounit) ns
  read(iunit, *) nv; write(ounit) nv
  read(iunit, *) nc; write(ounit) nc
  read(iunit, *) nk; write(ounit) nk
  allocate(kg (3,nk))
  read(iunit, *) ((kg(ii,ik),ii=1,3),ik=1,nk); write(ounit) ((kg(ii,ik),ii=1,3),ik=1,nk)
  if(allocated(kg))then;deallocate(kg);endif
  nmat = ns*nv*nc*nk
  allocate(Asvck (nmat))
  ! FIXME: there can be between 1 and nmat eigenvectors here, actually
  ! The code will crash if it is less than nmat, though what is written
  ! will be correct and complete nonetheless.
  do isvck = 1, nmat
    read(iunit, *) energy; write(ounit) energy
    read(iunit, *) (Asvck(ik),ik=1,nmat); write(ounit) (Asvck(ik),ik=1,nmat)
  enddo
  if(allocated(Asvck))then;deallocate(Asvck);endif
 
  return
end subroutine eigenvectors_ascbin
end module bse_convert_m
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
