!===============================================================================
!
! Routines:
!
! 1. read_cube() Originally By gsm Last Modified 9/3/2010 (gsm)
!
! Reads Gaussian Cube file fnam on unit unum. The result is placed
! into real (if ip = 1) or imaginary (if ip = 2) part of array boxr
! on node 0 (if fparafft = .true.) or array boxr_d distributed over
! nodes (if fparafft = .false.). Gaussian Cube file is tested using
! lattice vectors a and lattice constant al (in Bohr), FFTgrid is FFT
! grid size, Nplane is number of FFT xy-planes per node, ierr is
! return error code (0 means success).
!
!===============================================================================
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
module read_cube_m
  use global_m
  implicit none
  private
  public :: read_cube
  contains
subroutine read_cube(fparafft,unum,fnam,a,al,ip,Nplane,FFTgrid,boxr,boxr_d,ierr)
  logical, intent(in) :: fparafft
  integer, intent(in) :: unum
  character(len=256), intent(in) :: fnam
  real(DP), intent(in) :: al
  real(DP), intent(in) :: a(3,3)
  integer, intent(in) :: ip
  integer, intent(in) :: Nplane
  integer, intent(in) :: FFTgrid(3)
  complex(DPC), pointer, intent(inout) :: boxr(:,:,:) !< (FFTgrid(1),FFTgrid(2),FFTgrid(3))
  complex(DPC), pointer, intent(inout) :: boxr_d(:,:,:) !< (FFTgrid(1),FFTgrid(2),Nplane)
  integer, intent(out) :: ierr
  real(DP), parameter :: eps4 = 1.0d-4
  integer :: i,j,k,k0,k1,k2,jerr,knum,na,ngrid(3)
  real(DP) :: dr,origin(3),step(3,3)
  real(DP), allocatable :: buffer(:)
  character(len=256) :: tmpstr
 
  if (peinf%inode.eq.0) then
    na=0
    ngrid(:)=0
    origin(:)=0.0d0
    step(:,:)=0.0d0
    call open_file(unit=unum,file=fnam,status='old',form='formatted')
    read(unum,*,iostat=jerr)
    if (jerr.eq.0) read(unum,*,iostat=jerr)
    if (jerr.eq.0) read(unum,*,iostat=jerr)na,(origin(j),j=1,3)
    do i=1,3
      if (jerr.eq.0) read(unum,*,iostat=jerr)ngrid(i),(step(j,i),j=1,3)
    enddo
    if (jerr.eq.0) call close_file(unit=unum)
  endif
  if (ngrid(1).ne.FFTgrid(1).or.ngrid(2).ne.FFTgrid(2).or. &
    ngrid(3).ne.FFTgrid(3)) jerr=1
  dr=0.0d0
  do j=1,3
    dr=dr+abs(origin(j))
  enddo
  dr=dr/dble(3)
  if (dr.gt.eps4) jerr=1
  dr=0.0d0
  do i=1,3
    do j=1,3
      dr=dr+abs(dble(ngrid(i))*step(j,i)-al*a(j,i))
    enddo
  enddo
  dr=dr/dble(9)
  if (dr.gt.eps4) jerr=1
  if (jerr .ne. 0) write(0,*) 'WARNING: Inconsistency in your calculation and .cube file.'
  if (jerr .ne. 0) write(0,*) 'al',al
  if (jerr .ne. 0) write(0,*) 'a',a
  if (jerr .ne. 0) write(0,*) 'ngrid',ngrid
  if (jerr .ne. 0) write(0,*) 'FFTgrid',FFTgrid
  if (jerr .ne. 0) write(0,*) 'step',step
  if (jerr.eq.0) then
    allocate(buffer (FFTgrid(3)))
    if (peinf%inode.eq.0) then
      if (mod(FFTgrid(3),6).eq.0) then
        knum=FFTgrid(3)/6
      else
        knum=FFTgrid(3)/6+1
      endif
      call open_file(unit=unum,file=fnam,status='old',form='formatted')
      do i=1,6+na
        read(unum,*)
      enddo
    endif
    do i=1,FFTgrid(1)
      do j=1,FFTgrid(2)
        if (peinf%inode.eq.0) then
          do k=1,knum
            read(unum,103)tmpstr
            k1=6*(k-1)+1
            k2=6*(k-1)+6
            if (k2.gt.FFTgrid(3)) k2=FFTgrid(3)
            read(tmpstr,*)(buffer(k0),k0=k1,k2)
          enddo
        endif
        if (fparafft) then
          if (ip.eq.1) then
            do k=1,FFTgrid(3)
              if (k.ge.Nplane*peinf%inode+1.and.k.le.Nplane*(peinf%inode+1)) &
                boxr_d(i,j,k-Nplane*peinf%inode)=cmplx(buffer(k),aimag(boxr_d(i,j,k-Nplane*peinf%inode)),kind=DPC)
            enddo
          elseif (ip.eq.2) then
            do k=1,FFTgrid(3)
              if (k.ge.Nplane*peinf%inode+1.and.k.le.Nplane*(peinf%inode+1)) &
                boxr_d(i,j,k-Nplane*peinf%inode)=cmplx(dble(boxr_d(i,j,k-Nplane*peinf%inode)),buffer(k),kind=DPC)
            enddo
          endif
        else
          if (ip.eq.1) then
            do k=1,FFTgrid(3)
              boxr(i,j,k)=cmplx(buffer(k),aimag(boxr(i,j,k)),kind=DPC)
            enddo
          elseif (ip.eq.2) then
            do k=1,FFTgrid(3)
              boxr(i,j,k)=cmplx(dble(boxr(i,j,k)),buffer(k),kind=DPC)
            enddo
          endif
        endif
      enddo
    enddo
    if (peinf%inode.eq.0) then
      call close_file(unit=unum)
    endif
    if(allocated(buffer))then;deallocate(buffer);endif
  endif
  ierr=jerr
 
  return
103 format(a)
end subroutine read_cube
end module read_cube_m
