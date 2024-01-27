!==============================================================================================
!
! Utilities:
!
! (1) ploteps Originally By JRD Last Modified 6/30/2008 (JRD)
!
!==============================================================================================
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
program ploteps
  use global_m
  implicit none
  integer, allocatable :: isrtq(:),oldx(:),oldy(:),oldz(:)
  integer, allocatable :: isrtold(:)
  real(DP) :: qvec(3),qk(3),gmax_in,garbaged
  integer :: ii, jj, nold, nge, nmtx, nrq0
  real(DP), allocatable :: eps(:,:)
  character :: ajname*6, adate*11
  write(*,*) 'Welcome to Epsilon Plotting'
  call open_file(unit=10,file='eps0mat',form='unformatted',status='old')
! call open_file(unit=11,file='epsmat',form='unformatted',status='old')
  call open_file(unit=42,file='eps_surface.out',status='replace')
  call open_file(unit=43,file='eps_diagonal.out',status='replace')
!---------------------------
! Now read eps0mat
  read(10)
  read(10) ii
  if (ii.ne.0) then
    call die('Full frequency dependence not supported')
  endif
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10) nold
  read(10) nge
  rewind(10)
! deallocate(oldx)
  allocate(oldx (nold))
! deallocate(oldy)
  allocate(oldy (nold))
! deallocate(oldz)
  allocate(oldz (nold))
! deallocate(ekold)
! allocate(ekold(nge))
! deallocate(isrtold)
  allocate(isrtold (nge))
! deallocate(isrtq)
! allocate(isrtq(nge))
! deallocate(isrtq)
  allocate(isrtq (nge))
  read(10) ajname,adate
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10) gmax_in
  read(10) nrq0,(qvec(ii),ii=1,3)
  read(10) nold,(oldx(ii),oldy(ii),oldz(ii),ii=1,nold)
! JRD: This should be changed later
  if(nrq0.gt.1) then
    call die("There is more than one q-point in epsmat", only_root_writes = .true.)
  endif
! Read q->0 dielectric matrix
  read(10) nge,nmtx,(isrtold(ii),jj,ii=1,nge)
  write(*,*) 'nge, nmtx =', nge,nmtx
  allocate(eps (nmtx,nmtx))
  write(*,*) ' '
  write(*,*) 'Allocated eps!', nmtx
  write(*,*) ' '
  read(10) (garbaged,ii=1,nge)
  read(10) (qk(ii),ii=1,3)
  do jj = 1, nmtx
    read(10) (eps(ii,jj),ii=1,nmtx)
  enddo
! read(10) (totalreal,ii=1,nge)
  write(*,*) 'read in eps'
  isrtq=0
  call close_file(10)
! call close_file(11)
  do ii = 1, nmtx
    do jj = 1, nmtx
      write(42,*) ii, jj, dble(eps(ii,jj))
      if (ii .eq. jj) then
        write(43,*) ii,oldx(isrtold(ii)),oldy(isrtold(ii)), &
          oldz(isrtold(ii)),dble(eps(ii,jj))
      endif
    end do
    write(42,*) ' '
  end do
  write(*,*) 'Finalizing'
end program ploteps
