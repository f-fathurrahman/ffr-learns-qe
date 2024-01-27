!=========================================================================
!
! Utilities:
!
! (1) eps0sym Originally by MJ Last Modified 1/12/2010 (MJ)
!
! This utility symmetrizes eps0mat file.
!
!=========================================================================
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
program eps0sym
  use global_m
  use epsread_hdf5_m
  use epswrite_hdf5_m
  use write_matrix_m
  implicit none
  character :: ajname*6,adate*11,outfile*80,infile*80
  real(DP) :: ecuts
  real(DP), allocatable :: dFreqGrid(:), epsdiag(:,:,:)
  complex(DPC), allocatable :: dFreqBrd(:)
  integer :: i,ig,j, &
    ng,ngq,nmtx, &
    ii,qgrid(3),freq_dep,nFreq,nfreq_imag,&
    jj,kk,ll,nargs,imap,iout, &
    nq, gx, gy, gz, gmax, &
    ig0, igx, igy, igz, jgx, jgy, jgz, jout
  real(DP), allocatable :: ekin(:)
  real(DP) :: q(3,1), qk(3,1), errorMax
  real(DP), allocatable :: eps(:,:), tempeps(:,:)
  real(DP) :: errorMaxTemp
  integer, allocatable :: isort(:),isorti(:),kx(:),ky(:),kz(:)
  integer, allocatable :: minusgidx(:), map(:), old(:,:)
  integer :: nmtx0_of_q(1), isize, error
  logical :: use_hdf5
  write(6,*) "Real" // ' version is used to symmetrize file'
!------------------
! Get file names from command-line arguments
  nargs = command_argument_count()
  if (nargs .ne. 2) then
    call die('Usage: eps0sym eps0mat_in eps0mat_out')
  endif
  call get_command_argument(1,infile)
  call get_command_argument(2,outfile)
  use_hdf5 = .false.
    call open_file(unit=11,file=TRUNC(infile),form='unformatted',status='old')
!--------------------------
! Read in the eps0mat file
!
    read(11) ajname,adate
    read(11) freq_dep,nFreq
    if (freq_dep.ne.0) then
      call die('eps0sym: Full frequency not supported')
    endif
    read(11) (qgrid(ii),ii=1,3)
    if (freq_dep .eq. 2) then
      allocate(dFreqGrid (nFreq))
      allocate(dFreqBrd (nFreq))
      read(11) (dFreqGrid(i),i=1,nFreq),(dFreqBrd(i),i=1,nFreq)
    else
      read(11)
    endif
    read(11)
    read(11)
    read(11) ecuts
    read(11) nq,(q(j,1),j=1,3)
    !if (nq .ne. 1) then
    ! call die('This only works for the q->0 point.')
    !endif
    read(11) ng
    read(11) ngq,nmtx
    call close_file(11)
  allocate(kx (ng))
  allocate(ky (ng))
  allocate(kz (ng))
  allocate(isort (ngq))
  allocate(isorti (ngq))
  allocate(ekin (ngq))
  allocate(eps (nmtx,nmtx))
    call open_file(unit=11,file=TRUNC(infile),form='unformatted',status='old')
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11) ng,(kx(i),ky(i),kz(i),i=1,ng)
    read(11) ngq,nmtx,(isort(ig),isorti(ig),ig=1,ngq)
    read(11) (ekin(ig),ig=1,ngq)
    read(11) (qk(j,1),j=1,3)
    do jj = 1, nmtx
      read(11) (eps(ii,jj),ii=1,nmtx)
    enddo
    call close_file(11)
  ! Since we want the q=0 dielectric function but we have the
  ! q0<>0 but small dielectric, we can average over q0 and -q0
  ! to get a better dielectric (linear terms in q0 will be canceled)
! Calculate the maximum
  gmax = 0
  do ii = 1,ng
    if (abs(kx(ii)) > gmax) gmax = abs(kx(ii))
    if (abs(ky(ii)) > gmax) gmax = abs(ky(ii))
    if (abs(kz(ii)) > gmax) gmax = abs(kz(ii))
  enddo
! Create a map
! write(6,*) gmax
  allocate(map ((2*gmax+1)*(2*gmax+1)*(2*gmax+1)))
  map = 0
  do ii = 1, ngq
    iout = isort(ii)
    if (iout.eq.0) cycle
    gx = kx(iout)
    gy = ky(iout)
    gz = kz(iout)
    imap = ((gx+gmax)*(2*gmax+1)+gy+gmax)*(2*gmax+1)+gz+gmax+1
    map(imap) = ii
  enddo
! For each g, find -g in the gvector list
! and also mark which ii corresponds to G=0
  allocate(minusgidx (nmtx))
  minusgidx = 0
  do ii=1,nmtx
    iout = isort(ii)
    gx = kx(iout)
    gy = ky(iout)
    gz = kz(iout)
    imap = ((-gx+gmax)*(2*gmax+1)-gy+gmax)*(2*gmax+1)-gz+gmax+1
    minusgidx(ii) = map(imap)
    if (gx .eq. 0 .and. gy .eq. 0 .and. gz .eq. 0) ig0 = ii
  enddo
! do ii = 1, min(100,nmtx)
! iout = isort(ii)
! gx = kx(iout)
! gy = ky(iout)
! gz = kz(iout)
! mgx = kx(isort(minusgidx(ii)))
! mgy = ky(isort(minusgidx(ii)))
! mgz = kz(isort(minusgidx(ii)))
! write(6,'(7i4)') ii, gx, gy , gz, mgx, mgy, mgz
! enddo
! Set the wings to zero
! This is as per Baldereschi and Tosatti, PRB 17, 4710 (1978)
! This is perhaps not the correct thing to do. One should
! still symmetrize epsilon - but the wings should come out
! whatever they need to be automatically. What is given in that
! paper by Baldereschi and Tosatti is at q=0 and the dielectric
! function is weird at that point. What is needed in the GW
! code is the average in the minibz...
  !do ii=1,nmtx
  ! do jj=1,nmtx
  ! if (ii .eq. ig0 .and. jj .eq. ig0) cycle
  ! if (ii .eq. ig0) eps(ii,jj) = 0.0d0
  ! if (jj .eq. ig0) eps(ii,jj) = 0.0d0
  ! enddo
  !enddo
! Copy eps(q0) into a temporary
  allocate(tempeps (nmtx,nmtx))
  tempeps = eps
  errorMax=0D0
  errorMaxTemp=0D0
  write(6,'(5a4,3x,3a4,1x)',advance='no') "ig", "ig'", "Gx", "Gy", "Gz", "G'x", "G'y", "G'z"
  ! handle the fact that the two components of the complex numbers will be written in the complex case
  write(6,'(3a16)') "eps(G,G')", "eps(-G,-G')*", "difference"
! do ii = 1, min(100,nmtx)
! do jj = 1, min(100,nmtx)
  do ii = 1, nmtx
    do jj = 1, nmtx
      iout = isort(ii)
      igx = kx(iout)
      igy = ky(iout)
      igz = kz(iout)
      jout = isort(jj)
      jgx = kx(jout)
      jgy = ky(jout)
      jgz = kz(jout)
      kk = minusgidx(ii)
      ll = minusgidx(jj)
      if ((kk .le. nmtx) .and. (ll .le. nmtx)) then
        errorMaxTemp=eps(ii,jj)-(eps(kk,ll))
        if (abs(errorMaxTemp) .gt. errorMax) errorMax = abs(errorMaxTemp)
        if (ii .lt. 20 .and. jj .lt. 20) then
          write(6,'(5i4,3x,3i4,1x,6E16.5)') ii,jj,igx,igy,igz,jgx,jgy,jgz,eps(ii,jj),(eps(kk,ll)),errorMaxTemp
        endif
      endif
    enddo
  enddo
  write(6,*) "Error = eps(G,G')-eps*(-G,-G')"
  write(6,'("The max error in your matrix is",E12.5)') errorMax
  write(6,*) "Symmetrizing the matrix"
! Now add in contribution from -q0 which means conjg(eps(-g,-gp))
  do ii=1,nmtx
    do jj=1,nmtx
      kk = minusgidx(ii)
      ll = minusgidx(jj)
      if ((kk .le. nmtx) .and. (ll .le. nmtx)) then
        tempeps(ii,jj) = tempeps(ii,jj) + (eps(kk,ll))
      endif
    enddo
  enddo
! Average over q0 and -q0 and put back into eps
  eps = 0.5d0*tempeps
  if(allocated(tempeps))then;deallocate(tempeps);endif
  if(allocated(minusgidx))then;deallocate(minusgidx);endif
    call open_file(unit=20,file=TRUNC(outfile),form='unformatted',status='replace')
    ajname='chiGG0'
    write(20) ajname,adate
    write(20) freq_dep,nFreq
    write(20) (qgrid(ii),ii=1,3)
    if (freq_dep .eq. 2) then
      write(20) (dFreqGrid(i),i=1,nFreq),(dFreqBrd(i),i=1,nFreq)
    else
      write(20)
    endif
    write(20)
    write(20)
    write(20) ecuts
    write(20) nq,(q(j,1),j=1,3)
    write(20) ng,(kx(ig),ky(ig),kz(ig),ig=1,ng)
    write(20) ngq,nmtx,(isort(ig),isorti(ig),ig=1,ngq)
    write(20) (ekin(ig),ig=1,ngq)
    write(20) (qk(j,1),j=1,3)
    do jj = 1, nmtx
      write(20) (eps(ii,jj),ii=1,nmtx)
    enddo
    call close_file(20)
end program eps0sym
