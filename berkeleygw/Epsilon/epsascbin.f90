!=========================================================================
!
! Utilities:
!
! (1) epsascbin() Originally by JLL Last Modified 5/5/2008 (JRD)
!
! This utility converts a ascii epsmat file to an binary epsmat file.
! It uses the input file epsconv.inp
!
!=========================================================================

program epsascbin
  use global_m
  implicit none
! JRD adate should be 11 characters in v. 2
  character :: ajname*6,adate*11,outfile*20
  real(DP) :: ecuts1,ecuts2
  real(DP), allocatable :: dFreqGrid(:)
  Complex(DPC), allocatable :: dFreqBrd(:)
  integer :: iunit,ig,im,iq,istart,j,jm, &
    nfiles,ng1,ng2,ngq1,ngq2,nmtx1,nmtx2,nq,nqtot, &
    ii,qgrid(3),freq_dep,nFreq,ijk
  character*120, allocatable :: filename(:)
  real(DP), allocatable :: ekin(:),q2(:,:)
  real(DP), allocatable :: eps(:)
  real(DP) :: tmpc
  integer, allocatable :: isort1(:),isort2(:),kx(:),ky(:),kz(:)
  real(DP) :: tmp
  integer :: itmp
  logical :: file_exists
  write(6,*) 'This routine should only be used on non-HDF5-based epsmat files.'
  write(6,*) 'If you built with HDF5 support (which is ideal), or otherwise'
  write(6,*) 'have HDF5-based epsmat files, you can use the h5dump command'
  write(6,*) 'to see your data in ascii format.'
  call open_file(55,file='epsconv.inp',form='formatted',status='old')
  read(55,*) nqtot
  allocate(q2 (3,nqtot))
  read(55,'(a20)') outfile
  write(6,*) 'Output -> ',outfile
  write(6,*)
!-------------------------
! Find maximal values, and check consistency between
! the input file and the epsmat files...
  istart=1
  ng1=0
  ngq1=0
  nmtx1=0
  read(55,*) nfiles
  allocate(filename (nfiles))
  do iunit=1,nfiles
    read(55,'(a20)') filename(iunit)
    write(6,*) 'Checking file ',TRUNC(filename(iunit))
    call open_file(unit=11,file=filename(iunit),form='formatted',status='old')
    read(11,'(1x,a6,1x,a11)') ajname,adate
    if(ajname /= 'chiGG0') then
      call die("Incorrect header '" // ajname // "' (must be 'chiGG0') in file '" // TRUNC(filename(iunit)) // "'")
    endif
    read(11,*) freq_dep,nFreq
    if (freq_dep.ne.0) then
      call die('epsascbin: freq_dep')
    endif
    read(11,*) (qgrid(ii),ii=1,3)
    if (freq_dep.eq.2) then
      allocate(dFreqGrid (nFreq))
      allocate(dFreqBrd (nFreq))
      read(11,*) (dFreqGrid(ijk),ijk=1,nFreq),(dFreqBrd(ijk),ijk=1,nFreq)
    else
      read(11,*)
    endif
    read(11,*)
    read(11,*)
    read(11,*) ecuts2
    if(iunit == 1) then
      ecuts1 = ecuts2
    else
      if(ecuts2.ne.ecuts1) then
        write(0,*) 'The cut-off in previous file (',ecuts1,') does not match ', &
          'the one in file ',TRUNC(filename(iunit)),' (',ecuts2,').'
        call die('epsascbin cutoff mismatch')
      endif
    endif
    read(11,*) nq,((q2(j,iq),j=1,3),iq=istart,istart+nq-1)
    read(11,*) ng2, (tmp,ig=1,3*ng2)
    if(ng1.eq.0) ng1=ng2
    if(ng2.ne.ng1) then
      call die('The number of G-vectors differs in epsmat files')
    endif
    do iq=istart,istart+nq-1
      read(11,*) ngq2,nmtx2,(itmp, ig=1,2*ngq2)
      if(ngq1.lt.ngq2) ngq1=ngq2
      if(nmtx1.lt.nmtx2) nmtx1=nmtx2
      read(11,*) (tmp,ig=1,ngq2)
      read(11,*) (q2(j,iq),j=1,3)
      do jm = 1, nmtx2
        read(11,*) (tmpc,im=1,nmtx2)
      enddo
    enddo
    istart=istart+nq
    call close_file(11)
  enddo
  if(istart-1.ne.nqtot) then
    write(0,*) 'found = ', istart - 1, ' expected ', nqtot
    call die('Number of q-vectors found differs from number in input file.')
  endif
  allocate(kx (ng1))
  allocate(ky (ng1))
  allocate(kz (ng1))
  call open_file(unit=11,file=filename(1),form='formatted',status='old')
  read(11,*)
  read(11,*) freq_dep,nFreq
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*) itmp,(tmp,iq=1,3*itmp)
  read(11,*) ng1,(kx(ig),ky(ig),kz(ig),ig=1,ng1)
  call close_file(11)
  allocate(isort1 (ng1))
  allocate(isort2 (ng1))
  allocate(ekin (ngq1))
  allocate(eps (nmtx1))
  call open_file(unit=12,file=outfile,form='unformatted',status='replace')
  write(12) ajname,adate
  write(12) freq_dep,nFreq
  write(12) (qgrid(ii),ii=1,3)
  if (freq_dep .eq. 2) then
    write(12) (dFreqGrid(ijk),ijk=1,nFreq),(dFreqBrd(ijk),ijk=1,nFreq)
  else
    write(12)
  endif
  write(12)
  write(12)
  write(12) ecuts1
  write(12) nqtot,((q2(j,iq),j=1,3),iq=1,nqtot)
  write(12) ng1,(kx(ig),ky(ig),kz(ig),ig=1,ng1)
  write(6,*)
  istart=1
  do iunit=1,nfiles
    write(6,*) 'Dealing with file ',TRUNC(filename(iunit))
    call open_file(unit=11,file=filename(iunit),form='formatted',status='old')
    read(11,*)
    read(11,*) freq_dep,nFreq
    read(11,*)
    read(11,*)
    read(11,*)
    read(11,*)
    read(11,*)
    read(11,*) nq,(tmp,iq=1,3*nq)
! write(6,*) 'Number of qs', nq
    read(11,*) itmp,(tmp,iq=1,3*itmp)
    do iq=istart,istart+nq-1
      write(6,'(a,f9.6,3x,f9.6,3x,f9.6)') ' -> q=',(q2(j,iq),j=1,3)
      read(11,*) ngq2,nmtx2,(isort1(ig),isort2(ig),ig=1,ngq2)
      write(12) ngq2,nmtx2,(isort1(ig),isort2(ig),ig=1,ngq2)
      read(11,*) (ekin(ig),ig=1,ngq2)
      write(12) (ekin(ig),ig=1,ngq2)
      read(11,*)
      write(12) (q2(j,iq),j=1,3)
      do jm =1, nmtx2
        read(11,*) (eps(im),im=1,nmtx2)
        write(12) (eps(im),im=1,nmtx2)
      enddo
    enddo
    istart=istart+nq
    call close_file(11)
  enddo
  call close_file(12)
end program epsascbin
