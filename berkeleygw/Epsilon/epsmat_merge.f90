!=========================================================================
!
! Utilities:
!
! (1) epsmat_merge() Originally by ? Last Modified 5/5/2008 (JRD)
!
! This utility merges a list of epsmat files into one file.
! It uses the input file epsmat_merge.inp. See example in current directory.
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
program epsmat_merge
  use global_m
  implicit none
  character :: ajname*6,adate*11,atime*14
  real(DP) :: dtol,div,ecuts1,ecuts2
  real(DP), allocatable :: dFreqGrid(:)
  complex(DPC), allocatable :: dFreqBrd(:)
  integer i,ig,im,iq,istart,j,k,jm,jj, &
    nfiles,ng1,ng2,ngq1,ngq2,nmtx1,nmtx2,nq,nqtot, &
    freq_dep,nFreq,qgrid(3),ijk
  character*20, allocatable :: filename(:)
  real(DP), allocatable :: ekin(:),q1(:,:),q2(:,:)
  complex(DPC), allocatable :: epsDyn(:)
  real(DP), allocatable :: eps(:)
  integer, allocatable :: isort1(:),isort2(:),kx(:),ky(:),kz(:)
  logical :: subspace, matrix_in_subspace_basis, keep_full_eps_static
  logical :: subspace_new, matrix_in_subspace_basis_new, keep_full_eps_static_new
  integer :: neig_sub, ierr
  real(DP), allocatable :: vcoul(:)
  complex(DPC), allocatable :: epsDyn_aux(:)
  dtol = 1.0d-5
  write(6,*) 'This routine should only be used on non-HDF5-based epsmat files.'
  write(6,*) 'If you built with HDF5 support (which is ideal), or otherwise'
  write(6,*) 'have HDF5-based epsmat files, you should use the routine'
  write(6,*) 'epsmat_hdf5_merge.py instead.'
  call open_file(55,file='epsmat_merge.inp',form='formatted',status='old')
  read(55,*) ecuts1,nqtot
  allocate(q1 (3,nqtot))
  allocate(q2 (3,nqtot))
  do iq=1,nqtot
    read(55,*) q1(1,iq),q1(2,iq),q1(3,iq),div
    do j=1,3
      q1(j,iq)=q1(j,iq)/div
    enddo
  enddo
!---------------------------
! Find maximal values, and check consistency between
! the input file and the epsmat files...
  subspace = .FALSE.
  matrix_in_subspace_basis = .FALSE.
  keep_full_eps_static = .FALSE.
  istart=1
  ng1=0
  ngq1=0
  nmtx1=0
  read(55,*) nfiles
  allocate(filename (nfiles))
  do i=1,nfiles
    read(55,'(a20)') filename(i)
    write(6,*) 'Checking file ',filename(i)
    call open_file(unit=11,file=filename(i),form='unformatted',status='old')
    subspace_new = .FALSE.
    matrix_in_subspace_basis_new = .FALSE.
    keep_full_eps_static_new = .FALSE.
    read(11) ajname,adate
    ierr = 0
    read(11,IOSTAT=ierr) freq_dep,nFreq, subspace_new, matrix_in_subspace_basis_new, keep_full_eps_static_new
    IF(ierr .NE. 0) THEN
      subspace_new = .FALSE.
      matrix_in_subspace_basis_new = .FALSE.
      keep_full_eps_static_new = .FALSE.
    END IF
    !XXXX
    IF(i > 1) THEN
      IF((subspace .NEQV. subspace_new) .OR. &
         (matrix_in_subspace_basis .NEQV. matrix_in_subspace_basis_new) .OR. &
         (keep_full_eps_static .NEQV. keep_full_eps_static_new)) THEN
         write(0,*) 'Inconsistent definition of subspace parameters.'
         write(0,*) 'For file', filename(i)
         write(0,*) 'subspace', subspace_new
         write(0,*) 'matrix_in_subspace_basis', matrix_in_subspace_basis_new
         write(0,*) 'keep_full_eps_static', keep_full_eps_static_new
         write(0,*) 'For file', filename(i-1)
         write(0,*) 'subspace', subspace
         write(0,*) 'matrix_in_subspace_basis', matrix_in_subspace_basis
         write(0,*) 'keep_full_eps_static', keep_full_eps_static
         call die('epsmat_merge subspace parameters mismatch')
      END IF
    END IF
    subspace = subspace_new
    matrix_in_subspace_basis = matrix_in_subspace_basis_new
    keep_full_eps_static = keep_full_eps_static_new
    !XXXX
    read(11) (qgrid(j),j=1,3)
    if (freq_dep .ne. 0 .and. i .eq. 1) then
      allocate(dFreqGrid (nFreq))
      allocate(dFreqBrd (nFreq))
      read(11) (dFreqGrid(ijk),ijk=1,nFreq),(dFreqBrd(ijk),ijk=1,nFreq)
    else
      read(11)
    endif
    read(11)
    read(11)
    read(11) ecuts2
    if(ecuts2.ne.ecuts1) then
      write(0,*) 'The cut-off in input file (',ecuts1,') does not match the one in file ',filename(i),' (',ecuts2,').'
      call die('epsmat_merge cutoff mismatch')
    endif
    read(11) nq,((q2(j,iq),j=1,3),iq=istart,istart+nq-1)
    do iq=istart,istart+nq-1
!$$ tolerance
      if((abs(q2(1,iq)-q1(1,iq)).gt.dtol).or. &
        (abs(q2(2,iq)-q1(2,iq)).gt.dtol).or. &
        (abs(q2(3,iq)-q1(3,iq)).gt.dtol)) then
        write(0,*) 'The q-vector ',iq,' in input file (', (q1(j,iq),j=1,3), &
          ') does not match does the one in file ', filename(i),' (',(q2(j,iq),j=1,3),').'
        call die('epsmat_merge q-vector mismatch')
      endif
    enddo
    read(11) ng2
    if(ng1.eq.0) ng1=ng2
    if(ng2.ne.ng1) then
      call die('The number of G-vectors differs in epsmat files')
    endif
    do iq=istart,istart+nq-1
      read(11) ngq2,nmtx2
      if(ngq1.lt.ngq2) ngq1=ngq2
      if(nmtx1.lt.nmtx2) nmtx1=nmtx2
      read(11)
      read(11) (q2(j,iq),j=1,3)
      IF(subspace .AND. matrix_in_subspace_basis) THEN
        read(11) ! vcoul
        IF(keep_full_eps_static) THEN
          do j=1, nmtx2
            read(11) ! For epsRDyn only
          end do
        END IF
        ! eigenvectors
        read(11) neig_sub
        do j=1, nmtx2
          read(11) ! For epsRDyn only
        end do
        ! subspace matrices
        do j=1, neig_sub
          do k=1, neig_sub
            read(11)
          end do
        end do
      ELSE ! subspace
        if(freq_dep.eq.0) then
          do j = 1, nmtx2
            read(11)
          enddo
        else
          do j = 1, nmtx2
            do k = 1, nmtx2
              read(11)
            enddo
          enddo
        endif
      END IF ! subspace
!$$ tolerance
      if((abs(q2(1,iq)-q1(1,iq)).gt.dtol).or. &
        (abs(q2(2,iq)-q1(2,iq)).gt.dtol).or. &
        (abs(q2(3,iq)-q1(3,iq)).gt.dtol)) then
        write(0,*) 'The q-vector ',iq,' in input file (', (q1(j,iq),j=1,3), &
          ') does not match does the one in file ', filename(i),' (',(q2(j,iq),j=1,3),').'
        call die('epsmat_merge q-vector mismatch')
      endif
    enddo
    istart=istart+nq
    call close_file(11)
  enddo
  if(istart-1.ne.nqtot) then
    call die('Could not find all q-vectors in epsmat files')
  endif
  allocate(kx (ng1))
  allocate(ky (ng1))
  allocate(kz (ng1))
  call open_file(unit=11,file=filename(1),form='unformatted',status='old')
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11)
  read(11) ng1,(kx(i),ky(i),kz(i),i=1,ng1)
  call close_file(11)
  allocate(isort1 (ng1))
  allocate(isort2 (ng1))
  allocate(ekin (ngq1))
  if(freq_dep.eq.0) then
    allocate(eps (nmtx1))
  else
    allocate(epsDyn (nFreq))
    !XXXXX
    IF(subspace .AND. matrix_in_subspace_basis) THEN
      allocate(vcoul (nmtx1))
      allocate(epsDyn_aux (nmtx1))
    END IF
    !XXXXX
  endif
  ajname='chiGG0'
  call date_time(adate,atime)
  call open_file(unit=12,file='epsmat',form='unformatted',status='replace')
  write(12) ajname,adate
  write(12) freq_dep,nFreq, subspace, matrix_in_subspace_basis, keep_full_eps_static
  write(12) (qgrid(i),i=1,3)
  if (freq_dep .eq. 2) then
    write(12) (dFreqGrid(i),i=1,nFreq),(dFreqBrd(i),i=1,nFreq)
  else
    write(12)
  endif
  write(12)
  write(12)
  write(12) ecuts1
  write(12) nqtot,((q1(j,iq),j=1,3),iq=1,nqtot)
  write(12) ng1,(kx(ig),ky(ig),kz(ig),ig=1,ng1)
  write(6,*)
  istart=1
  do i=1,nfiles
    write(6,*) 'Dealing with file ',filename(i)
    call open_file(unit=11,file=filename(i),form='unformatted',status='old')
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11) nq
    read(11)
    do iq=istart,istart+nq-1
      write(6,'(a,f9.6,3x,f9.6,3x,f9.6)') ' -> q=',(q1(j,iq),j=1,3)
      read(11) ngq2,nmtx2,(isort1(ig),isort2(ig),ig=1,ngq2)
      write(12) ngq2,nmtx2,(isort1(ig),isort2(ig),ig=1,ngq2)
      read(11) (ekin(ig),ig=1,ngq2)
      write(12) (ekin(ig),ig=1,ngq2)
      read(11)
      write(12) (q1(j,iq),j=1,3)
      IF(subspace .AND. matrix_in_subspace_basis) THEN
        ! vcoul
        read(11) (vcoul(im),im=1,nmtx2)
        write(12) (vcoul(im),im=1,nmtx2)
        IF(keep_full_eps_static) THEN
          do jm = 1, nmtx2
            read(11) (epsDyn_aux(im),im=1,nmtx2) ! For epsRDyn only
            write(12) (epsDyn_aux(im),im=1,nmtx2)
          enddo
        END IF
        read(11) neig_sub
        write(12) neig_sub
        ! eigenvectors
        do jm = 1, neig_sub
          read(11) (epsDyn_aux(im),im=1,nmtx2)
          write(12) (epsDyn_aux(im),im=1,nmtx2)
        end do
        do jm = neig_sub + 1, nmtx2
          read(11)
          write(12)
        end do
        ! subspace matrices
        do jm = 1, neig_sub
          do im = 1, neig_sub
            read(11) (epsDyn(jj),jj=1,nFreq)
            write(12) (epsDyn(jj),jj=1,nFreq)
          enddo
        enddo
      ELSE ! subspace
        if(freq_dep.eq.0) then
          do jm = 1, nmtx2
            read(11) (eps(im),im=1,nmtx2)
            write(12) (eps(im),im=1,nmtx2)
          enddo
        else
          do jm = 1, nmtx2
            do im = 1, nmtx2
              read(11) (epsDyn(jj),jj=1,nFreq)
              write(12) (epsDyn(jj),jj=1,nFreq)
            enddo
          enddo
        endif
      END IF ! subspace
    enddo
    istart=istart+nq
    call close_file(11)
  enddo
end program epsmat_merge
