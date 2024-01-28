!==============================================================================
!
! Utilities:
!
! (1) summarize_eigenvectors() Originally By JRD,CHP Last Modified: 3/16/2018 (GKA)
!
! This program prints some useful information about the exciton wavefunction.
!
!==============================================================================
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
program summarize_eigenvectors
  use global_m
  use evecs_m
  implicit none
  type(evecs_t) :: evecs
  integer :: error
  integer :: ns, nc, nv, nk, ii, ic, iv, is, ik, ikmax
  integer :: ijk, nexc, ninfile
  integer :: up=10
  logical :: tda
  real(DP) :: energy, weight, this_weight, wmax, emax, emin
  real(DP), allocatable :: energies(:)
!
  character*20, allocatable :: filename(:)
  call open_file(unit=up,file='summarize_eigenvectors.inp',status='old')
  read(up,*) tda
  read(up,*) ninfile
  read(up,*) emin, emax
  read(up,*) nexc
  if (nexc .gt. 0) then
    allocate(energies (nexc))
    allocate(filename (nexc))
    do ijk =1, nexc
      read(up,*) energies(ijk)
    enddo
    do ijk =1, nexc
      write(filename(ijk),'(a8,i2.2)') 'exciton_',ijk
    enddo
  endif
  close (up)
  ! Set up evecs (without calling init)
  evecs%tda = tda
  if (evecs%use_hdf5) then
    INQUIRE(file=evecs%default_fname_hdf5, exist=evecs%use_hdf5)
    if (.not. evecs%use_hdf5) then
      write(6,'(a,a)') 'HDF5 file not found: ', trim(evecs%default_fname_hdf5)
      write(6,'(a,a)') 'Falling back on binary file: ', trim(evecs%default_fname_bin)
    end if
  end if
  ! Open eigenvector file and read header
  write(6,'(a)')
  write(6,'(a)') 'Reading eigenvectors'
  call evecs%open_read()
  call evecs%read_header(broadcast=.false.)
  ! Print out information found in file
  call evecs%print_out_header_info(6)
  ns=evecs%ns; nk=evecs%nk; nv=evecs%nv; nc=evecs%nc
  if (ninfile .eq. 0) ninfile = evecs%neig
  ! We will be reading a single eigenvector at a time
  evecs%meig = 1
  ! Allocate memory
  call evecs%alloc(with_Avc=.true.)
  write(6,'(a)')
  write(6,'(a)') 'exciton energies follow (eV)'
  write(6,'(a)') 'wtot = sum_k |A_vck|^2. wmax = max_k |A_vck|^2. |A_vc (ikmax)|^2 = wmax.'
  do ii=1,ninfile
    call evecs%read_next_eigenvector()
    energy = evecs%evals(1)
    if ((energy > emax) .or. (energy < emin)) cycle
    call evecs%reshape_Avc()
    write(6,'(a)')
    write(6,'(a)')
    write(6,'(a,i5,f10.5)') ' Special analysis for state ',ii,energy
    write(6,'(2a5,3a10)') 'c','v','wtot','wmax','ikmax'
    do ic=1,nc
      do iv=1,nv
        weight = 0.0d0
        wmax = 0.0d0
        ikmax = -1
        do ik=1,nk
          do is=1,ns
            if ( tda ) then
              this_weight = abs(evecs%Avc(is,iv,ic,ik,1))**2
            else
              this_weight = (evecs%Avc_l(is,iv,ic,ik,1))*evecs%Avc(is,iv,ic,ik,1) + &
                            (evecs%Bvc_l(is,iv,ic,ik,1))*evecs%Bvc_r(is,iv,ic,ik,1)
            end if
            weight=weight+this_weight
            if (this_weight > wmax) then
              wmax = this_weight
              ikmax = ik
            endif
          enddo
        enddo
        write(6,'(2i5,2f10.5,i10)') ic, iv, weight, wmax, ikmax
      enddo
    enddo
    do ijk = 1, nexc
      if (abs(energy - energies(ijk)) .le. 1d-5) then
        write(6,'(a,i6,f12.6)') 'Calculating A(k) for :',ijk,energies(ijk)
        call open_file(unit=200+ijk,file=filename(ijk),status='replace')
        write(200+ijk, '(a,a9,2a10,a18)') '#', 'kx', 'ky', 'kz', 'sum |A(k)|^2'
        do ik=1,nk
          weight = 0.0d0
          do ic=1,nc
            do iv=1,nv
              do is=1,ns
                if ( tda ) then
                  this_weight = abs(evecs%Avc(is,iv,ic,ik,1))**2
                else
                  this_weight = (evecs%Avc_l(is,iv,ic,ik,1))*evecs%Avc(is,iv,ic,ik,1) + &
                                (evecs%Bvc_l(is,iv,ic,ik,1))*evecs%Bvc_r(is,iv,ic,ik,1)
                end if
                weight = weight + this_weight
              enddo
            enddo
          enddo
          write(200+ijk,'(3f10.5,f18.5)') evecs%kpts(1:3,ik),weight
        enddo
        call close_file(200+ijk)
      endif
    enddo
  enddo
  write(6,'(a)')
  call evecs%close_file()
  call evecs%free()
end program summarize_eigenvectors
