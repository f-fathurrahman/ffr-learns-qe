!
! SIB: This routine looks the same as input().
! Except it reads from WFNq and
! writes to iunit_v='INT_VWFQ' and only valence bands.
! And k-point information is read into kpq.
!
! SUBROUTINE READS CRYSTAL DATA AND WAVEFUNCTIONS FROM TAPE26
! AND PARAMETERS FOR POLARIZABILITY CALCULATION FROM TAPE5
! TAPE10 (OUTPUT TAPE) IS INITIALIZED
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
module input_q_m
  use global_m
  use eqpcor_m
  use input_utils_m
  use misc_m
  use scissors_m
  use wfn_rho_vxc_io_m
  use io_utils_m
  use hdf5_io_m
  use wfn_io_hdf5_m
  implicit none
  private
  public :: input_q
contains
subroutine input_q(gvec,kpq,cwfn,vwfn,pol,intwfnvq)
  type (gspace), intent(in) :: gvec
  type (kpoints), intent(out) :: kpq
  type (conduction_wfns), intent(in) :: cwfn
  type (valence_wfns), intent(in) :: vwfn
  type (polarizability), intent(inout) :: pol
  type (int_wavefunction), intent(out) :: intwfnvq
  type (crystal) :: crys
  type (symmetry) :: syms
  real(DP) :: vcell
  integer :: dummygvec(1, 1)
  character :: fncor*32
  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvecq
 
  sheader = 'WFN'
  iflavor = 0
  if(pol%wfn_hdf5) then
  else
    if(peinf%inode == 0) then
      write(6,'(a)') ' Reading header of WFNq'
      call open_file(26,file='WFNq',form='unformatted',status='old')
    endif
    call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, syms, crys, &
      dont_warn_kgrid=pol%subsample, warn=.false.)
    call read_binary_gvectors(26, gvecq%ng, gvecq%ng, dummygvec, dont_read = .true.)
  endif
  call check_trunc_kpts(pol%icutv, kpq)
  call scissors_shift(kpq, pol%scis)
  call get_volume(vcell,crys%bdot)
  if (abs(crys%celvol-vcell).gt.TOL_Small) then
    call die('volume mismatch')
  endif
!-----------------------------------------------------------------
! If a quasi-particle correction file exists, read the corrected
! quasiparticle energies from file (in eV)
  if (pol%eqp_corrections) then
    fncor='eqp_q.dat'
    call eqpcor(fncor,peinf%inode,peinf%npes,kpq,1+vwfn%ncore_excl,vwfn%nband+vwfn%ncore_excl+pol%ncrit,0,0,kpq%el,kpq%el,&
                 kpq%el,1,0)
  endif
  call find_efermi(pol%rfermi, pol%efermi, pol%efermi_input, kpq, kpq%mnband, 1+vwfn%ncore_excl, &
    "shifted grid", should_search = .false., should_update = .false., write7 = .true.)
  if (peinf%inode==0) then
    if(any (kpq%ifmax(:,:) < vwfn%nband+vwfn%ncore_excl .or. kpq%ifmax(:,:) > vwfn%nband +vwfn%ncore_excl+ pol%ncrit)) then
      write(0,'(a,i6,a,i6,a)') 'epsilon.inp says there are ', vwfn%nband, ' fully occupied bands and ', &
        pol%ncrit, ' partially occupied.'
      write(0,'(a,2i6)') 'This is inconsistent with highest bands in WFNq file; min, max = ', minval(kpq%ifmax), maxval(kpq%ifmax)
      call die("band_occupation, number_partial_occup, and WFNq inconsistent.")
    endif
    if(maxval(kpq%ifmax) - minval(kpq%ifmax) > pol%ncrit) then
      write(0,'(a,i6,a)') 'epsilon.inp says there are ', pol%ncrit, ' partially occupied bands.'
      write(0,'(a,i6)') 'This is less than the number partially occupied in WFNq file: ', maxval(kpq%ifmax) - minval(kpq%ifmax)
      call die("number_partial_occup and WFNq inconsistent.")
    endif
  endif
  if (.not.pol%skip_chi) then
    if (pol%wfn_hdf5) then
    else
      call read_wavefunctions(kpq, gvec, pol, cwfn, vwfn, intwfnvq)
      if (peinf%inode.eq.0) then
        call close_file(26)
      endif
    endif
  endif
 
  return
end subroutine input_q
subroutine read_wavefunctions(kpq, gvec, pol, cwfn, vwfn, intwfnvq)
  type (kpoints), intent(in) :: kpq
  type (gspace), intent(in) :: gvec
  type (polarizability), intent(in) :: pol
  type (conduction_wfns), intent(in) :: cwfn
  type (valence_wfns), intent(in) :: vwfn
  type (int_wavefunction), intent(out) :: intwfnvq
  integer, allocatable :: isort(:)
  real(DP), allocatable :: zc(:,:)
  character :: filenamevq*20
  integer :: i,i2,j,k,ik,iiii,is
  integer :: iunit_v
  real(DP) :: qk(3)
  logical :: dont_read
  type(gspace) :: gvec_kpt
  type(progress_info) :: prog_info !< a user-friendly progress report
 
  allocate(intwfnvq%ng (kpq%nrk))
  allocate(intwfnvq%isort (kpq%ngkmax,kpq%nrk))
  allocate(intwfnvq%cg (kpq%ngkmax,kpq%nrk*peinf%nvownactual,kpq%nspin*kpq%nspinor))
  allocate(intwfnvq%qk (3,kpq%nrk))
  call progress_init(prog_info, 'reading wavefunctions (WFNq)', 'state', kpq%nrk*(vwfn%nband+pol%ncrit))
  do ik=1,kpq%nrk
    qk(1:3) = kpq%rk(1:3, ik)
    allocate(gvec_kpt%components (3, kpq%ngk(ik)))
    call read_binary_gvectors(26, kpq%ngk(ik), kpq%ngk(ik), gvec_kpt%components)
    allocate(isort (kpq%ngk(ik)))
    do i = 1, kpq%ngk(ik)
      call findvector(isort(i), gvec_kpt%components(:, i), gvec)
      if (isort(i) == 0) then
        if(peinf%inode == 0) write(0,*) 'ik = ', ik, 'ig = ', i, 'gvec = ', gvec_kpt%components(:, i)
        call die('input_q: could not find gvec')
      endif
    enddo
    if(associated(gvec_kpt%components))then;deallocate(gvec_kpt%components);nullify(gvec_kpt%components);endif
    intwfnvq%ng(ik)=kpq%ngk(ik)
    intwfnvq%isort(1:kpq%ngk(ik),ik)=isort(1:kpq%ngk(ik))
    intwfnvq%qk(:,ik)=qk(:)
!
! SIB: loop on max number of bands, and proc 0 reads the wave function
! from unit 26, checks normalization, and if the band is less than
! cwfn%nband (# of bands in total) **AND** is a valence band (so
! its index is <= vwfn%nband), then it is written to iunit_v.
!
    allocate(zc (kpq%ngk(ik), kpq%nspinor*kpq%nspin))
    do i=1,kpq%mnband
      dont_read = (i > cwfn%nband .or. i <= vwfn%ncore_excl)
      if(.not. dont_read) dont_read = i > vwfn%nband+vwfn%ncore_excl+pol%ncrit
      call read_binary_data(26, kpq%ngk(ik), kpq%ngk(ik), kpq%nspin*kpq%nspinor, zc, dont_read = dont_read)
      ! FHJ: the following lines were introduced in r6294 and are supposed to
      ! be a shortcut if we are past the last band of the last k-point. However,
      ! in light of a previous bug (#223), this feature is commented out for now.
      !! FHJ: shortcut if this is past the last band of the last k-point
      !if (dont_read .and. ik==kpq%nrk) exit
      if (.not.dont_read) then
        ! DVF: recall that we redefined the number of valence bands to exclude the
        ! core states. So, we have to subtract ncore_excl right here because ib is
        ! referenced to the full wavefunction file including all the core states.
        i2=i-vwfn%ncore_excl
        call progress_step(prog_info, (ik-1)*(vwfn%nband+pol%ncrit) + i)
        if (peinf%inode == 0) then
          do is = 1, kpq%nspin
            call checknorm('WFNq',i,ik,kpq%ngk(ik),is,kpq%nspinor,zc(:,:))
          enddo
        endif
        if (peinf%doiownv(i2)) then
          iiii=peinf%indexv(i2)+(ik-1)*peinf%nvownactual
          intwfnvq%cg(1:kpq%ngk(ik),iiii,1:kpq%nspin*kpq%nspinor)=zc(1:kpq%ngk(ik),1:kpq%nspinor*kpq%nspin)
        endif
      endif
    enddo
    if(allocated(isort))then;deallocate(isort);endif
    if(allocated(zc))then;deallocate(zc);endif
  enddo ! end loop over k+q points
  call progress_free(prog_info)
 
  return
end subroutine read_wavefunctions
end module input_q_m
