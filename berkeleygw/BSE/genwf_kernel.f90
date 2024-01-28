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
module genwf_kernel_m
  use blas_m
  use global_m
  use gmap_m
  use input_utils_m
  use misc_m
  use sort_m
  use susymmetries_m
  implicit none
  private
  public :: &
    genwf_kernel
contains
!-----------------------------------------------------------------------
subroutine genwf_kernel(crys,gvec,kg,kgq,syms, &
  wfnc,wfnv,nspin,ik,ic,iv,indexq,xct,intwfnv,intwfnc)
!
! input: crys, gvec, kg, syms types
! nspin (number of spins)
! ik (index of k-point in full BZ)
! ib (index of ck block)
! output: wfnc (conduction states at point ik, block ib)
! wfnv (valence states at point ik)
!
!-----------------------------------------------------------------------
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (grid), intent(in) :: kg
  type (grid), intent(in) :: kgq
  type (symmetry), intent(in) :: syms
  type (wavefunction), intent(out) :: wfnc,wfnv
  integer, intent(in) :: nspin,ik,ic,iv
  integer, intent(in) :: indexq(kg%nf)
  type (int_wavefunction), intent(in) :: intwfnc,intwfnv
  type (xctinfo), intent(in) :: xct
  character :: errmsg*100
  integer :: cnb,cng,cns,cnsp,vng,vns,vnsp
  integer :: ig,ii,jj,kk,info
  integer :: ivown,icown,ikown,invband,incband,ijk
  real(DP) :: xnorm,qk(3)
  integer, save :: ikold=0
  integer, save :: ikold2=0
  integer, save :: vngold=0
  integer, save :: cngold=0
  integer, save :: ifirst=0
  integer, allocatable :: isorti(:)
  integer, allocatable :: ind(:),isort(:)
  integer, allocatable :: indq(:),isortq(:)
  integer, allocatable, save :: ind_old(:),isort_old(:)
  integer, allocatable, save :: indq_old(:),isortq_old(:)
  integer, allocatable, save :: ind_old2(:),isort_old2(:)
  integer, allocatable, save :: indq_old2(:),isortq_old2(:)
  real(DP), allocatable :: ekin(:)
  real(DP), allocatable :: ccg(:,:,:),vcg(:,:,:),ph(:),phq(:)
  real(DP), allocatable, save :: ph_old(:),phq_old(:)
  real(DP), allocatable, save :: ph_old2(:),phq_old2(:)
 
!-----------------------------------------------------------------------
! Deal with the valence wavefunctions
! write(6,*) peinf%inode, 'in genwf',iv,ic,ik
  if (xct%qflag .eq. 0) then
    ivown = peinf%ipev(peinf%inode+1,iv,kgq%indr(indexq(ik)))
    ikown = peinf%ipekq(peinf%inode+1,kgq%indr(indexq(ik)))
  elseif (xct%qflag .eq. 2) then
    ! DYQ: For finite center-of-mass Q, the valence wavefunctions
    ! will be at k+Q, and the conduction band will be at k
    ivown = peinf%ipev(peinf%inode+1,iv,kg%indr(indexq(ik)))
    ikown = peinf%ipekq(peinf%inode+1,kg%indr(indexq(ik)))
  else
    ivown = peinf%ipev(peinf%inode+1,iv,kg%indr(ik))
    ikown = peinf%ipek(peinf%inode+1,kg%indr(ik))
  endif
  vng = intwfnv%ng(ikown)
  vns = intwfnv%nspin
  vnsp = intwfnv%nspinor
  if (xct%ivpar .eq. 1) then
    allocate(vcg (vng,1,vns*vnsp))
  else
    allocate(vcg (vng,xct%nvb_co,vns*vnsp))
  endif
  allocate(indq (vng))
  allocate(phq (vng))
  allocate(isortq (gvec%ng))
  if (ik .ne. ikold .and. ik .ne. ikold2 .and. ifirst .ne. 0) then
    if(allocated(indq_old2))then;deallocate(indq_old2);endif
    if(allocated(phq_old2))then;deallocate(phq_old2);endif
    if(allocated(isortq_old2))then;deallocate(isortq_old2);endif
    allocate(indq_old2 (vngold))
    allocate(phq_old2 (vngold))
    allocate(isortq_old2 (gvec%ng))
    indq_old2 = indq_old
    phq_old2 = phq_old
    isortq_old2 = isortq_old
    if(allocated(indq_old))then;deallocate(indq_old);endif
    if(allocated(phq_old))then;deallocate(phq_old);endif
    if(allocated(isortq_old))then;deallocate(isortq_old);endif
    allocate(indq_old (vng))
    allocate(phq_old (vng))
    allocate(isortq_old (gvec%ng))
  endif
  if (ifirst .eq. 0) then
    allocate(indq_old (vng))
    allocate(phq_old (vng))
    allocate(isortq_old (gvec%ng))
    allocate(indq_old2 (vng))
    allocate(phq_old2 (vng))
    allocate(isortq_old2 (gvec%ng))
  endif
! Initalize parameters in variable wfnv
  wfnv%ng=vng
  wfnv%nband=1
  wfnv%nspin=vns
  wfnv%nspinor=vnsp
  if (intwfnv%nspin.ne.nspin) then
    write(errmsg,'(2a,2i2)') 'spin number mismatch: ', nspin, vns
    call die(errmsg, only_root_writes = .true.)
  endif
  if (xct%ivpar .eq. 1) then
    allocate(wfnv%cg (wfnv%ng,1,wfnv%nspin*wfnv%nspinor))
    vcg(1:vng,1,:) = intwfnv%cg(1:vng,ivown,:)
  else
    allocate(wfnv%cg (wfnv%ng,xct%nvb_co,wfnv%nspin*wfnv%nspinor))
    vcg(1:vng,:,:) = intwfnv%cg(1:vng,ivown:ivown+xct%nvb_co-1,:)
  endif
  allocate(wfnv%isort (gvec%ng))
  isortq(:) = intwfnv%isort(:,ikown)
! Compute inverse index array of Fourier components around rk-kpoint
  if (ik .ne. ikold .and. ik .ne. ikold2) then
    allocate(isorti (gvec%ng))
    isorti(:)=0
    do ii=1,wfnv%ng
      isorti(isortq(ii))=ii
    enddo
  endif
! Compute index array of Fourier components around fk-kpoint
  if (ik .ne. ikold .and. ik .ne. ikold2) then
    if (peinf%verb_debug) then
      write(6,*) 'ikneikold',peinf%inode,ik,ikold
    endif
    allocate(ekin (gvec%ng))
    if (xct%qflag .eq. 0) then
      qk(1:3)=kgq%f(1:3,indexq(ik))
    elseif (xct%qflag .eq. 2) then
      qk(1:3)=kg%f(1:3,indexq(ik))
    else
      qk(1:3)=kg%f(1:3,ik)
    endif
    call kinetic_energies(gvec, crys%bdot, ekin, qvec = qk)
    call sortrx(gvec%ng,ekin,isortq)
    if(allocated(ekin))then;deallocate(ekin);endif
    isortq_old(:)=isortq(:)
  elseif (ik .eq. ikold) then
    isortq(:)=isortq_old(:)
  else
    isortq(:)=isortq_old2(:)
  endif
! Find ind and ph relating wavefunctions in fk to rk-kpoint
  if (ik .ne. ikold .and. ik .ne. ikold2) then
    indq=0
    phq=0.0d0
    if (xct%qflag .eq. 0) then
      call gmap(gvec,syms,wfnv%ng,kgq%itran(indexq(ik)), &
        kgq%kg0(:,indexq(ik)),isortq,isorti,indq,phq,.true.)
    elseif (xct%qflag .eq. 2) then
      call gmap(gvec,syms,wfnv%ng,kg%itran(indexq(ik)), &
        kgq%kg0(:,ik),isortq,isorti,indq,phq,.true.)
    else
      call gmap(gvec,syms,wfnv%ng,kg%itran(ik), &
        kg%kg0(:,ik),isortq,isorti,indq,phq,.true.)
    endif
    if(allocated(isorti))then;deallocate(isorti);endif
    indq_old(:)=indq(:)
    phq_old(:)=phq(:)
  else if (ik .eq. ikold) then
    indq(:)=indq_old(:)
    phq(:)=phq_old(:)
  else
    indq(:)=indq_old2(:)
    phq(:)=phq_old2(:)
  endif
! Compute and renormalize valence wavefunctions in fk-kpoint
  if (xct%ivpar .eq. 1) invband = 1
  if (xct%ivpar .eq. 0) invband = xct%nvb_co
  do ijk = 1, invband
    do kk=1,wfnv%nspin
      do jj=1,wfnv%nspinor
        do ii=1,wfnv%ng
          if (indq(ii) .gt. 0) then
            wfnv%cg(ii,ijk,kk*jj)=phq(ii)*vcg(indq(ii),ijk,kk*jj)
          else
            wfnv%cg(ii,ijk,kk*jj) = 0.0d0
          endif
        enddo !ii
      enddo ! jj
      call checknorm('wfnv%cg',ijk,ik,wfnv%ng,kk,wfnv%nspinor,wfnv%cg(1:wfnv%ng,ijk,:))
    enddo ! kk
  enddo ! ijk
  vcg=wfnv%cg
  ! In spinor case, we must rotate spinors according to spinor rotation matrix umtrx
  wfnv%cg=vcg
  wfnv%isort=isortq
  if(allocated(vcg))then;deallocate(vcg);endif
  if(allocated(indq))then;deallocate(indq);endif
  if(allocated(phq))then;deallocate(phq);endif
  if(allocated(isortq))then;deallocate(isortq);endif
!-----------------------------------------------------------------------
! Deal with the conduction wavefunctions
  icown = peinf%ipec(peinf%inode+1,ic,kg%indr(ik))
  ikown = peinf%ipek(peinf%inode+1,kg%indr(ik))
  cng = intwfnc%ng(ikown)
  cnb = 1
  cns = intwfnc%nspin
  cnsp = intwfnc%nspinor
  if (xct%icpar .eq. 1) then
    allocate(ccg (cng,1,cns*cnsp))
  else
    allocate(ccg (cng,xct%ncb_co,cns*cnsp))
  endif
  allocate(ind (cng))
  allocate(ph (cng))
  allocate(isort (gvec%ng))
  if (ik .ne. ikold .and. ik .ne. ikold2 .and. ifirst .ne. 0) then
    if(allocated(ind_old2))then;deallocate(ind_old2);endif
    if(allocated(ph_old2))then;deallocate(ph_old2);endif
    if(allocated(isort_old2))then;deallocate(isort_old2);endif
    allocate(ind_old2 (cngold))
    allocate(ph_old2 (cngold))
    allocate(isort_old2 (gvec%ng))
    ind_old2 = ind_old
    ph_old2 = ph_old
    isort_old2 = isort_old
    if(allocated(ind_old))then;deallocate(ind_old);endif
    if(allocated(ph_old))then;deallocate(ph_old);endif
    if(allocated(isort_old))then;deallocate(isort_old);endif
    allocate(ind_old (cng))
    allocate(ph_old (cng))
    allocate(isort_old (gvec%ng))
  endif
  if (ifirst .eq. 0) then
    allocate(ind_old (cng))
    allocate(ph_old (cng))
    allocate(isort_old (gvec%ng))
    allocate(ind_old2 (cng))
    allocate(ph_old2 (cng))
    allocate(isort_old2 (gvec%ng))
  endif
  wfnc%ng=cng
  wfnc%nband=1
  wfnc%nspin=cns
  wfnc%nspinor=cnsp
  if (cns.ne.nspin) then
    write(errmsg,'(2a,2i2)') 'spin number mismatch: ', nspin, cns
    call die(errmsg, only_root_writes = .true.)
  endif
  allocate(wfnc%isort (gvec%ng))
  isort(:)=intwfnc%isort(:,ikown)
  if (xct%icpar .eq. 1) then
    allocate(wfnc%cg (wfnc%ng,1,wfnc%nspin*wfnc%nspinor))
    ccg(1:cng,1,:) = intwfnc%cg(1:cng,icown,:)
  else
    allocate(wfnc%cg (wfnc%ng,xct%ncb_co,wfnc%nspin*wfnc%nspinor))
    ccg(1:cng,:,:) = intwfnc%cg(1:cng,icown:icown+xct%ncb_co-1,:)
  endif
! JRD: Below is now necessary because kg might be different from kgq
! Compute inverse index array of Fourier components around rk-kpoint
  if (ik .ne. ikold .and. ik .ne. ikold2) then
    allocate(isorti (gvec%ng))
    isorti(:)=0
    do ii=1,wfnc%ng
      isorti(isort(ii))=ii
    enddo
  endif
! Compute index array of Fourier components around fk-kpoint
  if (ik .ne. ikold .and. ik .ne. ikold2) then
    allocate(ekin (gvec%ng))
    call kinetic_energies(gvec, crys%bdot, ekin, qvec = kg%f(:, ik))
    call sortrx(gvec%ng,ekin,isort)
    if(allocated(ekin))then;deallocate(ekin);endif
    isort_old(:) = isort(:)
  else if (ik .eq. ikold) then
    isort(:)=isort_old(:)
  else
    isort(:)=isort_old2(:)
  endif
! Find ind and ph relating wavefunctions in fk to rk-kpoint
  if (ik .ne. ikold .and. ik .ne. ikold2) then
    ind=0
    ph=0.0d0
    call gmap(gvec,syms,wfnc%ng,kg%itran(ik), &
      kg%kg0(:,ik),isort,isorti,ind,ph,.true.)
    if(allocated(isorti))then;deallocate(isorti);endif
    ind_old(:)=ind(:)
    ph_old(:)=ph(:)
  else if (ik .eq. ikold) then
    ind(:)=ind_old(:)
    ph(:)=ph_old(:)
  else
    ind(:)=ind_old2(:)
    ph(:)=ph_old2(:)
  endif
! Compute and renormalize conduction wavefunctions
  if (xct%icpar .eq. 1) incband =1
  if (xct%icpar .eq. 0) incband = xct%ncb_co
  do ijk =1, incband
    do kk=1,wfnc%nspin
      do jj=1,wfnc%nspinor
        do ii=1,wfnc%ng
          if (ind(ii) .gt. 0) then
            wfnc%cg(ii,ijk,kk*jj)=ph(ii)*ccg(ind(ii),ijk,kk*jj)
          else
            wfnc%cg(ii,ijk,kk*jj)=0.0d0
          endif
        enddo
      enddo ! jj
      call checknorm('wfnc%cg',ijk,ik,wfnc%ng,kk,wfnc%nspinor,wfnc%cg(1:wfnc%ng,ijk,:))
    enddo ! kk
  enddo ! ijk
  ccg=wfnc%cg
  ! In spinor case, we must rotate spinors according to spinor rotation matrix umtrx
  wfnc%cg=ccg
  wfnc%isort=isort
  if(allocated(ccg))then;deallocate(ccg);endif
  if(allocated(ind))then;deallocate(ind);endif
  if(allocated(ph))then;deallocate(ph);endif
  if(allocated(isort))then;deallocate(isort);endif
  if (ik .ne. ikold .and. ik .ne. ikold2) then
    ikold2=ikold
    ikold=ik
    cngold=cng
    vngold=vng
  endif
  ifirst=-1
 
  return
end subroutine genwf_kernel
end module genwf_kernel_m
