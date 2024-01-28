!==============================================================================
!
! Routines:
!
!
! (1) write_eigenvalues Originally by MLT Last Edited: 4/13/2016 (GKA)
!
! Write the eigenvalues to a text file.
!
! (2) read_eigenvalues Originally by MLT Last Edited: 4/13/2016 (GKA)
!
! Read the eigenvalues from a text file.
!
! (3) write_eigenvalues_noeh Originally by MLT Last Edited: 4/13/2016 (GKA)
!
! Write the eigenvalues without e-h interaction to a text file.
!
! (4) read_eigenvalues_noeh Originally by GKA Last Edited: 4/13/2016 (GKA)
!
! Read the eigenvalues without e-h interaction from a text file.
!
! (5) write_eigenvectors Originally By MLT Last Edited: 6/6/2008 (JRD)
!
! Sohrab Ismail-Beigi: May 4, 2001
!
! This routine is called when we want to write eigenvectors of the BSE
! equations to file (i.e. the eigenvectors coefficients). It will open
! a file called eigenvectors and write to it.
!
! (6) write_vmtxel Originally by MLT Last Edited: 4/13/2016 (GKA)
!
! Write the matrix elements of the momentum operator.
!
! (7) read_vmtxel Originally by MLT Last Edited: 4/13/2016 (GKA)
!
! Read the matrix elements of the momentum operator.
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
module absp_io_m
  use global_m
  use misc_m
  implicit none
  public :: write_eigenvalues, read_eigenvalues, &
            write_eigenvalues_noeh, read_eigenvalues_noeh, &
            write_eigenvectors, &
            write_vmtxel, read_vmtxel
  private
contains
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine write_eigenvalues(xct,flag,neig,vol,evals,cs,dipoles_r,dipoles_l)
  ! Arguments
  type (xctinfo), intent(in) :: xct
  type (flags), intent(in) :: flag
  integer, intent(in) :: neig
  real(DP), intent(in) :: vol
  real(DP), intent(in) :: evals(neig), cs(neig,xct%npol)
  real(DP), intent(in) :: dipoles_r(neig,xct%npol)
  real(DP), intent(in), optional :: dipoles_l(neig,xct%npol)
  ! Local variables
  integer :: ii,ipol
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
  ! ------------------------------------------------------------
 
  if (any(evals(1:neig)<-TOL_Zero).and.xct%tda) then
    write(0,'(a)') "WARNING: There are negative excitation energies."
  end if
  if (peinf%inode.eq.0) then
    do ipol = 1,xct%npol
      if (xct%npol==1) then
        fname = 'eigenvalues.dat'
      else
        fname = 'eigenvalues_'//suffix(ipol)//'.dat'
      endif
      call open_file(unit=14,file=trim(fname),form='formatted',status='replace')
      write(14,'(a,i8)') '# neig  = ', neig
      write(14,'(a,e16.9)') '# vol   = ', vol
      write(14,'(a,2i8)') '# nspin, nspinor = ', xct%nspin, xct%nspinor
      write(14,'(a)',advance='no') '#       eig (eV)   abs(dipole)^2'
      if (flag%krnl .ne. 0) then
        write(14,'(a)',advance='no') '    '
      else ! triplet transition matrix element = 0 if we consider spin overlap
        write(14,'(a)',advance='no') ' mg ' !FHJ: what`s the purpose of this "magn" comment?
      endif
      if (xct%tda) then
        write(14,'(a)') '      dipole'
      else
        write(14,'(a)') '    dipole_l        dipole_r'
      endif
      if (xct%tda) then
        do ii=1,neig
          write(14,'(4e16.8)') evals(ii),cs(ii,ipol),dipoles_r(ii,ipol)
        enddo
      else
        do ii=1,neig
          write(14,'(6e16.8)') evals(ii),cs(ii,ipol),dipoles_l(ii,ipol),dipoles_r(ii,ipol)
        enddo
      endif
      call close_file(14)
    enddo
  endif
 
  return
end subroutine write_eigenvalues
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine read_eigenvalues(xct,neig,vol,evals,cs0,ipol)
  ! Arguments
  type (xctinfo), intent(inout) :: xct
  integer, intent(out) :: neig
  real(DP), intent(out) :: vol
  real(DP), allocatable, intent(out) :: evals(:), cs0(:)
  integer, intent(in) :: ipol
  ! Local variables
  integer :: ii
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
  ! ------------------------------------------------------------
 
  if (peinf%inode.eq.0) then
    if (xct%npol==1) then
      fname = 'eigenvalues.dat'
    else
      fname = 'eigenvalues_'//suffix(ipol)//'.dat'
    endif
    call open_file(unit=14,file=trim(fname),form='formatted',status='old')
    read(14,'(10x,i8)') neig
    read(14,'(10x,e16.9)') vol
    read(14,'(19x,2i8)') xct%nspin, xct%nspinor
    write(6,'(a,i8)') '# neig  = ', neig
    write(6,'(a,e16.9)') '# vol   = ', vol
    write(6,'(a,2i8)') '# nspin, nspinor = ', xct%nspin, xct%nspinor
    read(14,*)
    allocate(cs0 (neig))
    allocate(evals (neig))
    do ii=1,neig
      read(14,*) evals(ii), cs0(ii)
    enddo
    call close_file(14)
  endif
 
  return
end subroutine read_eigenvalues
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine write_eigenvalues_noeh(xct,neig,vol,eqp,s0,ipol)
  ! Arguments
  type (xctinfo), intent(in) :: xct
  integer, intent(in) :: neig
  real(DP), intent(in) :: vol
  type (eqpinfo), intent(in) :: eqp
  real(DP), intent(in) :: s0(xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin)
  integer, intent(in) :: ipol
  ! Local variables
  integer :: ii
  integer :: ic,iv,ik,ikcvs,is
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
  ! ------------------------------------------------------------
 
  if (xct%npol==1) then
    fname = 'eigenvalues_noeh.dat'
  else
    fname = 'eigenvalues_'//suffix(ipol)//'_noeh.dat'
  endif
  call open_file(10,file=trim(fname),form='formatted',status='replace')
  write(10,'(a,i8)') '# neig  = ', neig
  write(10,'(a,e16.9)') '# vol   = ', vol
  write(10,'(a,5i8)') '# nspin, nspinor, nkpt, ncb, nvb = ', &
    xct%nspin, xct%nspinor, xct%nkpt_fi, xct%ncb_fi, xct%nvb_fi
  write(10,'(a)',advance='no') '#   ik    ic    iv    is         ec (eV)         ev (eV)        eig (eV)   abs(dipole)^2'
  write(10,'(a)') '          dipole'
  do ik=1,xct%nkpt_fi
    do ic=1,xct%ncb_fi
      do iv=1,xct%nvb_fi
        do is=1,xct%nspin
          ikcvs = bse_index(ik, ic, iv, is, xct)
          if (xct%qflag.ne.2) then
            write(10,'(4i6,6e16.8)') ik, ic, iv, is, eqp%ecqp(ic,ik,is)*ryd, eqp%evqp(iv,ik,is)*ryd, &
              (eqp%ecqp(ic,ik,is) - eqp%evqp(iv,ik,is))*ryd, abs(s0(ikcvs))**2, s0(ikcvs)
          else
            if (xct%indexq_fi(ik).eq.0 .and. xct%patched_sampling) cycle
            write(10,'(4i6,6e16.8)') ik, ic, iv, is, eqp%ecqp(ic,ik,is)*ryd, eqp%evqp(iv,xct%indexq_fi(ik),is)*ryd, &
              (eqp%ecqp(ic,ik,is) - eqp%evqp(iv,xct%indexq_fi(ik),is))*ryd, abs(s0(ikcvs))**2, s0(ikcvs)
          endif
        enddo
      enddo
    enddo
  enddo
  call close_file(10)
 
  return
end subroutine write_eigenvalues_noeh
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine read_eigenvalues_noeh(xct,neig,vol,eqp,s0,ipol)
  ! Arguments
  type (xctinfo), intent(inout) :: xct
  integer, intent(out) :: neig
  real(DP), intent(out) :: vol
  type (eqpinfo), intent(inout) :: eqp
  real(DP), allocatable, intent(out) :: s0(:)
  integer, intent(in) :: ipol
  ! Local variables
  integer :: ii,ic,iv,ik,ikcvs,is
  real(DP) :: ec, ev, eig, cs, cr, ci
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
  ! ------------------------------------------------------------
 
  if (xct%npol==1) then
    fname = 'eigenvalues_noeh.dat'
  else
    fname = 'eigenvalues_'//suffix(ipol)//'_noeh.dat'
  endif
  call open_file(10,file=trim(fname),form='formatted',status='old')
  read(10,'(10x,i8)') neig
  read(10,'(10x,e16.9)') vol
  read(10,'(35x,5i8)') xct%nspin, xct%nspinor, xct%nkpt_fi, xct%ncb_fi, xct%nvb_fi
  write(6,'(a,i8)') '# neig  = ', neig
  write(6,'(a,e16.9)') '# vol   = ', vol
  write(6,'(a,5i8)') '# nspin, nspinor, nkpt, ncb, nvb = ', &
    xct%nspin, xct%nspinor, xct%nkpt_fi, xct%ncb_fi, xct%nvb_fi
  read(10,*)
  allocate(eqp%evqp (xct%nvb_fi,xct%nkpt_fi,xct%nspin))
  allocate(eqp%ecqp (xct%ncb_fi,xct%nkpt_fi,xct%nspin))
  allocate(eqp%evlda (xct%nvb_fi,xct%nkpt_fi,xct%nspin))
  allocate(eqp%eclda (xct%ncb_fi,xct%nkpt_fi,xct%nspin))
  allocate(s0 (neig))
  do ii=1,neig
    read(10,*) ik, ic, iv, is, ec, ev, eig, cs, cr, ci
    ikcvs = bse_index(ik, ic, iv, is, xct)
    s0(ikcvs) = cmplx(cr,ci,kind=DPC)
    if (xct%qflag.ne.2) then
      eqp%ecqp(ic,ik,is) = ec / ryd
      eqp%evqp(iv,ik,is) = ev / ryd
      eqp%eclda(ic,ik,is) = ec / ryd
      eqp%evlda(iv,ik,is) = ev / ryd
    else
      eqp%ecqp(ic,ik,is) = ec / ryd
      eqp%evqp(iv,xct%indexq_fi(ik),is) = ev / ryd
      eqp%eclda(ic,ik,is) = ec / ryd
      eqp%evlda(iv,xct%indexq_fi(ik),is) = ev / ryd
    endif
  enddo
  call close_file(10)
 
  return
end subroutine read_eigenvalues_noeh
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine write_eigenvectors(xct,kg,ld_evecs,pblock,neig,evals,evecs_r,nwrite,evecs_l)
  type(xctinfo), intent(in) :: xct
  type(grid), intent(in) :: kg
  integer, intent(in) :: ld_evecs, pblock, neig
  real(DP), intent(in) :: evals(neig)
  real(DP), intent(in) :: evecs_r(ld_evecs, pblock)
  integer, intent(inout) :: nwrite
  real(DP), intent(in), optional :: evecs_l(ld_evecs, pblock)
  integer :: ieig,ii,jj,kk,peadd
  real(DP) :: A_r(ld_evecs), A_l(ld_evecs)
  integer :: rank_r, rank_l
  logical :: io_r, io_l, full_bse
  character(len=128) :: fname_r, fname_l
 
! Who can do io?
! FHJ: we do i/o for the right and left eigenvectors using different MPI ranks.
! The hope is to get better load balance on a lustre FS. The optimal solution,
! however, is to use HDF5.
  rank_r = 0
  io_r = peinf%inode==rank_r
  full_bse = present(evecs_l)
  rank_l = peinf%npes-1
  io_l = (peinf%inode==rank_l) .and. full_bse
  if (nwrite.lt.0) nwrite=neig
! Open the file we will write to and write header information
  if (io_r) then
    write(6,*)
    if (full_bse) then
      fname_r = "eigenvectors_r"
      write(6,'(1x,a,i0,a)') 'Writing ',nwrite, ' right eigenvectors to file "'//trim(fname_r)//'"'
    else
      fname_r = "eigenvectors"
      write(6,'(1x,a,i0,a)') 'Writing ',nwrite, ' eigenvectors to file "'//trim(fname_r)//'"'
    endif
    write(6,'(1x,a,i0)') 'Length of each vector: ', ld_evecs
    write(6,*)
    call open_file(unit=200,file=trim(fname_r),form='unformatted',status='replace')
    write(200) xct%nspin
    write(200) xct%nvb_fi
    write(200) xct%ncb_fi
    write(200) xct%nkpt_fi
    write(200) ((kg%f(jj,kk),jj=1,3),kk=1,xct%nkpt_fi)
  endif
  if (io_l) then
    fname_l = "eigenvectors_l"
    write(6,*)
    write(6,'(1x,a,i0,a)') 'Writing ',nwrite, ' left eigenvectors to file "eigenvectors_l"'
    write(6,'(1x,a,i0)') 'Length of each vector: ', ld_evecs
    write(6,*)
    call open_file(unit=201,file=trim(fname_l),form='unformatted',status='replace')
    write(201) xct%nspin
    write(201) xct%nvb_fi
    write(201) xct%ncb_fi
    write(201) xct%nkpt_fi
    write(201) ((kg%f(jj,kk),jj=1,3),kk=1,xct%nkpt_fi)
  endif
! Loop over states to be written
  do ieig=1,nwrite
! Figure out which processor (peadd) and column
! state ieig belongs to
    peadd_loop: do peadd=1,peinf%npes
      do jj=1,peinf%neig(peadd)
        if (peinf%peig(peadd,jj).eq.ieig) then
! Get the coeffs for state ieig into A (on all processors)
          if (peinf%inode==peadd-1) then
            A_r(:) = evecs_r(:,jj)
            if (full_bse) then
              A_l(:) = evecs_l(:,jj)
            endif
          endif
! Write to file
          if (io_r) then
            if (peinf%verb_debug) then
              if (full_bse) then
                write(6,'(1x,a,i0,a,f0.6)') 'Writing right state ',ieig,' energy = ',evals(ieig)
              else
                write(6,'(1x,a,i0,a,f0.6)') 'Writing state ',ieig,' energy = ',evals(ieig)
              endif
            endif
            write(200) evals(ieig)
            write(200) (A_r(ii),ii=1,ld_evecs)
          endif
          if (io_l) then
            if (peinf%verb_debug) then
              write(6,'(1x,a,i0,a,f0.6)') 'Writing left state ',ieig,' energy = ',evals(ieig)
            endif
            write(201) evals(ieig)
            write(201) (A_l(ii),ii=1,ld_evecs)
          endif
        endif
      enddo
    enddo peadd_loop
  enddo ! of loop over states (ieig)
  if (io_r) call close_file(200)
  if (io_l) call close_file(201)
 
  return
end subroutine write_eigenvectors
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine write_vmtxel(xct,flag,nmat,s1)
  ! Arguments
  type (xctinfo), intent(in) :: xct
  type (flags), intent(in) :: flag
  integer, intent(in) :: nmat
  real(DP), intent(in) :: s1(nmat,xct%npol)
  ! Local variables
  integer :: ipol
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
  ! ------------------------------------------------------------
 
  if (peinf%inode.eq.0) then
    write(6,'(1x,a)') 'Writing matrix elements into vmtxel'
    do ipol=1,xct%npol
      if (xct%npol==1) then
        fname = 'vmtxel'
      else
        fname = 'vmtxel_'//suffix(ipol)
      endif
      call open_file(16, file=trim(fname), form='unformatted', status='replace')
      write(16) xct%nkpt_fi,xct%ncb_fi,xct%nvb_fi,xct%nspin,flag%opr
      write(16) s1(:,ipol)
      call close_file(16)
    enddo
  ! If you want this file, you can get it by bsebinasc anyway.
  ! call open_file(17,file='vmtxel.dat',form='formatted',status='replace')
  ! write(17,*) xct%nkpt_fi,xct%ncb_fi,xct%nvb_fi,xct%nspin,flag%opr
  ! write(17,*) (s1(ikcvs),ikcvs=1,nmat)
  ! call close_file(17)
  endif
 
  return
end subroutine write_vmtxel
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine read_vmtxel(xct,flag,nmat,s1)
  ! Arguments
  type (xctinfo), intent(inout) :: xct
  type (flags), intent(in) :: flag
  integer, intent(in) :: nmat
  real(DP), intent(out) :: s1(nmat,xct%npol)
  ! Local variables
  integer :: ii,ipol
  integer :: ic,iv,ik,is
  character(len=128) :: fname
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
  ! ------------------------------------------------------------
 
    if (peinf%inode.eq.0) then
      write(6,'(1x,a)') 'Reading matrix elements from vmtxel'
      do ipol=1,xct%npol
        if (xct%npol==1) then
          fname = 'vmtxel'
        else
          fname = 'vmtxel_'//suffix(ipol)
        endif
        call open_file(16, file=trim(fname), form='unformatted', status='old')
        read(16) ik,ic,iv,is,ii
        if (ik.ne.xct%nkpt_fi.or.ic.ne.xct%ncb_fi.or.iv.ne.xct%nvb_fi &
          .or.is.ne.xct%nspin.or.ii.ne.flag%opr) then
          write(0,'(a,5i6)') 'read  : ', ik,ic,iv,is,ii
          write(0,'(a,5i6)') 'needed: ', xct%nkpt_fi,xct%ncb_fi,xct%nvb_fi,xct%nspin,flag%opr
          call die('parameter mismatch in vmtxel')
        endif
        read(16) s1(:,ipol)
        call close_file(16)
      enddo
    endif
 
  return
end subroutine read_vmtxel
end module absp_io_m
