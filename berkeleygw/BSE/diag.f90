!===============================================================================
!
! Routines:
!
! (1) diag main Originally by MLT Last Edited: 5/12/2008 (JRD)
!
! For more details see the README_absorption file.
!
! Calculates the real and imaginary parts of the macroscopic dielectric
! function starting from interaction matrix elements calculated by
! the Kernel code. It uses interpolation in the matrix elements and
! direct diagonalization of the Bethe-Salpeter equation.
! Spin-polarized case implemented.
!
! For more details, see:
! Rohlfing & Louie, PRB 62:(8), 4927 (2000)
! G. Strinati, Rivista del Nuovo Cimento 11:(12), 1 (1988)
!
! Please report bugs to: jdeslip@civet.berkeley.edu
!
!================================================================================
module diag_m
  use absp0_m
  use absp_io_m
  use absp_lanczos_m
  use absp_m
  use diagonalize_m
  use evecs_m
  use fullbz_m
  use genwf_m
  use global_m
  use input_fi_m
  use input_q_m
  use intkernel_m
  use intwfn_m
  use misc_m
  use timing_m, only: timing => bse_timing
  use vmtxel_m
  use wfn_rho_vxc_io_m, only: init_mf_header_from_types
  implicit none
  private
  public :: &
    diag
contains
subroutine diag(eqp,xct,flag,neig,nmax)
  type (eqpinfo), intent(inout) :: eqp
  type (xctinfo), intent(inout) :: xct
  type (flags), intent(inout) :: flag
  integer, intent(inout) :: neig
  integer, intent(in) :: nmax
  type (crystal) :: crys
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (grid) :: kg_fi, kgq_fi,kg_co,kgq_co
  type (kpoints) :: kp_fi, kpq_fi
  type (wavefunction) :: wfnc_fi
  type (wavefunction) :: wfnvq_fi
  type (work_genwf) :: work, workq
  type (int_wavefunction) :: intwfnc
  type (int_wavefunction) :: intwfnv
  type (evecs_t) :: evecs
  type (vmtxel_t) :: dip
  type (mf_header_t) :: mf_header
  character :: tmpstr*128,filename*20
  integer :: ii,ipol,jj,ncount,nmat,ncvs_fi,pblock
  integer :: ikb, icb, ivb
  integer :: ik,ikq,ikt,iblock,ikcvs,ikcvsd,ic,iv,is
  integer :: iunit_c,iunit_v
  integer :: version
  real(DP) :: vol,omega_plasma,en_diff_max
  real(DP) :: tsec(2),tmin(2),tmax(2)
  character*16, allocatable :: routnam(:)
  integer, allocatable :: routsrt(:)
  integer, allocatable :: fi2co_wfn(:,:),indexq_fi(:)
  real(DP), allocatable :: evals(:), kco(:,:), cs(:,:), cs0(:), rdum(:,:)
  real(DP), allocatable :: dcc(:,:,:,:,:),dvv(:,:,:,:,:),s0(:), &
                         hqpcc(:,:,:), hqpvv(:,:,:), rdum2(:,:)
  real(DP), allocatable :: dipoles_l(:,:), dipoles_r(:,:), cs_full(:,:)
  !> (kcvs, k`c`v`s`), "A" block of BSE Hamiltonian
  real(DP), allocatable :: hbse_a(:,:)
  !> (kcvs, k`c`v`s`), "B" block of BSE Hamiltonian, only if tda=.false.
  real(DP), allocatable :: hbse_b(:,:)
  real(DP), allocatable :: intp_coefs(:,:)
  character(len=2) :: suffix(3) = (/'b1', 'b2', 'b3'/)
  ! DYQ: Variables used in finite Q
  integer, allocatable :: kg0_temp(:,:)
  integer :: umk
  real(DP) :: kq(3),qq(3)
  real :: delta
  !DYQ: Variables used for clustered subsampling
  integer ::ik_sub,nsub_files,nk_sub
  type(grid) :: kg_sub_co
  real(DP), allocatable :: dvv_sub(:,:,:,:,:,:),dcc_sub(:,:,:,:,:,:)
  integer,allocatable :: closepts_sub(:,:)
 
! JRD: A Check for Finite Q
! if (peinf%inode .eq. 0) then
! write(6,*) 'nkpt_co = ', xct%nkpt_co
! endif
  if(flag%vm.eq.2) then
    if (peinf%inode.eq.0) then
      write(0,*) 'WARNING: read_eps2_moments not supported in this diagonalization code. Ignoring keyword.'
    endif
    flag%vm=0
  endif
!--------------------------
! If eigenvalues.dat is available, read them and go straight to
! calculation of the absorption spectrum
  if (flag%spec.eq.1) then
    if (peinf%inode .eq. 0) then
      omega_plasma = 0.d0
      write(6,*) 'Create absorption_noeh.dat from eigenvalues_noeh.dat'
      do ipol=1,xct%npol
        call read_eigenvalues_noeh(xct,neig,vol,eqp,s0,ipol)
        call absp0(eqp,xct,s0,vol,omega_plasma,flag,ipol)
        if(associated(eqp%evqp))then;deallocate(eqp%evqp);nullify(eqp%evqp);endif
        if(associated(eqp%ecqp))then;deallocate(eqp%ecqp);nullify(eqp%ecqp);endif
        if(associated(eqp%evlda))then;deallocate(eqp%evlda);nullify(eqp%evlda);endif
        if(associated(eqp%eclda))then;deallocate(eqp%eclda);nullify(eqp%eclda);endif
        if(allocated(s0))then;deallocate(s0);endif
      enddo
      if (xct%iabsorp0 .eq. 0) then
        write(6,*) 'Create absorption_eh.dat from eigenvalues.dat'
        do ipol=1,xct%npol
          call read_eigenvalues(xct,neig,vol,evals,cs0,ipol)
          call absp(xct,neig,cs0,evals,vol,omega_plasma,flag,ipol)
          if(allocated(evals))then;deallocate(evals);endif
          if(allocated(cs0))then;deallocate(cs0);endif
        enddo
      endif
    endif
    call diag_end()
   
    return
  endif
!--------------------------
! Read wavefunctions on the fine grid
  call logit('Calling input')
  call timing%start(timing%input)
  call input_fi(crys,gvec,kg_fi,kp_fi,syms,eqp,xct,flag, &
    omega_plasma,.true.,intwfnc)
! If there is no specified number of eigenvalues, calculate
! all eigenvalues
  nmat = xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin
  if (xct%algo==BSE_ALGO_DIAG .or. xct%algo==BSE_ALGO_DIAG_PRIMME) then
    if (xct%tda) then
      if (neig==0) neig = nmat
    else
      if (peinf%inode==0.and.neig/=0) then
        write(0,'(/,a,/)') 'WARNING: BSE calculations ignore the `number_eigenvalues` flag.'
      endif
      ! FHJ: generic solver
      neig = 2*nmat
    endif
    if ((neig<=0).or.(neig>2*nmat).or.(neig>nmat.and.xct%tda)&
    ) then
      write(tmpstr,'(a,2i6)') 'Incomprehensible request of eigenvalues : ', neig, nmat
      call die(tmpstr, only_root_writes = .true.)
    endif
  else
    neig = nmat
  endif
  vol = xct%nktotal*crys%celvol
  if (peinf%inode.eq.0) then
    write(6,'(/1x,a)') 'More job parameters:'
    write(6,'(1x,a,es9.3e2)') '- Crystal volume (bohr): ', vol
    write(6,'(1x,a,f0.3)') '- Broadening (eV): ', xct%eta
    write(6,'(1x,a,i0)') '- Number of valence bands: ', xct%nvb_fi
    write(6,'(1x,a,i0)') '- Number of cond. bands: ', xct%ncb_fi
    write(6,'(1x,a,i0)') '- Number of spins: ', xct%nspin
    write(6,'(1x,a,i0)') '- Number of eigenvalues to be computed: ', neig
    write(6,'()')
  endif
  call timing%stop(timing%input)
  allocate(indexq_fi (xct%nkpt_fi))
  allocate(xct%indexq_fi (xct%nkpt_fi))
  if (flag%vm.ne.1.or. .not. flag%read_dtmat) then
    call timing%start(timing%input_q)
    call logit('Calling input_q')
    call input_q(kp_fi,crys,gvec,kg_fi,kgq_fi,kpq_fi,syms,xct,indexq_fi,eqp,flag,intwfnv)
    call timing%stop(timing%input_q)
  endif
! JRD: Don`t do this next section if only want absorption_noeh.dat
  if (xct%iabsorp0 .eq. 0) then
!------------------------------
! Calculate the transformation matrices from coarse grid wavefunctions
! FHJ: These are the final transformation coefs that will be used to interpolate
! the kernel. However, we might use an unrestricted version of dvv/dcc to
! interpolate eqp if xct%unrestricted_transf==.true..
    allocate(dvv (xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
    allocate(dcc (xct%ncb_fi,xct%n2b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
    allocate(kco (3,xct%nkpt_co))
    allocate(fi2co_wfn (xct%npts_intp_kernel,xct%nkpt_fi))
    allocate(intp_coefs (xct%npts_intp_kernel, xct%nkpt_fi))
    call logit('Calling intwfn')
    call timing%start(timing%intwfn)
    call intwfn(kp_fi,crys,syms,xct,flag,gvec,kg_fi,kgq_fi,kg_co,kgq_co,dcc,dvv,&
      kco,fi2co_wfn,indexq_fi,eqp,intwfnv,intwfnc,intp_coefs)
    call timing%stop(timing%intwfn)
    if (xct%subsample_line) then
      call logit('Calling intwfn_sub')
      call open_file(unit=77, file='subsample.inp',form='formatted',status='old')
      read(77,*) nk_sub
      read(77,*) nsub_files
      if (nsub_files.ne.xct%nkpt_co) call die("Number of subsampled points must be the same as the number of coarse points",&
        only_root_writes=.true.)
      kg_sub_co%nr = nsub_files
      allocate(kg_sub_co%r (3,xct%nkpt_co))
      do ik_sub=1,xct%nkpt_co
        read(77,*) (kg_sub_co%r(ii,ik_sub),ii=1,3)
      enddo
      call close_file(77)
      allocate(dcc_sub (xct%ncb_fi, xct%n2b_co, xct%nspin, xct%nkpt_fi,xct%idimensions+1, nk_sub))
      allocate(dvv_sub (xct%nvb_fi, xct%n1b_co, xct%nspin, xct%nkpt_fi,xct%idimensions+1, nk_sub))
      allocate(closepts_sub (4,xct%nkpt_fi))
      call intwfn_sub(kp_fi,crys,syms,xct,flag,gvec,kg_fi,kgq_fi, &
        dcc_sub,dvv_sub,kg_sub_co,kg_co,indexq_fi,nk_sub,intwfnv,intwfnc,closepts_sub)
    endif
  endif
  if(associated(xct%ifmax))then;deallocate(xct%ifmax);nullify(xct%ifmax);endif
  if ((flag%vm.ne.1.or. .not. flag%read_dtmat) .and..not. xct%no_mtxel) then
    ! otherwise, we did not call input_q to allocate it
    if(associated(xct%ifmaxq))then;deallocate(xct%ifmaxq);nullify(xct%ifmaxq);endif
  endif
!------------ Calculate the velocity (or momentum) matrix elements -------------
! Each PE calculates a small number of them. At the end, share data
!
! If flag%vm.eq.1, skip this part and just read the matrix elements
! from "vmtxel".
!
! peinf%block_sz = size of a distributed column in hbse_a
  call logit('Calculating v/p matrix elememts')
  ncvs_fi = xct%ncb_fi*xct%nvb_fi*xct%nspin
  if (xct%ipar .eq. 1) then
    peinf%block_sz = ncvs_fi
  else if (xct%ipar .eq. 2) then
    peinf%block_sz = xct%nvb_fi*xct%nspin
  else
    peinf%block_sz = xct%nspin
  endif
  nmat = xct%nkpt_fi*ncvs_fi
  ! Initialize dipole operator and allocate dipole matrix elements
  call dip%init_from_xctinfo(xct, opr=flag%opr)
  call dip%alloc()
  ! Copy list of k-points
  do ik=1, dip%nk
    dip%kpt(:,ik) = kg_fi%f(:,ik)
  end do
  ! FIXME: It was decided to hold on the use of hdf5 for vmtxel
  ! until the next release, so as not to break backward compatibility
  ! Since using hdf5 for vmtxel is highly desirable, let us not
  ! forget to reactivate it as soon as possible.
  dip%use_hdf5 = .false.
  call timing%start(timing%vmtxel)
  if (.not.xct%no_mtxel) then
  if (flag%vm.eq.0) then
    do ikt=1, peinf%ikt(peinf%inode+1)
      ik = peinf%ik(peinf%inode+1,ikt)
      if (xct%qflag .eq. 1) then
        ikq = indexq_fi(ik)
      else
        ! genwf will retrieve the valence bands at k+q+Q
        ikq = xct%indexq_fi(ik)
        if (ikq.eq.0 .and. xct%patched_sampling) then
          if (peinf%inode.eq.0) write(6,*) "Skipping genwf for ik,ikq:",ik,ikq
          cycle
        endif
      endif
      call genwf(crys,gvec,kg_fi,syms,wfnc_fi,ik,ik,xct%nspin,xct%ncb_fi,&
                 work,intwfnc,xct%iwriteint,is_cond=.true.)
      call genwf(crys,gvec,kgq_fi,syms,wfnvq_fi,ik,ikq,xct%nspin,xct%nvb_fi,&
                 workq,intwfnv,xct%iwriteint,is_cond=.false.)
      call dip%compute_ik_vmtxel(ik, wfnc_fi, wfnvq_fi, gvec, xct%qshift, crys, eqp)
      if(associated(wfnc_fi%cg))then;deallocate(wfnc_fi%cg);nullify(wfnc_fi%cg);endif
      if(associated(wfnc_fi%isort))then;deallocate(wfnc_fi%isort);nullify(wfnc_fi%isort);endif
      if(associated(wfnvq_fi%cg))then;deallocate(wfnvq_fi%cg);nullify(wfnvq_fi%cg);endif
      if(associated(wfnvq_fi%isort))then;deallocate(wfnvq_fi%isort);nullify(wfnvq_fi%isort);endif
    enddo
    ! typedefs initializes all of these ikolds to 0
    if(work%ikold.ne.0) then
      if(associated(work%cg))then;deallocate(work%cg);nullify(work%cg);endif
      if(associated(work%ph))then;deallocate(work%ph);nullify(work%ph);endif
      if(associated(work%ind))then;deallocate(work%ind);nullify(work%ind);endif
      if(associated(work%isort))then;deallocate(work%isort);nullify(work%isort);endif
    endif
    if(workq%ikold.ne.0) then
      if(associated(workq%cg))then;deallocate(workq%cg);nullify(workq%cg);endif
      if(associated(workq%ph))then;deallocate(workq%ph);nullify(workq%ph);endif
      if(associated(workq%ind))then;deallocate(workq%ind);nullify(workq%ind);endif
      if(associated(workq%isort))then;deallocate(workq%isort);nullify(workq%isort);endif
    endif
    ! Share matrix elements
    call dip%reduce()
    ! Write file
    call dip%write_vmtxel()
  else
    ! Read dipole matrix elements from file if they were already computed
    call dip%read_vmtxel()
  endif
  if (flag%vm.ne.1.or. .not. flag%read_dtmat) then
    call dealloc_grid(kgq_fi)
  endif
  if (flag%vm == 0 .and. .not. flag%read_dtmat) then
    if(associated(intwfnc%cgk))then;deallocate(intwfnc%cgk);nullify(intwfnc%cgk);endif
    if(associated(intwfnv%cgk))then;deallocate(intwfnv%cgk);nullify(intwfnv%cgk);endif
    if(associated(intwfnc%isort))then;deallocate(intwfnc%isort);nullify(intwfnc%isort);endif
    if(associated(intwfnv%isort))then;deallocate(intwfnv%isort);nullify(intwfnv%isort);endif
  endif
  else ! xct%no_mtxel
    xct%npol = 1
    xct%pol = 1.0
  endif
  call timing%stop(timing%vmtxel)
! End Calculating Matrix Elements
!-------------------------------------------------------------------------------
!----------------------------
! Calculate the non-interacting spectrum. Only one PE works
  call timing%start(timing%absorp0)
  call logit('Calling absp0')
  if (peinf%inode.eq.0) then
    do ipol=1,xct%npol
      call absp0(eqp,xct,dip%s1(:,ipol),vol,omega_plasma,flag,ipol)
    enddo
  endif
  call timing%stop(timing%absorp0)
  if(associated(eqp%eclda))then;deallocate(eqp%eclda);nullify(eqp%eclda);endif
  if(associated(eqp%evlda))then;deallocate(eqp%evlda);nullify(eqp%evlda);endif
  if(allocated(indexq_fi))then;deallocate(indexq_fi);endif
! JRD If we only want absorp0 - we finish here
  if (xct%iabsorp0 .eq. 1) then
    call diag_end()
   
    return
  endif
!------------ Build Hamiltonian Matrix --------------------------------------------
! Build Hamiltonian matrix. Diagonal part comes from the "non-interacting"
! quasiparticle Hamiltonians. If the QP Greens function is diagonal,
! then these are just the quasiparticle energy differences on the
! diagonal. The more general case is:
!
! <cvk|H0|c'v'k'> = delta(k,k') *
! [ <c|hqp|c'>*delta(v',v) - delta(c,c')*<v'|hqp|v> ]
!
! The rest of the Hamiltonian, which is the electron-hole interaction,
! comes from interpolation further below.
  call logit('Building non-interacting Hamiltonian')
  allocate(hbse_a (xct%nkpt_fi*ncvs_fi, peinf%nblocks*peinf%block_sz))
  hbse_a(:,:) = 0.0d0
  if (.not.xct%tda) then
    allocate(hbse_b (xct%nkpt_fi*ncvs_fi, peinf%nblocks*peinf%block_sz))
    hbse_b(:,:) = 0.0d0
  endif
! iblock loop. This loop is over proc owned k if ipar = 1, (k,c) if
! ipar = 2 and (k,c,v) if ipar = 3
  en_diff_max = 0d0
  do iblock=1,peinf%ibt(peinf%inode+1)
    ik=peinf%ikb(iblock)
    if (ik .eq. 0) then
      write(0,*) "Illegal value for ik",peinf%inode, iblock, ik
      call die("internal error in diag, ik = 0")
    endif
! Build <c|hqp|c'> and <v|hqp|v'> for this kpoint
    allocate(hqpcc (xct%ncb_fi,xct%ncb_fi,xct%nspin))
    allocate(hqpvv (xct%nvb_fi,xct%nvb_fi,xct%nspin))
    hqpcc = 0.0d0
    hqpvv = 0.0d0
! Put QP energies on diagonals of hqpcc and hqpvv to start
    do is=1,xct%nspin
      do ic=1,xct%ncb_fi
        ! Skip if k+Q falls outside the patch
        if (xct%indexq_fi(ik).eq.0 .and. xct%patched_sampling) cycle
        hqpcc(ic,ic,is) = eqp%ecqp(ic,ik,is)
      enddo
      do iv=1,xct%nvb_fi
        if (xct%qflag.ne.2) then
          !DYQ: for qflag.eq.0, eqp%evqp(iv,ik,is) corresponds to the k-point kgq_fi%f(xct%indexq_fi(ik))
          hqpvv(iv,iv,is) = eqp%evqp(iv,ik,is)
        else
          ! Skip if k+Q falls outside the patch
          if (xct%indexq_fi(ik).eq.0 .and. xct%patched_sampling) cycle
          hqpvv(iv,iv,is) = eqp%evqp(iv,xct%indexq_fi(ik),is)
        endif
      enddo
    enddo
! Read possible offdiagonal QP elements from "hqp.<ik>" file
! if it exists. JRD: This is broken for now. Someone should fix
! it in the future if they want to use it
    !if (ik.lt.10) then
    ! write(tmpstr,'(a,i1)') 'hqp.',ik
    !else if (ik.lt.100) then
    ! write(tmpstr,'(a,i2)') 'hqp.',ik
    !else if (ik.lt.1000) then
    ! write(tmpstr,'(a,i3)') 'hqp.',ik
    !else if (ik.lt.100000) then
    ! write(tmpstr,'(a,i5)') 'hqp.',ik
        !else
    ! write(0,*) 'too many kpoints for reading hqp'
    !endif
    !call open_file(9,file=tmpstr,form='formatted',status='old')
    !if (is.eq.0) then
    ! if (peinf%inode.eq.0) then
    ! write(6,*) 'Reading offdiagonal hqp from file ',tmpstr
    ! write(6,*) 'All values in eV'
    ! endif
    ! do
    ! read(9,*,end=999) nocc,ii,jj,x,y
    ! if ii and jj both refer to valence, states, put
    ! matrix element into hqpvv
    ! if ((ii<=nocc).and.(ii>nocc-xct%nvb_fi).and. &
    ! (jj<=nocc).and.(jj>nocc-xct%nvb_fi)) then
    ! if (peinf%inode.eq.0) write(6,'(a,2i5,2f20.10)') ' hqp(v,vp) = ',ii,jj,x,y
    ! ii=nocc-ii+1
    ! jj=nocc-jj+1
    ! is = 1
    ! hqpvv(ii,jj,is) = x/ryd
    ! else if ((ii>nocc).and.(ii<=nocc+xct%ncb_fi).and. &
    ! (jj>nocc).and.(jj<=nocc+xct%ncb_fi)) then
    ! if (peinf%inode.eq.0) write(6,'(a,2i5,2f20.10)') ' hqp(c,cp) = ',ii,jj,x,y
    ! ii=ii-nocc
    ! jj=jj-nocc
    ! is = 1
    ! hqpcc(ii,jj,is) = x/ryd
    ! endif
    ! enddo
    !999 call close_file(9)
    ! write(6,*)
    !endif ! if hqp.<ik> was found
    ! Now build hamiltonian from hqcc and hqvv
    ! Consider only diagonal elements in k,v,c
    ! FHJ: Note: here, iv and ic are dummy indices, and the actual band indices
    ! are ivb/icb. When should we use the dummy or real band indices?
    ! - Use the dummy indices iv/ic to index a column of hbse_a (which is distributed)
    ! - Use the real indices ivb/icb to index a row of hbse_a (which is not distributed)
    ikb = ik
    do is=1,xct%nspin
      do iv=1,peinf%nv_block !1 for ipar==2
        do ic=1,peinf%nc_block !1 for ipar==2 or ipar==3
          if (xct%ipar==1) then
            ivb=iv
            icb=ic
          else if (xct%ipar==2) then
            icb=peinf%icb(iblock)
            ivb=iv
          else if (xct%ipar==3) then
            ivb=peinf%ivb(iblock)
            icb=peinf%icb(iblock)
          endif
          ikcvs = bse_index(ikb, icb, ivb, is, xct)
          ikcvsd = bse_index(iblock, ic, iv, is, xct, ncband=peinf%nc_block, nvband=peinf%nv_block)
          en_diff_max = max(en_diff_max, dble(hqpcc(icb,icb,is) - hqpvv(ivb,ivb,is)))
          hbse_a(ikcvs,ikcvsd) = hqpcc(icb,icb,is) - hqpvv(ivb,ivb,is)
        enddo
      enddo
    enddo
    if(allocated(hqpcc))then;deallocate(hqpcc);endif
    if(allocated(hqpvv))then;deallocate(hqpvv);endif
  enddo ! loop on k-points on this processor
!----------------------------
! Define the mapping of eigenvectors: the ii-th column of the matrix
! evecs%u_r(:,:) stored in PE #ipe will contain the eigenvector of order
! peinf%peig(ipe,ii). The total number of eigenvectors stored in
! each processor is given by peinf%neig(1:peinf%npes).
! pblock >= maxval(peinf%neig(1:peinf%npes))
  ! FHJ: Note: pblock gets ~doubled automatically in full BSE calculations.
  ! In BLACS terms, the following lines would set up a 1d block-column
  ! distribution for the eigenvectors with:
  ! M=(2*)nmat, N=neig, MB=M, NB=peinf%block_sz, LLD=M
  ! We should really get rid of this manual distribution and use BLACS.
  pblock = neig/(peinf%npes*peinf%block_sz)
  if (pblock*peinf%npes*peinf%block_sz.lt.neig) pblock = pblock + 1
  pblock = pblock*peinf%block_sz
  allocate(peinf%neig (peinf%npes))
  allocate(peinf%peig (peinf%npes,pblock))
  peinf%neig=0
  peinf%peig=0
  ii=1
  do jj=1,neig
    if (ii.eq.peinf%npes+1) ii=1
    peinf%neig(ii)=peinf%neig(ii)+1
    peinf%peig(ii,peinf%neig(ii))=jj
    if (mod(jj,peinf%block_sz).eq.0) ii=ii+1
  enddo
!-----------------------------
! Interpolation scheme in the Kernel
  call logit('Calling intkernel')
  call timing%start(timing%intkernel)
  if (xct%tda) then
    if (xct%subsample_line) then
      call intkernel(crys,kg_fi,kp_fi,syms,xct,hbse_a,dcc,dvv,kco,fi2co_wfn,flag,gvec,intp_coefs,&
        dcc_sub=dcc_sub,dvv_sub=dvv_sub,closepts_sub=closepts_sub)
    else
      call intkernel(crys,kg_fi,kp_fi,syms,xct,hbse_a,dcc,dvv,kco,fi2co_wfn,flag,gvec,intp_coefs)
    endif
  else
    if (xct%subsample_line) then
      call intkernel(crys,kg_fi,kp_fi,syms,xct,hbse_a,dcc,dvv,kco,fi2co_wfn,flag,gvec,intp_coefs,hbse_b=hbse_b,&
        dcc_sub=dcc_sub,dvv_sub=dvv_sub,closepts_sub=closepts_sub)
    else
      call intkernel(crys,kg_fi,kp_fi,syms,xct,hbse_a,dcc,dvv,kco,fi2co_wfn,flag,gvec,intp_coefs,hbse_b=hbse_b)
    endif
  endif
  call logit('Done intkernel')
  call timing%stop(timing%intkernel)
  if(allocated(fi2co_wfn))then;deallocate(fi2co_wfn);endif
  if(allocated(dcc))then;deallocate(dcc);endif
  if(allocated(dvv))then;deallocate(dvv);endif
  if(allocated(kco))then;deallocate(kco);endif
  !FHJ: This is for debugging purposes. Please keep this.
  if (.not.xct%tda.and.peinf%npes==1) then
    write(6,'(1x,a)') 'Dumping BSE Hamiltonian'
    call open_file(unit=666, file='hbse_a.dat', status='replace', form='formatted')
    call open_file(unit=667, file='hbse_b.dat', status='replace', form='formatted')
    write(666,'(2(i0,1x))') xct%nkpt_fi*ncvs_fi, peinf%nblocks*peinf%block_sz
    write(667,'(2(i0,1x))') xct%nkpt_fi*ncvs_fi, peinf%nblocks*peinf%block_sz
    do jj=1,peinf%nblocks*peinf%block_sz
      write(666,'(es16.8)') hbse_a(:, jj)
      write(667,'(es16.8)') hbse_b(:, jj)
    enddo
    call close_file(666)
    call close_file(667)
    write(6,'(1x,a)') 'Done dumping Hamiltonian.'
  endif
!--------------------------------
! Construct mean field header
  version = VER_BSE_FORT
  call init_mf_header_from_types(mf_header, 'EVC', 1, version, kp_fi, gvec, syms, crys)
!--------------------------------
! Exact diagonalization
  if (xct%algo==BSE_ALGO_DIAG .or. xct%algo==BSE_ALGO_DIAG_PRIMME) then
    ! Setup eigenvectors and allocate memory
    call evecs%init_from_xctinfo(xct, kg_fi, nmat, neig, pblock, mf_header=mf_header)
    call evecs%alloc()
    call logit('Calling diagonalize')
    call timing%start(timing%diagonalize)
    if (xct%tda) then
      call diagonalize(xct, neig, nmat, hbse_a, evecs%evals, evecs%u_r)
    else
      call diagonalize(xct, neig, nmat, hbse_a, evecs%evals, evecs%u_r, &
                       hbse_b=hbse_b, evecs_l=evecs%u_l)
      if (peinf%inode==0) then
        write(6,'(/,1x,a,i0)') 'Number of positive eigenvalues: ', count(evecs%evals>0d0)
      endif
    endif
    call timing%stop(timing%diagonalize)
    !--------------------------------
    ! Calculate transition matrix elements
    ! oscillator strength = 2 * cs / omega, as defined below
    call timing%start(timing%trans_mtxel)
    call logit('Computing transition matrix elements')
    allocate(cs (neig,xct%npol))
    cs = 0d0
    allocate(dipoles_r (neig,xct%npol))
    dipoles_r = 0.0d0
    if (.not.xct%tda) then
      allocate(dipoles_l (neig,xct%npol))
      dipoles_l = 0.0d0
      allocate(cs_full (neig,xct%npol))
      cs_full = 0.0d0
    endif
    ! The factor of 1/sqrt(dble(xct%nspin)) below is required to obtain the same
    ! transition matrix elements for the singlet excitons for nspin = 1 and nspin = 2,
    ! See just after eq. (25) in Rohlfing and Louie, PRB 62, 4927 (2000)
    do ipol=1,xct%npol
      if (xct%tda) then
        do ii=1,peinf%neig(peinf%inode+1)
          jj = peinf%peig(peinf%inode+1,ii)
          dipoles_r(jj,ipol) = sum(evecs%u_r(1:nmat,ii)*(dip%s1(1:nmat,ipol))) / sqrt(dble(xct%nspin))
          cs(jj,ipol) = ((dipoles_r(jj,ipol))**2)
        enddo
      else
        do ii=1,peinf%neig(peinf%inode+1)
          jj = peinf%peig(peinf%inode+1,ii)
          ! FHJ: contributions from positive transitions
          dipoles_l(jj,ipol) = sum(evecs%u_l(1:nmat,ii)*(dip%s1(1:nmat,ipol))) / sqrt(dble(xct%nspin))
          dipoles_r(jj,ipol) = sum(evecs%u_r(1:nmat,ii)*(dip%s1(1:nmat,ipol))) / sqrt(dble(xct%nspin))
          ! FHJ: contributions from negative transitions. Some notes:
          ! (1) For the velocity matrix elements: s_(c->v) = - s_(v->c)^* . Proof:
          ! 1st-order expand the wfns |ck> -> |ck+q> and |vk> -> |vk+q> and project
          ! onto <vk| and <ck|. Note that there`s a sign flip in the energy denom.
          ! (2) There`s a negative sign in dipoles_r for (c->v) transitions, which
          ! originates from Fermi factors. See eqn 8 in PRL 80, 4510 (1998).
          dipoles_l(jj,ipol) = dipoles_l(jj,ipol) &
            - sum(evecs%u_l(nmat+1:2*nmat,ii)*(-dip%s1(1:nmat,ipol))) / sqrt(dble(xct%nspin))
          dipoles_r(jj,ipol) = dipoles_r(jj,ipol) &
            + sum(evecs%u_r(nmat+1:2*nmat,ii)*(-dip%s1(1:nmat,ipol))) / sqrt(dble(xct%nspin))
          cs_full(jj,ipol) = (dipoles_l(jj,ipol)) * dipoles_r(jj,ipol)
        enddo
      endif
    enddo
    if (.not.xct%tda) then
      if (peinf%inode==0) then
        cs = dble(cs_full)
      endif
    endif
    call timing%stop(timing%trans_mtxel)
    ! Convert eigenvalues to eV and write them out
    evecs%evals(:) = evecs%evals(:)*ryd
    if (.not.xct%tda) then
      call write_eigenvalues(xct,flag,neig,vol,evecs%evals,cs,dipoles_r,dipoles_l=dipoles_l)
    else
      call write_eigenvalues(xct,flag,neig,vol,evecs%evals,cs,dipoles_r)
    endif
    if(allocated(dipoles_r))then;deallocate(dipoles_r);endif
    if (.not.xct%tda) then
      if(allocated(dipoles_l))then;deallocate(dipoles_l);endif
    endif
    !------------------------------
    ! Calculate the absorption and density of excitonic states
    call timing%start(timing%absorp)
    call logit('Calling absp')
    if (peinf%inode==0) then
      do ipol=1,xct%npol
        call absp(xct, neig, cs(:, ipol), evecs%evals, vol, omega_plasma, flag, ipol)
      enddo
    endif
    call timing%stop(timing%absorp)
  !------------------------------
  ! Write out eigenvectors to file if needed
    call timing%start(timing%write_eig)
    if (flag%eig/=0) then
      call logit('Calling write_eigenvectors')
      call evecs%write_file(nwrite=flag%eig)
    endif
    call timing%stop(timing%write_eig)
    if(allocated(cs))then;deallocate(cs);endif
    if (.not.xct%tda) then
      if(allocated(cs_full))then;deallocate(cs_full);endif
    endif
    call evecs%free()
  elseif (xct%algo==BSE_ALGO_LANCZOS) then
  endif
  ! xct%algo==BSE_ALGO_DIAG
!-------------------------------
! Deallocate stuff
  call dip%free()
  if(allocated(peinf%neig))then;deallocate(peinf%neig);endif
  if(allocated(peinf%peig))then;deallocate(peinf%peig);endif
  call dealloc_grid(kg_fi)
  call kp_fi%free()
  call kpq_fi%free()
  if (eqp%spl_tck%n>0) then
    if(associated(eqp%spl_tck%t))then;deallocate(eqp%spl_tck%t);nullify(eqp%spl_tck%t);endif
    if(associated(eqp%spl_tck%c))then;deallocate(eqp%spl_tck%c);nullify(eqp%spl_tck%c);endif
  endif
  if(allocated(hbse_a))then;deallocate(hbse_a);endif
  if (.not.xct%tda) then
    if(allocated(hbse_b))then;deallocate(hbse_b);endif
  endif
  if(associated(eqp%ecqp))then;deallocate(eqp%ecqp);nullify(eqp%ecqp);endif
  if(associated(eqp%evqp))then;deallocate(eqp%evqp);nullify(eqp%evqp);endif
  if(allocated(intp_coefs))then;deallocate(intp_coefs);endif
  if (xct%iwritecoul .eq. 1 .and. peinf%inode .eq. 0) then
    call close_file(19) ! file vcoul
  endif
  call diag_end()
 
  return
contains
  subroutine diag_end()
   
   
    return
  end subroutine diag_end
end subroutine diag
end module diag_m
