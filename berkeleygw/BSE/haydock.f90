!=======================================================================
!
! Routines:
!
! 1. haydock Originally By MLT Last Modified 6/24/08 (JRD)
!
! Calculates the real and imaginary parts of the macroscopic dielectric
! function starting from interaction matrix elements calculated by
! the Kernel code. It uses interpolation in the matrix elements and
! the Haydock recursion method. Spin-polarized case implemented.
!
! For more details, see:
! Rohlfing & Louie, PRB 62:(8), 4927 (2000)
! Benedict & Shirley, PRB 59:(8), 5441 (1999)
! R. Haydock, Comput. Phys. Commun. 20, 11 (1980)
! G. Strinati, Rivista del Nuovo Cimento 11:(12), 1 (1988)
!
! Please report bugs to: jdeslip@civet.berkeley.edu
!
!========================================================================

module haydock_m
  use absh_m
  use absp0_m
  use blas_m
  use fullbz_m
  use genwf_m
  use global_m
  use input_fi_m
  use input_q_m
  use intkernel_m
  use intwfn_m
  use iterate_m
  use misc_m
  use mtxel_jdos_m
  use random_m
  use timing_m, only: timing => bse_timing
  use vmtxel_m
  implicit none
  private
  public :: &
    haydock
contains
subroutine haydock(eqp,xct,flag,nmax)
  type (eqpinfo), intent(inout) :: eqp
  type (xctinfo), intent(inout) :: xct
  type (flags), intent(inout) :: flag
  integer, intent(in) :: nmax
  type (crystal) :: crys
  type (mmtsinfo) :: mmts
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (grid) :: kg_fi, kgq_fi,kg_co,kgq_co
  type (kpoints) :: kp_fi, kpq_fi
  type (wavefunction) :: wfnc_fi
  type (wavefunction) :: wfnvq_fi
  type (work_genwf) :: work, workq
  type (int_wavefunction) :: intwfnc
  type (int_wavefunction) :: intwfnv
  type (vmtxel_t) :: dip
  character :: tmpstr*128,filename*20
  integer :: ii,jj,ncount,nmat,ncvs_fi
  integer :: ikb, icb, ivb
  integer :: ik,ikq,ikt,iblock,ikrq,ikcvs,ikcvsd,ic,iv,is
  integer :: seed, iunit_c,iunit_v
  integer :: ih,ith,n0
  real(DP) :: vol,vol1,sum_,sum1,dtemp
  real(DP) :: tsec(2),tmin(2),tmax(2),omega_plasma
  character*16, allocatable :: routnam(:)
  integer, allocatable :: routsrt(:)
  integer, allocatable :: fi2co_wfn(:,:),indexq_fi(:)
  real(DP), allocatable :: kco(:,:)
  real(DP), allocatable :: &
    dcc(:,:,:,:,:),dvv(:,:,:,:,:),s1(:),s1_l(:),s1k(:,:,:),dummy(:), &
    s0(:),s0_l(:),hqpcc(:,:,:),hqpvv(:,:,:)
  !> (kcvs, k`c`v`s`), "A" block of BSE Hamiltonian
  real(DP), allocatable :: hbse_a(:,:)
  real(DP), allocatable :: intp_coefs(:,:)
 
  mmts%nmax = nmax
  if (xct%read_kpoints .and. peinf%inode.eq.0) then
    write(0,*) 'WARNING: you are using partial sampling with the Haydock iterative scheme!'
    write(0,*) 'The optical spectrum may not be meaningful...'
  endif
  if(mmts%nmax.eq.0.and.flag%spec.ne.1) then
    call die('missing number of iterations')
  endif
! If flag%spec.eq.1, just calculate the spectrum from input coefficients
  ith=21
  if (flag%spec.eq.1) then
    if (peinf%inode.eq.0) then
      call open_file(ith,file='eps2_moments',form='unformatted',status='old')
      read(ith) mmts%nmax,mmts%norm,mmts%vol,ii,xct%nspin
      allocate(mmts%an (mmts%nmax))
      allocate(mmts%bn (mmts%nmax))
      read(ith) (mmts%an(ii),ii=1,mmts%nmax)
      read(ith) (mmts%bn(ii),ii=1,mmts%nmax)
      write(6,*)
      write(6,*) 'Reading coefficients from file; nmax = ',mmts%nmax, ' norm= ',mmts%norm
      write(6,*) '     n    a(n)        b(n)  '
      do ii=1,mmts%nmax
        write(6,120) ii, mmts%an(ii),mmts%bn(ii)
      enddo
      call close_file(ith)
    endif
120 format(3x,i4,3x,2e12.5)
  else
!-----------------------
! Read wavefunctions on the fine grid
  call logit('Calling input')
  call timing%start(timing%input)
  call input_fi(crys,gvec,kg_fi,kp_fi,syms,eqp,xct,flag, &
    omega_plasma,.false.,intwfnc)
  nmat = xct%nkpt_fi*xct%ncb_fi*xct%nvb_fi*xct%nspin
  vol = xct%nktotal*crys%celvol
  if (peinf%inode.eq.0) then
    write(6,'(/1x,a)') 'More job parameters:'
    write(6,'(1x,a,es9.3e2)') '- Crystal volume (bohr): ', vol
    write(6,'(1x,a,f0.3)') '- Broadening (eV): ', xct%eta
    write(6,'(1x,a,i0)') '- Number of valence bands: ', xct%nvb_fi
    write(6,'(1x,a,i0)') '- Number of cond. bands: ', xct%ncb_fi
    write(6,'(1x,a,i0)') '- Number of spins: ', xct%nspin
    write(6,'(1x,a,i0)') '- Number of Haydock iterations: ', mmts%nmax
    write(6,'()')
  endif
  call timing%stop(timing%input)
  allocate(indexq_fi (xct%nkpt_fi))
  allocate(xct%indexq_fi (xct%nkpt_fi))
  if (flag%vm.eq.0 .or. .not. flag%read_dtmat) then
    call timing%start(timing%input_q)
    call logit('Calling input_q')
    call input_q(kp_fi,crys,gvec,kg_fi,kgq_fi,kpq_fi, &
      syms,xct,indexq_fi,eqp,flag,intwfnv)
    call timing%stop(timing%input_q)
  endif
!----------------------------
! Calculate the transformation matrices from coarse-grid wavefunctions
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
  call intwfn(kp_fi,crys,syms,xct,flag,gvec,kg_fi,kgq_fi,kg_co, kgq_co, &
    dcc,dvv,kco,fi2co_wfn,indexq_fi,eqp,intwfnv,intwfnc,intp_coefs)
  call timing%stop(timing%intwfn)
  if(associated(xct%ifmax))then;deallocate(xct%ifmax);nullify(xct%ifmax);endif
  if (flag%vm.eq.0 .or. .not. flag%read_dtmat) then
    ! otherwise, we did not call input_q to allocate it
    if(associated(xct%ifmaxq))then;deallocate(xct%ifmaxq);nullify(xct%ifmaxq);endif
  endif
!----------------------------
! Initialize recursion states and coefficients.
! If there is a previous run, read previous states from unit ith and
! skip the calculation of velocity/momentum matrix elements.
!
! peinf%block_sz = size of a distributed column in hbse_a
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
! Since this code always compute an odd number of an's and bn's,
! it will calculate nmax+1 coefficients if nmax is even
! solution: overshoot the array size
  mmts%nmaxp=mmts%nmax+1
  allocate(mmts%an (mmts%nmaxp))
  allocate(mmts%bn (mmts%nmaxp))
  mmts%an(:) = 0.d0
  mmts%bn(:) = 0.d0
  mmts%vol = vol
  allocate(s0 (nmat))
  s0 = 0.d0
  if (flag%vm.eq.2) then
    if (peinf%inode.eq.0) then
      call open_file(ith,file='eps2_moments',form='unformatted',status='old')
      read(ith) n0,mmts%norm,vol1,ii,is
      if ((abs(vol1-mmts%vol).gt.1.d-10).or.ii.ne.nmat.or.is.ne.xct%nspin) then
        call die(' Parameter mismatch in old file eps2_moments')
      endif
      read(ith) (mmts%an(ii),ii=1,n0)
      read(ith) (mmts%bn(ii),ii=1,n0)
      read(ith) (dip%s1(ii,1),ii=1,nmat)
      read(ith) (s0(ii),ii=1,nmat)
      call close_file(ith)
      write(6,*) 'Reading old data from file eps2_moments'
      sum_ = (4.d0*PI_D*ryd)**2*mmts%norm*mmts%an(1)/(mmts%vol*xct%nspin)
      write(6,*) 'Sum rule (excitons) : ',sum_,' eV^2'
      write(6,*) '      n         a(n)        b(n)'
      do ih=1,n0
        write(6,120) ih,mmts%an(ih),mmts%bn(ih)
      enddo
    endif
  else
    n0=1
  endif
!-------------------------
! Calculate the velocity (or momentum) matrix elements.
! Each PE calculates a small number of them. At the end, share data.
!
! If flag%vm.eq.1, skip this part and just read the matrix elements
! from "vmtxel".
  call logit('Calling v/p matrix elements')
  if (flag%vm.eq.0) then
    if (flag%opr .eq. 2) then
      seed =5000
      call genrand_init(put=seed)
      ! Just set the s1 vector to be random numbers on all processors
      ! The seed has to be the same on all processors otherwise we need
      ! to set this on one and broadcast..
      call mtxel_jdos(dip%s1,nmat)
    else
      do ikt=1, peinf%ikt(peinf%inode+1)
        ik = peinf%ik(peinf%inode+1,ikt)
        ikq = indexq_fi(ik)
        ikrq = kg_fi%indr(ik)
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
      if(work%ikold /= 0) then
        if(associated(work%cg))then;deallocate(work%cg);nullify(work%cg);endif
        if(associated(work%isort))then;deallocate(work%isort);nullify(work%isort);endif
        if(associated(work%ph))then;deallocate(work%ph);nullify(work%ph);endif
        if(associated(work%ind))then;deallocate(work%ind);nullify(work%ind);endif
      endif
      if(workq%ikold /= 0) then
        if(associated(workq%cg))then;deallocate(workq%cg);nullify(workq%cg);endif
        if(associated(workq%isort))then;deallocate(workq%isort);nullify(workq%isort);endif
        if(associated(workq%ph))then;deallocate(workq%ph);nullify(workq%ph);endif
        if(associated(workq%ind))then;deallocate(workq%ind);nullify(workq%ind);endif
      endif
      ! Share matrix elements
      call dip%reduce()
    endif ! Whether jdos or (velocity/momentum)
    ! Write file
    call dip%write_vmtxel()
  else if(flag%vm == 1) then
    ! Read dipole matrix elements from file if they were already computed
    call dip%read_vmtxel()
  endif
  if (flag%vm.eq.0 .or. .not. flag%read_dtmat) then
    call dealloc_grid(kgq_fi)
  endif
  if(allocated(indexq_fi))then;deallocate(indexq_fi);endif
! JRD: Now close the no longer needed wavefunction files
  if(associated(intwfnc%cgk))then;deallocate(intwfnc%cgk);nullify(intwfnc%cgk);endif
  if(associated(intwfnv%cgk))then;deallocate(intwfnv%cgk);nullify(intwfnv%cgk);endif
  if(associated(intwfnc%isort))then;deallocate(intwfnc%isort);nullify(intwfnc%isort);endif
  if(associated(intwfnv%isort))then;deallocate(intwfnv%isort);nullify(intwfnv%isort);endif
!-------------------------
! Calculate the non-interacting spectrum. Only one PE works
  call logit('Calling absp0')
  if (peinf%inode.eq.0) call absp0(eqp,xct,dip%s1,vol,omega_plasma,flag,1)
!---------------------------
! Build Hamiltonian matrix. Diagonal part comes from the "non-interacting"
! quasiparticle Hamiltonians. If the QP Greens func is diagonal,
! then these are just the quasiparticle energy differences on the
! diagonal. The more general case is:
!
! <cvk|H0|c'v'k'> = delta(k,k') *
! [ <c|hqp|c'>*delta(v',v) - delta(c,c')*<v'|hqp|v> ]
!
! The rest of the Hamiltonian, which is the electron-hole interaction,
! comes from interpolation further below.
  call logit('Building non-interacting Hamiltonian')
  allocate(hbse_a (xct%nkpt_fi*ncvs_fi,peinf%nblocks*peinf%block_sz))
  hbse_a(:,:) = 0.d0
! Loop over kpoints
  do iblock=1,peinf%ibt(peinf%inode+1)
    ik = peinf%ikb(iblock)
! Build <c|hqp|c'> and <v|hqp|v'> for this kpoint
    allocate(hqpcc (xct%ncb_fi,xct%ncb_fi,xct%nspin))
    allocate(hqpvv (xct%nvb_fi,xct%nvb_fi,xct%nspin))
    hqpcc = 0.0d0
    hqpvv = 0.0d0
! Put QP energies on diagonals of hqpcc and hqpvv to start
    do is=1,xct%nspin
      do ic=1,xct%ncb_fi
        hqpcc(ic,ic,is) = eqp%ecqp(ic,ik,is)
      enddo
      do iv=1,xct%nvb_fi
        hqpvv(iv,iv,is) = eqp%evqp(iv,ik,is)
      enddo
    enddo
! Read possible offdiagonal QP elements from "hqp.<ik>" file
! if it exists JRD: This is broken for now. Someone should fix
! it in the future if they want to use it
        !if (ik.lt.10) then
        ! write(tmpstr,'(a,i1)') 'hqp.',ik
        !else if (ik.lt.100) then
        ! write(tmpstr,'(a,i2)') 'hqp.',ik
        !else if (ik.lt.1000) then
        ! write(tmpstr,'(a,i3)') 'hqp.',ik
        !else if (ik.lt.10000) then
        ! write(tmpstr,'(a,i4)') 'hqp.',ik
        !else
        ! write(0,*) 'too many kpoints for reading hqp'
        !endif
        !call open_file(9,file=tmpstr,form='formatted',status='old',iostat=is)
        !if (is.eq.0) then
        ! if (peinf%inode.eq.0) then
        ! write(6,*) 'Reading offdiagonal hqp from file ',tmpstr
        ! write(6,*) 'All values in eV'
        ! endif
        ! do
        ! read(9,*,end=999) nocc,ii,jj,x,y
!
!! if ii and jj both refer to valence, states, put
!! matrix element into hqpvv
!
! if ((ii<=nocc).and.(ii>nocc-xct%nvb_fi).and. &
! (jj<=nocc).and.(jj>nocc-xct%nvb_fi)) then
! if (peinf%inode.eq.0) write(6,'(a,2i5,2f20.10)') ' hqp(v,vp) = ',ii,jj,x,y
! ii=nocc-ii+1
! jj=nocc-jj+1
! is = 1
!#ifdef CPLX
! hqpvv(ii,jj,is) = cmplx(x,y,kind=DPC)/ryd
!#else
! hqpvv(ii,jj,is) = x/ryd
!#endif
! else if ((ii>nocc).and.(ii<=nocc+xct%ncb_fi).and. &
! (jj>nocc).and.(jj<=nocc+xct%ncb_fi)) then
! if (peinf%inode.eq.0) write(6,'(a,2i5,2f20.10)') ' hqp(c,cp) = ',ii,jj,x,y
! ii=ii-nocc
! jj=jj-nocc
! is = 1
!#ifdef CPLX
! hqpcc(ii,jj,is) = cmplx(x,y,kind=DPC)/ryd
!#else
! hqpcc(ii,jj,is) = x/ryd
!#endif
! endif
! enddo
! 999 call close_file(9)
! write(6,*)
! endif ! if hqp.<ik> was found
! Now build hamiltonian from hqpcc and hqpvv
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
          hbse_a(ikcvs,ikcvsd) = hqpcc(icb,icb,is) - hqpvv(ivb,ivb,is)
        enddo
      enddo
    enddo
    if(allocated(hqpcc))then;deallocate(hqpcc);endif
    if(allocated(hqpvv))then;deallocate(hqpvv);endif
  enddo ! loop on k-points on this processor
!-----------------------------------------------------------------------
! Define the mapping of eigenvectors: the ii-th column of the matrix
! evecs_r(:,:) stored in PE #ipe will contain the eigenvector of order
! peinf%peig(ipe,ii). The total number of eigenvectors stored in
! each processor is given by peinf%neig(1:peinf%npes).
! pblock >= maxval(peinf%neig(1:peinf%npes))
! NOTE: this is taken from diag code, but it also holds for the
! mapping of hbse_a, that has the same layout as evecs_r (Murilo)
  allocate(peinf%neig (peinf%npes))
  allocate(peinf%peig (peinf%npes,peinf%nblocks*peinf%block_sz))
  peinf%neig=0
  peinf%peig=0
  ii=1
  do jj=1,nmat
    if (ii.eq.peinf%npes+1) ii=1
    peinf%neig(ii)=peinf%neig(ii)+1
    peinf%peig(ii,peinf%neig(ii))=jj
    if (mod(jj,peinf%block_sz).eq.0) ii=ii+1
  enddo
!---------------------------------
! Initialize local states : s1_l --> local part of s1
! s0_l --> local part of s0
! Whenever s1 and s0 change, their local parts must be updated
! Initialize s1 from the dipole matrix elements.
  allocate(s1 (nmat))
  s1(:) = dip%s1(:,1)
  call dip%free()
  allocate(s1_l (peinf%nblocks*peinf%block_sz))
  !MJ : CHANGE THIS
  call local_s(nmat,s1,s1_l)
  allocate(s0_l (peinf%nblocks*peinf%block_sz))
  call local_s(nmat,s0,s0_l)
!--------------------------------
! Interpolation scheme in the Kernel
  call logit('Calling intkernel')
  call timing%start(timing%intkernel)
  call intkernel(crys,kg_fi,kp_fi,syms,xct,hbse_a,dcc,dvv,kco,fi2co_wfn,flag,gvec,intp_coefs)
  call timing%stop(timing%intkernel)
  call dealloc_grid(kg_fi)
  if(allocated(fi2co_wfn))then;deallocate(fi2co_wfn);endif
  if(allocated(dcc))then;deallocate(dcc);endif
  if(allocated(dvv))then;deallocate(dvv);endif
  if(allocated(kco))then;deallocate(kco);endif
!------------------------------
! Normalize s1 to get the first-order Haydock state: |1>= (e.v) |Ground>
! Calculate the first a coefficient
  if (flag%vm .ne. 2) then
    ! FHJ: We have to use BLAS because Intel incorrectly optimizes this
    dtemp = blas_nrm2(nmat, s1, 1)
    mmts%norm = dtemp**2
    s1(:) = s1(:) / dtemp
    call local_s(nmat,s1,s1_l)
    sum1 = dble( DOT_PRODUCT(s1,MATMUL(hbse_a,s1_l)) )
    mmts%an(1) = sum1
!--------------------------
! Calculate the second order state, s0, and the first b coefficient
    s0 = MATMUL(hbse_a,s1_l)
    s0(:) = s0(:) - s1 * mmts%an(1)
    ! FHJ: We have to use BLAS because Intel incorrectly optimizes this
    dtemp = blas_nrm2(nmat, s0, 1)
    mmts%bn(1) = dtemp**2
    s0(:) = s0(:) / dtemp
    call local_s(nmat,s0,s0_l)
  endif
  if (peinf%inode.eq.0) then
    write(6,750) mmts%nmax
750 format(/,'Performing Haydock recursion with ',i5,' steps. ')
    write(6,*) 'Norm of first state : ',mmts%norm
! Check plasmon sum rule
!
! Exact value of sum rule:
! sum_ = (pi/2.d0)*( plasma frequency )^2
    sum_ = (4.d0*PI_D*ryd)**2*mmts%norm*mmts%an(1)/ &
      (mmts%vol*xct%nspin)
    write(6,*) 'Sum rule (excitons) : ',sum_,' eV^2'
    write(6,*) '      n         a(n)        b(n)'
    write(6,120) n0,mmts%an(n0),mmts%bn(n0)
  endif
!-----------------------------------------------------------------------
! Start Haydock recursion method. At this point, n0 is always an odd number
! After each pair of iterations, we have the states:
! | ih+1 > --- s1
! | ih+2 > --- s0
! Lower states are lost.
  call timing%start(timing%iterate)
  do ih=n0+1,mmts%nmax,2
    call iterate(mmts,xct,nmat,ih,hbse_a,s1,s0)
  enddo
  call timing%stop(timing%iterate)
  if(allocated(hbse_a))then;deallocate(hbse_a);endif
  if(allocated(s1))then;deallocate(s1);endif
  if(allocated(s0))then;deallocate(s0);endif
  if(allocated(peinf%neig))then;deallocate(peinf%neig);endif
  if(allocated(peinf%peig))then;deallocate(peinf%peig);endif
!---------------------------
! Calculate absorption spectrum using Haydock recursion
  endif
  if (peinf%inode.eq.0) call absh(mmts,xct%nspin,xct%nspinor,xct%eta)
  if(associated(mmts%an))then;deallocate(mmts%an);nullify(mmts%an);endif
  if(associated(mmts%bn))then;deallocate(mmts%bn);nullify(mmts%bn);endif
  if(associated(eqp%eclda))then;deallocate(eqp%eclda);nullify(eqp%eclda);endif
  if(associated(eqp%evlda))then;deallocate(eqp%evlda);nullify(eqp%evlda);endif
  if(allocated(intp_coefs))then;deallocate(intp_coefs);endif
  if (xct%iwritecoul .eq. 1 .and. peinf%inode .eq. 0) then
    call close_file(19) ! file vcoul
  endif
  if (eqp%spl_tck%n>0) then
    if(associated(eqp%spl_tck%t))then;deallocate(eqp%spl_tck%t);nullify(eqp%spl_tck%t);endif
    if(associated(eqp%spl_tck%c))then;deallocate(eqp%spl_tck%c);nullify(eqp%spl_tck%c);endif
  endif
  call kp_fi%free()
  call kpq_fi%free()
 
  return
end subroutine haydock
end module haydock_m
