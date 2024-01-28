!===========================================================================
!
! Routines:
!
! (1) epsdiag Originally By MLT Last Modified 5/6/2008 (JRD)
!
! input: crys, gvec, syms, xct types
! read_epsdiag determine which file to be read
!
! output: epsi type, contains the head of the dielectric matrix
! epsi%nq = number of q-vectors stored
! epsi%eps = dielec. matrix at a given q-vector
! epsi%q = coordinates of a q-vector
! epsi%emax = maximum length of the stored q-vectors
!
! Subroutine reads dielectric matrices from unit 10 ('eps0mat')
! and unit 11 ('epsmat') in formatted form and stores the
! head part only in epsi%eps.
!
! The stored q-vectors span a sphere of radius given by epsi%emax
! (that, in general, contains at least the first BZ).
!
! no symmetrization around q=0 is performed
!
! If read_epsdiag, read epsi%eps from unit 12 ("epsdiag.dat") instead
!
! Non-parallelized subroutine
!
!===========================================================================

module epsdiag_m
  use checkbz_m
  use fullbz_m
  use global_m
  use gmap_m
  use input_utils_m
  use misc_m
  use sort_m
  implicit none
  private
  public :: read_epsdiag, write_epsdiag, calc_emax, epsdiag
contains
  !---------------------------------------------------------------------------------------------------
  !> Read data from "epsdiag.dat"
  subroutine read_epsdiag(epsi)
    type(epsinfo), intent(out) :: epsi
    integer :: iq
   
    write(6,'(/1x,a)') 'Reading epsdiag.dat'
    call open_file(12,file='epsdiag.dat',form='formatted',status='old')
    read(12,*) epsi%nq, epsi%emax, epsi%q0vec(:), epsi%epshead
    allocate(epsi%eps (epsi%nq))
    allocate(epsi%q (3,epsi%nq))
    do iq = 1, epsi%nq
      read(12,*) epsi%eps(iq), epsi%q(1:3, iq)
    enddo
    call close_file(12)
    write(6,'(1x,a,f0.6)') '- Dielectric constant: ', 1.d0/epsi%eps(1)
   
    return
  end subroutine read_epsdiag
  !---------------------------------------------------------------------------------------------------
  !> Write data to "epsdiag.dat"
  subroutine write_epsdiag(epsi)
    type(epsinfo), intent(in) :: epsi
    integer :: iq
   
    call open_file(12,file='epsdiag.dat',form='formatted',status='replace')
    write(12,*) epsi%nq, epsi%emax, epsi%q0vec(:), epsi%epshead
    do iq = 1, epsi%nq
      write(12,*) epsi%eps(iq), epsi%q(1:3, iq)
    enddo
    call close_file(12)
   
    return
  end subroutine write_epsdiag
  !---------------------------------------------------------------------------------------------------
  !> emax is some large energy, not smaller than the length squared of any unit lattice vector.
  real(DP) function calc_emax(bdot)
    real(DP), intent(in) :: bdot(:,:)
    integer :: ii
   
    calc_emax = 0.d0
    do ii = 1, 3
      if (calc_emax < bdot(ii,ii)) calc_emax = bdot(ii,ii)
    enddo
    calc_emax = 1.5d0 * calc_emax
   
  end function calc_emax
  !---------------------------------------------------------------------------------------------------
  !> Determine if a G vector is in a periodic or aperiodic direction.
  logical function keep_gvec(gvec, xct, ig)
    type(gspace), intent(in) :: gvec
    type(xctinfo), intent(in) :: xct
    integer, intent(in) :: ig
    ! no push/pop, called too often
    keep_gvec = all(gvec%components(1:3,ig)==0 .or. xct%is_periodic(1:3))
  end function keep_gvec
!---------------------------------------------------------------------------------------------------
! This routine is called only by rank==0
subroutine epsdiag(crys,gvec,syms,epsi,xct)
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (symmetry), intent(in) :: syms
  type (epsinfo), intent(out) :: epsi
  type (xctinfo), intent(in) :: xct
  type (grid) :: qg
  integer :: igamma,nrq0,nrq1,nmtx,ng,ngmax
  integer, allocatable :: isrtq(:),isrtqi(:)
  real(DP) :: q0t(3,1),qk(3)
  real(DP), allocatable :: q1(:,:),eknq(:),eknqa(:,:)
  real(DP), allocatable :: eps(:)
  character :: ajname*6,adate*10
  character :: tmpfn*16
  character :: filenameh5*80
  character :: filenameh50*80
  integer :: ii,jj,nold,itot
  integer :: ifq,ifqt,ig,irq,nqtmax,qgrid(3),idummy,idummy2,idummy3,freq_dep
  real(DP) :: gmax_in,qlength,eaux,qshift(3)
  integer, allocatable :: ind(:), ifqta(:), inda(:,:), isrtqia(:,:)
  integer, allocatable :: oldx(:),oldy(:),oldz(:),old(:,:),isrtinvdummy(:),isrtold(:)
  real(DP), allocatable :: ekold(:),qt(:,:),epsdiagt(:,:)
  real(DP), allocatable :: ph(:),epst(:)
  logical :: skip_checkbz
  logical :: file_exists
  integer, allocatable :: nmtx_of_q(:),nmtx0_of_q(:)
 
  qgrid(:)=0
  epsi%epshead = 0.0d0
  allocate(ind (gvec%ng))
  allocate(ph (gvec%ng))
  allocate(isrtq (gvec%ng))
  allocate(isrtqi (gvec%ng))
  allocate(eknq (gvec%ng))
  ! Read information for inverse dielectric matrix for q->0 unit10
  call open_file(unit=10,file='eps0mat',form='unformatted',status='old')
  call open_file(unit=11,file='epsmat',form='unformatted',status='old',iostat=igamma)
  read(10)
  read(10) ii
  if (ii.ne.0) call die('eps0mat is full-frequency')
  read(10) (qgrid(ii),ii=1,3)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10) nrq0,(epsi%q0vec(ii),ii=1,3)
  read(10) nold
  read(10) ng
  call close_file(10)
  call open_file(unit=10,file='eps0mat',form='unformatted',status='old')
  if(igamma.eq.0) then
    read(11)
    read(11) ii
    if (ii.ne.0) call die('epsmat is full-frequency')
    read(11) (qgrid(ii),ii=1,3)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11) nrq1
    call close_file(11)
    call open_file(unit=11,file='epsmat',form='unformatted',status='old')
    allocate(q1 (3,nrq1))
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11) nrq1,((q1(ii,jj),ii=1,3),jj=1,nrq1)
    call close_file(11)
    call open_file(unit=11,file='epsmat',form='unformatted',status='old')
  else
    nrq1=0
  endif
  ! Store the coordinates of q-vectors in the IBZ
  qg%nr=1+nrq1
  allocate(qg%r (3,qg%nr))
  ! JRD: Remove q0 vec for purpose of generating full zone
  ! we don`t want all the vectors q0 would generate.
  ! qg%r(1:3,1)=q0(1:3)
  qg%r(1:3,1)=0d0
  if(nrq1.ne.0) then
    qg%r(1:3,2:qg%nr)=q1(1:3,1:nrq1)
    if(allocated(q1))then;deallocate(q1);endif
  endif
  !-------------------
  ! Generate full brillouin zone from irreducible wedge
  ! rq -> fq
  call fullbz(crys,syms,qg,syms%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  qshift(:)=0.0d0
  if (igamma.ne.0) then
    tmpfn='eps0mat'
  else
    tmpfn='epsmat'
  endif
  if (.not. skip_checkbz) then
    call checkbz(qg%nf,qg%f,qgrid,qshift,crys%bdot, &
      tmpfn,'q',.true.,xct%freplacebz,xct%fwritebz)
  endif
  ! JRD: Leave q0 vector as exactly 0 or else interpolation will ignore it.
  ! qg%f(1:3,1)=q0(1:3)
  !----------- Read q->0 dielectric matrix -------------------------------
  allocate(oldx (nold))
  allocate(oldy (nold))
  allocate(oldz (nold))
  allocate(isrtold (ng))
  allocate(ekold (ng))
  read(10) ajname,adate
  read(10)
  read(10)
  read(10)
  read(10)
  read(10)
  read(10) gmax_in
  read(10) nrq0,(epsi%q0vec(ii),ii=1,3)
  read(10) nold,(oldx(ii),oldy(ii),oldz(ii),ii=1,nold)
  if(xct%use_hdf5) then
    write(6,'(/1x,a)') 'Epsilon matrix for q->0 read from eps0mat.h5:'
  else
    write(6,'(/1x,a)') 'Epsilon matrix for q->0 read from eps0mat:'
  endif
  write(6,'(1x,a,i0)') '- Number of q-points: ', nrq0
  write(6,'(1x,a,f0.2)') '- Dielectric cutoff (Ry): ', gmax_in
  if (peinf%verb_high) then
    write(6,'(1x,a,3(1x,f10.6))') '- q0 vector:', epsi%q0vec(1:3)
  endif
  read(10) ng,nmtx,(isrtold(ii),jj,ii=1,ng)
  read(10) (ekold(ii),ii=1,ng)
  read(10) (qk(ii),ii=1,3)
  ! Build the sorting/inverse sorting arrays
  isrtq= 0
  isrtqi= 0
  epsi%emax = calc_emax(crys%bdot)
  eaux= 3.d0*epsi%emax
  ! eaux= gmax_in
  call get_ngmax(is_q0=.true.)
  ! Does nqtmax really have to be so big?
  nqtmax = 4*ngmax*qg%nf
  if(allocated(isrtold))then;deallocate(isrtold);endif
  if(allocated(ekold))then;deallocate(ekold);endif
  if(allocated(oldx))then;deallocate(oldx);endif
  if(allocated(oldy))then;deallocate(oldy);endif
  if(allocated(oldz))then;deallocate(oldz);endif
  ! Initializes arrays epst,qt
  allocate(epst (nqtmax))
  allocate(qt (3,nqtmax))
  epst=0.d0
  epsi%nq=0
  qt=0.d0
  ifqt = -1
  do ifq=1,qg%nf
    if (qg%indr(ifq).eq.1) then
      if (ifqt .ne. -1) call die('ifqt multiply defined?')
      ifqt = ifq
    endif
  enddo
  call kinetic_energies(gvec, crys%bdot, eknq, qvec = qg%f(:, ifqt))
  ind=0
  ph=0.0d0
  call gmap(gvec,syms,ngmax,qg%itran(ifqt), &
    qg%kg0(:,ifqt),isrtq,isrtqi,ind,ph,.true.)
  allocate(eps (nmtx))
  do jj=1,nmtx
    read(10) (eps(ii),ii=1,nmtx)
    if (jj .eq. 1) epsi%epshead = eps(1)
    ! If (G+q0) has length greater than emax, store the corresponding
    ! (G+q0),epsinv(G+q0,G+q0) in qt,epst
    !
    ! At the end of the loop, epst will have epsi%nq stored matrix elements
    ! gmap is used to find out the right g-vector ind
    ! JRD: Time HAZARD we have a loop over ig inside loop over jj, probably bad!
    ! DAS: why do we need to do this loop rather than just keep the head, G=G`=0?
    ! JRD: well we need eps_q(0,0), but to interpolate near the edge of the zone we might need eps_q=0(1,1)
    ! because that is equivalent to eps_q=1(0,0) etc... eps_q(G,G`) = eps(q+G,q+G`)
    ! FHJ: but we only have to do this in the periodic directions!
    ! We don`t need the extra G`s from the vacuum direction.
    do ig=1,ngmax
      if (eknq(ig)<epsi%emax .and. isrtqi(ig)/=0) then
        if (keep_gvec(gvec, xct, ig) .and. ind(isrtqi(ig))==jj) then
          epsi%nq = epsi%nq + 1
          if (epsi%nq.gt.nqtmax) then
            write(0,999) epsi%nq,nqtmax,epsi%emax,eaux
            call die('nqtmax too small! No space to store epsdiag.')
          endif
          epst(epsi%nq) = eps(jj)
          qt(:,epsi%nq) = gvec%components(:,ig) + qg%f(:,ifqt)
        endif
      endif
    enddo
  enddo
  if(allocated(eps))then;deallocate(eps);endif
  call close_file(10)
999 format('nqt =',i8,1x,'nqtmax =',i8,1x,'emax =',f8.3,1x,'eaux =',f8.3)
  if (nrq0.gt.1) call die('There is more than one q-point in eps0mat.')
  !---------- Read data for dielectric matrices from unit11 for q<>0 ------------------------
  if(igamma.ne.0) then
    nrq1=0
  else
    if(xct%use_hdf5) then
    ! We do this above.
    ngmax = nold
    else
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11) nrq1
    read(11) nold
    ngmax= 0
    do ii=1,nrq1
      read(11) ng,nmtx
      read(11)
      read(11)
      do jj = 1,nmtx
        read(11)
      end do
      if (ng.gt.ngmax) ngmax= ng
    enddo
    call close_file(11)
    endif
    allocate(oldx (nold))
    allocate(oldy (nold))
    allocate(oldz (nold))
    allocate(isrtold (ngmax))
    allocate(ekold (ngmax))
    call open_file(unit=11,file='epsmat',form='unformatted',status='old')
    read(11) ajname,adate
    read(11)
    read(11)
    read(11)
    read(11)
    read(11)
    read(11) gmax_in
    read(11) nrq1
    read(11) nold,(oldx(ii),oldy(ii),oldz(ii),ii=1,nold)
    if(xct%use_hdf5) then
      write(6,'(/1x,a)') 'Epsilon matrix for q/=0 read from epsmat.h5:'
    else
      write(6,'(/1x,a)') 'Epsilon matrix for q/=0 read from epsmat:'
    endif
    write(6,'(1x,a,i0)') '- Number of q-points: ', nrq1
    write(6,'(1x,a,i0)') '- Number of q-points in the full BZ: ', qg%nf
    write(6,'(1x,a,f0.2)') '- Dielectric cutoff (Ry): ', gmax_in
    write(6,'(1x,a,i0)') '- Number of G-vectors for interpolation of epsinv_00(q+G): ', ngmax
    if (peinf%verb_high) then
      write(6,'(1x,a)') '- Q-points:'
      write(6,'(1(3x,3(1x,f10.6)))') qg%r(1:3,2:nrq1+1)
    endif
    !--------------------------
    ! Read inverse dielectric matrices from unit11 for q<>0
    do irq=2,nrq1+1
      isrtold=0
      ekold=0.d0
      read(11) ng,nmtx,(isrtold(ii),jj,ii=1,ng)
      read(11) (ekold(ii),ii=1,ng)
      read(11) (qk(ii),ii=1,3)
      isrtqi=0
      isrtq=0
      call get_ngmax(is_q0=.false.)
      itot=0
      do ifq=2,qg%nf
        if (irq.eq.qg%indr(ifq)) then
          itot=itot+1
        endif
      enddo
      allocate(ifqta (itot))
      allocate(eknqa (gvec%ng,itot))
      allocate(isrtqia (gvec%ng,itot))
      allocate(inda (gvec%ng,itot))
      eknqa= 0.d0
      itot = 0
      inda=0
      do ifq=2,qg%nf
        if (irq.eq.qg%indr(ifq)) then
          itot=itot+1
          ifqta(itot)=ifq
          call kinetic_energies(gvec, crys%bdot, eknqa(:, itot), qvec = qg%f(:, ifq))
          ph=0.0d0
          ind = 0
          call gmap(gvec,syms,ngmax,qg%itran(ifq), &
            qg%kg0(:,ifq),isrtq,isrtqi,ind,ph,.true.)
          inda(:,itot)=ind
          isrtqia(:,itot)=isrtqi
        endif
      enddo
      allocate(eps (nmtx))
      do jj = 1,nmtx
        read (11) (eps(ii),ii=1,nmtx)
        ! If (G+q0) has length greater than emax, store the corresponding
        ! (G+q0),epsinv(G+q0,G+q0) in qt,epst
        !
        ! At the end of the loop, epst will have epsi%nq stored matrix elements
        ! gmap is used to find out the right g-vector
        ! JRD: Possible Time HAZARD
        do ii = 1 , itot
          do ig=1,ngmax
            if (eknqa(ig,ii)<epsi%emax .and. isrtqia(ig,ii)/=0) then
              if (keep_gvec(gvec, xct, ig) .and. inda(isrtqia(ig,ii),ii)==jj) then
                epsi%nq = epsi%nq + 1
                if (epsi%nq.gt.nqtmax) then
                  write(0,*) epsi%nq,nqtmax,epsi%emax,eaux
                  call die('nqtmax too small! Not enough space to store epsdiag.')
                endif
                epst(epsi%nq)= eps(jj)
                qt(:,epsi%nq) = gvec%components(:,ig) + qg%f(:,ifqta(ii))
              endif
            endif
          enddo
        enddo
      enddo
      if(allocated(eps))then;deallocate(eps);endif
      if(allocated(ifqta))then;deallocate(ifqta);endif
      if(allocated(eknqa))then;deallocate(eknqa);endif
      if(allocated(isrtqia))then;deallocate(isrtqia);endif
      if(allocated(inda))then;deallocate(inda);endif
    enddo
    if(.not. xct%use_hdf5) call close_file(11)
    if(allocated(isrtold))then;deallocate(isrtold);endif
    if(allocated(ekold))then;deallocate(ekold);endif
    if(allocated(oldx))then;deallocate(oldx);endif
    if(allocated(oldy))then;deallocate(oldy);endif
    if(allocated(oldz))then;deallocate(oldz);endif
  endif
  call dealloc_grid(qg)
  if(allocated(isrtq))then;deallocate(isrtq);endif
  if(allocated(isrtqi))then;deallocate(isrtqi);endif
  if(allocated(eknq))then;deallocate(eknq);endif
  if(allocated(ind))then;deallocate(ind);endif
  if(allocated(ph))then;deallocate(ph);endif
  if(xct%use_hdf5) then
  if(allocated(nmtx0_of_q))then;deallocate(nmtx0_of_q);endif
  if (igamma .eq. 0) then
    if(allocated(nmtx_of_q))then;deallocate(nmtx_of_q);endif
  endif
  endif
  !------------------------
  ! Transfer data from epst to epsi%eps. epsi%eps stores the real
  ! part of epst (its imaginary part should be zero!).
  ! Sorting according to length of epsi%q
  allocate(epsi%eps (epsi%nq))
  allocate(epsi%q (3,epsi%nq))
  allocate(eknq (epsi%nq))
  allocate(isrtq (epsi%nq))
  isrtq=0
  do ii = 1, epsi%nq
    eknq(ii) = DOT_PRODUCT(qt(1:3, ii), MATMUL(crys%bdot(1:3, 1:3), qt(1:3, ii)))
  enddo
  call sortrx(epsi%nq, eknq, isrtq)
  ! Store data in epsi and write it out
  do ii=1,epsi%nq
    epsi%eps(ii) = dble( epst(isrtq(ii)) )
    epsi%q(:,ii) = qt(:,isrtq(ii))
    qlength = sqrt(DOT_PRODUCT(qt(1:3, isrtq(ii)), MATMUL(crys%bdot(1:3, 1:3), qt(1:3, isrtq(ii)))))
  enddo
  write(6,'(1x,a,f0.6)') '- Dielectric constant: ', 1.d0/epsi%eps(1)
  write(6,'(1x,a,f0.2)') '- Max. |q+G|^2 to interp. epsinv (Ry): ', epsi%emax
  write(6,'(1x,a,i0)') '- Number of q+G points to interp. epsinv: ', epsi%nq
  call write_epsdiag(epsi)
  if(allocated(qt))then;deallocate(qt);endif
  if(allocated(epst))then;deallocate(epst);endif
  if(allocated(eknq))then;deallocate(eknq);endif
  if(allocated(isrtq))then;deallocate(isrtq);endif
 
  return
contains
  !> ngmax = number of G-vectors within the sphere of radius emax
  !! (if nqtmax is too small, increase it)
  !! FHJ: Calculate the number of G vectors so that |q+G|^2 < emax.
  !! See DAS/JRD` discussion why we do this.
  subroutine get_ngmax(is_q0)
    logical, intent(in) :: is_q0
    integer :: gg(3),iout
   
    ngmax = 0
    do ii=1,ng
      if (ekold(isrtold(ii)).lt.eaux) then
        gg(1)=oldx(isrtold(ii))
        gg(2)=oldy(isrtold(ii))
        gg(3)=oldz(isrtold(ii))
        call findvector(iout,gg,gvec)
        if (iout.gt.gvec%ng) then
          write(0,'(a,2i6)') 'WARNING: eps ',iout,gvec%ng
        endif
        if (iout.le.0) then
          write(0,'(a,4i6)') 'WARNING: eps ',iout,gg(1:3)
        endif
        isrtq(ii)=iout
        isrtqi(iout)=ii
        if (.not.is_q0) eknq(ii)=ekold(isrtold(ii))
        ngmax=ii
      else
        exit
      endif
    enddo
   
  end subroutine get_ngmax
end subroutine epsdiag
end module epsdiag_m
