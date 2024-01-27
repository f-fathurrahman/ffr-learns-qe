!===============================================================================
!
! Utilities:
!
! (1) epsomega Originally By gsm Last Modified 11/18/2010 (gsm)
!
! Parallel code that plots frequency dependence of Plasmon-Pole, Retarded,
! Advanced epsilon for a given q, G, G` vectors. Input parameters are read
! from file epsomega.inp in the working directory.
!
! Note that epsPP constructed with full-frequency epsmat will be
! different from epsPP constructed with static epsmat because of
! finite broadening used in full-frequency epsmat.
!
! FIXME: use WFN_RHO_VXC for reading RHO
!
! FIXME: [2013-01-15 gsm] compact form of GPP model doesn`t work when
! wtilde2 < 0 which is often the case for G .ne. G`; replace with the
! explicit expression from Hybertsen & Louie (see commented out code
! in the GPP section of epsinvomega.f90)
!
!===============================================================================
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
program epsomega
  use global_m
  use input_utils_m
  use inversion_m
  use read_matrix_m
  use scalapack_m
  use epsread_hdf5_m
  implicit none
  logical :: sflag
  integer :: i,j,k,iq,mq,ig,igp,jg,jgp,iunit,ierr
  integer :: icurr,icol,irow,icolm,irowm,icols,irows
  integer :: iout,nFFTgridpts,FFTgrid(3)
  real(DP) :: omega
  character*256, parameter :: fninp = "epsomega.inp"
  character*256 :: fneps,fnrho,fngpp,fnffr,fnffa
  real(DP) :: q(3)
  integer :: g(3),gp(3)
  character :: infile*80
  type(gspace) :: gvec
  integer :: freq_dep,nFreq,nfreq_imag,ii,nq,ng,nmtx,kgrid(3)
  real(DP) :: dDeltaFreq,dBrdning,ecuts,delta,qvec(3)
  real(DP), allocatable :: qpt(:,:)
  real(DP), allocatable :: ekin(:)
  real(DP), allocatable :: dFreqGrid(:)
  integer, allocatable :: isrtx(:)
  integer, allocatable :: isorti(:)
  integer, allocatable :: nmtx_of_q(:)
  real(DP), allocatable :: eps(:,:)
  complex(DPC), allocatable :: dFreqBrd(:)
  complex(DPC), allocatable :: epsPP(:,:,:)
  complex(DPC), allocatable :: epsR(:,:,:)
  complex(DPC), allocatable :: epsA(:,:,:)
  complex(DPC), allocatable :: epsAux(:,:)
  real(DP), allocatable :: rhog(:)
  real(DP) :: bvec(3,3),bdot(3,3),blat,celvol,reccelvol
  integer :: nvecs,nspin
  integer :: nproc_para,num_gvec,gx,gy,gz,qgrid(3),error
  real(DP) :: rho0
  integer, allocatable :: gvec_rho(:,:)
  real(DP), allocatable :: xcdum(:,:)
  real(DP) :: wp2,qg(3),qgp(3),qgqg,qgqgp,lambda,phi
  real(DP) :: Omega2,wtilde2,epsggp,I_epsggp
  complex(DPC) :: wtilde2_temp
  type (scalapack) :: scal
  call peinfo_init()
!-----------------------------
! read input file
  if (peinf%inode.eq.0) then
    write(6,'(/,1x,"reading",1x,a,1x,"file",/)')trim(fninp)
    call open_file(55,file=trim(fninp),form='formatted',status='old')
    read(55,'(a)') fneps
    read(55,'(a)') fnrho
    read(55,*) (q(i),i=1,3)
    read(55,*) (g(i),i=1,3)
    read(55,*) (gp(i),i=1,3)
    read(55,*) nFreq
    read(55,*) dDeltaFreq
    read(55,*) dBrdning
    read(55,'(a)') fngpp
    read(55,'(a)') fnffr
    read(55,'(a)') fnffa
    call close_file(55)
    write(6,'(2a)') "     eps  file  = ", trim(fneps)
    write(6,'(2a)') "     rho  file  = ", trim(fnrho)
    write(6,'(a, 3f15.12)') "     q  vector  = ", q(1:3)
    write(6,'(a, 3i4)') "     G  vector  = ", g(1:3)
    write(6,'(a, 3i4)') "     G' vector  = ", gp(1:3)
    write(6,'(a, i6)') "     nFreq      = ", nFreq
    write(6,'(a, f7.3)') "     dDeltaFreq = ", dDeltaFreq
    write(6,'(a, f7.3)') "     dBrdning   = ", dBrdning
    write(6,'(2a)') "     GPP file   = ", trim(fngpp)
    write(6,'(2a)') "     FFR file   = ", trim(fnffr)
    write(6,'(2a,/)') "     FFA file   = ", trim(fnffa)
  endif
!-----------------------------
! read eps file
  if (peinf%inode.eq.0) then
    write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fneps)
    iunit=12
    call open_file(unit=iunit,file=trim(fneps),form='unformatted',status='old')
    read(iunit)
    read(iunit) freq_dep,ii
    if (freq_dep.ne.0) then
      nFreq=ii
    endif
    read(iunit) (kgrid(i),i=1,3)
    allocate(dFreqGrid (nFreq))
    allocate(dFreqBrd (nFreq))
    if (freq_dep.ne.0) then
      read(iunit) (dFreqGrid(i),i=1,nFreq),(dFreqBrd(i),i=1,nFreq)
      if (nFreq.gt.1) dDeltaFreq=dFreqGrid(2)-dFreqGrid(1)
      dBrdning=aimag(dFreqBrd(1))
    else
      do i=1,nFreq
        dFreqGrid(i)=dDeltaFreq*dble(i-1)
        dFreqBrd(i)=dBrdning
      enddo
      read(iunit)
    endif
    read(iunit)
    read(iunit)
    read(iunit) ecuts
    read(iunit) nq
    read(iunit) ng
    rewind(iunit)
    allocate(qpt (3,nq))
    allocate(gvec%components (3,ng))
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit) nq,((qpt(j,i),j=1,3),i=1,nq)
    read(iunit) ng,((gvec%components(j,i),j=1,3),i=1,ng)
    ig=0
    igp=0
    do i=1,ng
      if (gvec%components(1,i).eq.g(1).and.gvec%components(2,i).eq.g(2).and. &
        gvec%components(3,i).eq.g(3)) ig=i
      if (gvec%components(1,i).eq.gp(1).and.gvec%components(2,i).eq.gp(2).and. &
        gvec%components(3,i).eq.gp(3)) igp=i
    enddo
    mq=-1
    do iq=1,nq
      if (abs(q(1)-qpt(1,iq)).lt.TOL_Zero.and. &
        abs(q(2)-qpt(2,iq)).lt.TOL_Zero.and. &
        abs(q(3)-qpt(3,iq)).lt.TOL_Zero) then
        mq=iq-1
        exit
      endif
    enddo
    if(allocated(qpt))then;deallocate(qpt);endif
    write(6,'(a,i6)') "     omega num  = ", nFreq
  endif
  if (mq.eq.-1) then
    call die("cannot find q vector in file " // trim(fneps))
  endif
  if (ig.eq.0) then
    call die("cannot find G vector in file " // trim(fneps))
  endif
  if (igp.eq.0) then
    call die("cannot find G' vector in file " // trim(fneps))
  endif
  if (peinf%inode.eq.0) then
    do iq=1,mq
      read(iunit) ng,nmtx
      read(iunit)
      read(iunit)
      if (freq_dep.eq.0) then
        do j=1,nmtx
          read(iunit)
        enddo
      else
        do j=1,nmtx
          do i=1,nmtx
            read(iunit)
          enddo
        enddo
      endif
    enddo
    allocate(isrtx (ng))
    allocate(isorti (ng))
    isrtx(:)=0
    isorti(:)=0
    read(iunit) ng,nmtx,(isrtx(i),isorti(i),i=1,ng)
    read(iunit)
    read(iunit)
    ig=isorti(ig)
    igp=isorti(igp)
    if(allocated(isorti))then;deallocate(isorti);endif
  endif
  if (ig.eq.0) then
    call die("cannot find G vector in file " // trim(fneps))
  endif
  if (igp.eq.0) then
    call die("cannot find G' vector in file " // trim(fneps))
  endif
  call blacs_setup(scal, nmtx, .true.)
!-----------------------------
! read distributed matrices
  allocate(eps (scal%npr,scal%npc))
  allocate(epsPP (nFreq,scal%npr,scal%npc))
  if (freq_dep.ne.0) then
    allocate(epsR (nFreq,scal%npr,scal%npc))
    allocate(epsA (nFreq,scal%npr,scal%npc))
  endif
  allocate(epsAux (scal%npr,scal%npc))
  if (freq_dep.eq.0) then
    call read_matrix_d(scal,eps,nmtx,iunit)
  else
    peinf%rank_mtxel=0
    call read_matrix_f(scal, nFreq, nFreq, epsR, nmtx, 1, iunit, advanced=epsA)
    do icolm=1,scal%npc
      do irowm=1,scal%npr
        do i=1,nFreq
          epsA(i,irowm,icolm)=conjg(epsR(i,irowm,icolm))
        enddo
      enddo
    enddo
    do icolm=1,scal%npc
      do irowm=1,scal%npr
        eps(irowm,icolm)=dble(epsR(1,irowm,icolm))
      enddo
    enddo
  endif
  if (peinf%inode.eq.0) then
    call close_file(iunit)
  endif
  if (peinf%inode.eq.0) then
    write(6,'(a,i6)') "     omega num  = ", nFreq
    write(6,'(a,f7.3)') "     omega step = ", dDeltaFreq
    write(6,'(a,f7.3,/)') "     omega brd  = ", dBrdning
  endif
!-----------------------------
! read rho file
  allocate(rhog (ng))
  rhog(:)=0.0d0
  if (peinf%inode.eq.0) then
    write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fnrho)
    iunit=95
    call open_file(unit=iunit,file=trim(fnrho),form='unformatted',status='old')
    rho0=0.0d0
    read(iunit)
    read(iunit) nspin,nvecs
    read(iunit) (FFTgrid(i),i=1,3)
  endif
! (gsm) this is ugly but the only way I see to quickly fix the mess
! with gvec structure and gvec_index subroutine
  gvec%ng=ng
  gvec%FFTgrid=FFTgrid
  call gvec_index(gvec)
  nFFTgridpts=gvec%nFFTgridpts
  if (peinf%inode.eq.0) then
    read(iunit) celvol
    read(iunit) reccelvol, blat, ((bvec(i,j),i=1,3),j=1,3), ((bdot(i,j),i=1,3),j=1,3)
    read(iunit)
    read(iunit)
    read(iunit)
    allocate(gvec_rho (3, nvecs))
    allocate(xcdum (nvecs, nspin))
    read(iunit) nproc_para
    jg = 1
    do i=1,nproc_para
      read(iunit) num_gvec
      read(iunit) ((gvec_rho(k, j), k = 1, 3), j = jg, jg + num_gvec - 1)
      jg = jg + num_gvec
    enddo
    read(iunit) nproc_para
    jg = 1
    do i=1,nproc_para
      read(iunit) num_gvec
      read(iunit) ((xcdum(j, k), j = jg, jg + num_gvec - 1), k = 1, nspin)
      jg = jg + num_gvec
    enddo
    ierr=0
    do j=1,nvecs
      if (gvec_rho(1,j).eq.0.and.gvec_rho(2,j).eq.0.and.gvec_rho(3,j).eq.0) then
        do k=1,nspin
          rho0=rho0+dble(xcdum(j,k))
        enddo
      endif
      iout=((gvec_rho(1,j)+FFTgrid(1)/2)*FFTgrid(2)+gvec_rho(2,j)+FFTgrid(2)/2)* &
        FFTgrid(3)+gvec_rho(3,j)+FFTgrid(3)/2+1
      if (iout.ge.1.and.iout.le.nFFTgridpts) then
        iout=gvec%index_vec(iout)
        if (iout.ge.1.and.iout.le.ng) then
          rhog(iout)=0.0d0
          do k=1,nspin
            rhog(iout)=rhog(iout) + xcdum(j,k)
          enddo
        else
          ierr=1
        endif
      else
        ierr=1
      endif
    enddo
    if(allocated(gvec_rho))then;deallocate(gvec_rho);endif
    if(allocated(xcdum))then;deallocate(xcdum);endif
    call close_file(iunit)
    write(6,'(5x,"cel vol =",e20.12)') celvol
    write(6,'(5x,"rho(0) =",f10.3,/)') rho0
  endif
  if (ierr.ne.0) then
    call die("unknown G vector in file " // trim(fnrho))
  endif
  if (abs(rho0).le.TOL_Zero) then
    call die("cannot find rho(0) in file " // trim(fnrho))
  endif
!-----------------------------
! construct generalized plasmon pole model
  if (peinf%inode.eq.0) then
    write(6,'(1x,"constructing GPP model",/)')
  endif
  epsPP(:,:,:)=(0.0d0,0.0d0)
  wp2=ryd*ryd*16.0d0*PI_D*rho0/celvol
  sflag=.false.
  irows=0
  icols=0
  icurr=0
  do jgp=1,nmtx
    icol=mod(int(((jgp-1)/scal%nbl)+TOL_Small),scal%npcol)
    do jg=1,nmtx
      irow=mod(int(((jg-1)/scal%nbl)+TOL_Small),scal%nprow)
      if (irow.eq.scal%myprow.and.icol.eq.scal%mypcol) then
        icurr=icurr+1
        icolm=int((icurr-1)/scal%npr+TOL_Small)+1
        irowm=mod((icurr-1),scal%npr)+1
        if (jg.eq.ig.and.jgp.eq.igp) then
          sflag=.true.
          irows=irowm
          icols=icolm
        endif
        qg(:)=q(:)+dble(gvec%components(:,isrtx(jg)))
        qgp(:)=q(:)+dble(gvec%components(:,isrtx(jgp)))
        qgqg=dot_product(qg,matmul(bdot,qg))
        qgqgp=dot_product(qg,matmul(bdot,qgp))
        if (abs(qgqg).lt.TOL_Zero) cycle
        gx=gvec%components(1,isrtx(jg))-gvec%components(1,isrtx(jgp))
        gy=gvec%components(2,isrtx(jg))-gvec%components(2,isrtx(jgp))
        gz=gvec%components(3,isrtx(jg))-gvec%components(3,isrtx(jgp))
        iout=((gx+FFTgrid(1)/2)*FFTgrid(2)+gy+FFTgrid(2)/2)* &
          FFTgrid(3)+gz+FFTgrid(3)/2+1
        if (iout.lt.1.or.iout.gt.nFFTgridpts) cycle
        iout=gvec%index_vec(iout)
        if (iout.lt.1.or.iout.gt.ng) cycle
        Omega2=wp2*qgqgp/qgqg*rhog(iout)/rho0
        epsggp=eps(irowm,icolm)
        if (jg.eq.jgp) then
          delta=1.0d0
        else
          delta=0.0d0
        endif
        I_epsggp=delta-epsggp
        if (abs(I_epsggp).lt.TOL_Small) cycle
! Real GPP [PRB 34, 5390 (1986)]
        wtilde2 = Omega2 / I_epsggp
        if (abs(wtilde2) .lt. TOL_Small) cycle
        do i=1,nFreq
          omega=dFreqGrid(i)
          epsPP(i,irowm,icolm)=delta+Omega2/(omega**2-wtilde2-cmplx(0D0,dBrdning,kind=DPC))
        enddo
      endif
    enddo
  enddo
  if (peinf%inode.eq.0) then
    write(6,'(5x,"plasma frequency =",f10.3," eV",/)') sqrt(wp2)
  endif
!-----------------------------
! invert matrices
  if (peinf%inode.eq.0) then
    write(6,'(1x,"inverting matrices",/)')
  endif
  do i=1,nFreq
    epsAux(:,:) = epsPP(i,:,:)
    call zinvert_serial(nmtx, epsAux)
    epsPP(i,:,:) = epsAux(:,:)
  enddo
  if(freq_dep .ne. 0) then
    do i=1,nFreq
      epsAux(:,:) = epsR(i,:,:)
      call zinvert_serial(nmtx, epsAux)
      epsR(i,:,:) = epsAux(:,:)
      epsAux(:,:) = epsA(i,:,:)
      call zinvert_serial(nmtx, epsAux)
      epsA(i,:,:) = epsAux(:,:)
    enddo
  endif
!-----------------------------
! write generalized plasmon pole file
  if (sflag) then
    write(6,'(1x,"writing",1x,a,1x,"file",/)')trim(fngpp)
    call open_file(unit=7,file=trim(fngpp),form='formatted',status='replace')
    do i=1,nFreq
      omega=dFreqGrid(i)
      write(7,100)omega,epsPP(i,irows,icols)
    enddo
    call close_file(7)
  endif
!-----------------------------
! write full frequency files
  if (sflag) then
    if (freq_dep.ne.0) then
      write(6,'(1x,"writing",1x,a,1x,"and",1x,a,1x,"files",/)') &
        trim(fnffr),trim(fnffa)
      call open_file(unit=8,file=trim(fnffr),form='formatted',status='replace')
      call open_file(unit=9,file=trim(fnffa),form='formatted',status='replace')
      do i=1,nFreq
        omega=dFreqGrid(i)
        write(8,100)omega,epsR(i,irows,icols)
        write(9,100)omega,epsA(i,irows,icols)
      enddo
      call close_file(8)
      call close_file(9)
    endif
  endif
!-----------------------------
! deallocate and finish
  if(associated(gvec%components))then;deallocate(gvec%components);nullify(gvec%components);endif
  if(allocated(isrtx))then;deallocate(isrtx);endif
  if(associated(gvec%index_vec))then;deallocate(gvec%index_vec);nullify(gvec%index_vec);endif
  if(allocated(eps))then;deallocate(eps);endif
  if(allocated(epsPP))then;deallocate(epsPP);endif
  if(allocated(dFreqGrid))then;deallocate(dFreqGrid);endif
  if(allocated(dFreqBrd))then;deallocate(dFreqBrd);endif
  if (freq_dep.ne.0) then
    if(allocated(epsR))then;deallocate(epsR);endif
    if(allocated(epsA))then;deallocate(epsA);endif
  endif
  if(allocated(epsAux))then;deallocate(epsAux);endif
  if(allocated(rhog))then;deallocate(rhog);endif
100 format(3f25.15)
end program epsomega
