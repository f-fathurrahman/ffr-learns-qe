!==========================================================================
!
! Routines:
!
! (1) input_co_q() By ? Last Modified 4/6/2015 (DYQ)
!
! input: crys, gvec, syms, xct, kg, flagbz types
!
! output: indexq_co, kgq, distgwfcoq types
!
! Reads in the coarse grid valence wavefunctions from file WFNq_co
! and distributes them between processors. The k-point grid is stored in kgq.
!
! WARNING: Since this routine is not used by a working part of the code,
! it has not been tested since implementation of new wfn format. --DAS
! eqp corrections, Fermi level check, scissor shift probably needed.
! UPDATED WARNING: This routine is now used by a working part of the code
! and is compatible with the new WFN format. Fermi level check and scissor
! shift still needed. --DYQ
!==========================================================================

module input_co_q_m
  use checkbz_m
  use eqpcor_m
  use fullbz_m
  use global_m
  use input_utils_m
  use io_utils_m
  use misc_m
  use scissors_m
  use wfn_rho_vxc_io_m
  implicit none
  private
  public :: &
    input_co_q
contains
subroutine input_co_q(kp,kp_co,crys,gvec,kg,kgq,syms,xct,flagbz,distgwfcoq,eqp)
  type (kpoints), intent(inout) :: kp_co
  type (kpoints), intent(in) :: kp
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (symmetry), intent(in) :: syms
  type (xctinfo), intent(inout) :: xct
  type (grid), intent(in) :: kg
  integer, intent(in) :: flagbz
  type (tdistgwf), intent(out) :: distgwfcoq
  type (grid), intent(inout) :: kgq
  type (eqpinfo), intent(inout) :: eqp
  type (crystal) :: crysq
  type (gspace) :: gvecq
  type (symmetry) :: symsq
  type(kpoints) :: kpq
  type (wavefunction) :: wfnv
  character :: filenamev*20
  character :: fncor*32
  integer :: iunit_v
  integer :: irk,irks,ik,ikq,umk
  integer :: ii,jj,kk,is,isp,minband,maxband
  real(DP) :: qq_temp(3),delta,diffvol,vcell
  integer, allocatable :: indxkq(:)
  real(DP), allocatable :: cg(:,:), cgarray(:)
  character(len=3) :: sheader
  integer :: iflavor
  type(gspace) :: gvec_kpt
  logical :: skip_checkbz, broken_degeneracy
  type(progress_info) :: prog_info
  ! WFN HDF5 stuff
  integer, allocatable :: gvec_kpt_components_all(:,:)
  integer, allocatable :: ib_size_array(:)
  integer :: istart, ngktot
  integer :: iband_min, iband_max, tot_nbands, band_block, min_band_block
  integer :: ib_first, ib_last, ib_size, ib_size_max
  integer :: ipe, npes_hdf5, ii_start, ii_end, ii_loc
  real(DP), allocatable :: my_cg_all(:,:,:), ipe_cg_all(:,:,:)
  integer :: error
 
!-------------------------------------------------------------------------
! Read the Wavefunction header and check that it`s compatible with WFN_co
!
  sheader = 'WFN'
  iflavor = 0
  if ( xct%use_wfn_hdf5 ) then
    if (peinf%inode==0) write(6,'(a)') ' Reading header of WFNq_co.h5'
  else
    if (peinf%inode == 0) call open_file(26,file='WFNq_co',form='unformatted',status='old')
    call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, symsq, crysq, dont_warn_kgrid=xct%patched_sampling_co)
  end if
  call check_trunc_kpts(xct%icutv, kpq)
  call check_header('WFN_fi', kp, gvec, syms, crys, 'WFNq_co', kpq, gvecq, symsq, crysq, is_wfn = .true.)
  allocate(gvecq%components (3, gvecq%ng))
  if( xct%use_wfn_hdf5 ) then
  else
    call read_binary_gvectors(26, gvecq%ng, gvecq%ng, gvecq%components)
  end if
  if(associated(gvecq%components))then;deallocate(gvecq%components);nullify(gvecq%components);endif
  if(any(kp_co%kgrid(1:3) /= kpq%kgrid(1:3))) then
    if(peinf%inode == 0) then
      write(0,*) 'WFN_co  kgrid = ', kp_co%kgrid(1:3)
      write(0,*) 'WFNq_co kgrid = ', kpq%kgrid(1:3)
    endif
    call die('kgrids for WFN_co and WFNq_co must be the same', only_root_writes = .true.)
  endif
  call get_volume(vcell,crysq%bdot)
  diffvol=abs(crysq%celvol-vcell)
  if (diffvol.gt.TOL_Small) then
    call die('volume mismatch', only_root_writes = .true.)
  endif
  call assess_degeneracies(kpq, kpq%el(kpq%mnband, :, :), kpq%mnband - 1, xct%efermi, TOL_Degeneracy)
  if(any(kpq%ifmax(:,:) == 0)) &
    call die("BSE codes cannot handle a system where some k-points have no occupied bands.", only_root_writes = .true.)
  kpq%nvband=minval(kpq%ifmax(:,:)-kpq%ifmin(:,:))+1
  kpq%ncband=kpq%mnband-maxval(kpq%ifmax(:,:))
  kgq%nr=kpq%nrk
  allocate(kgq%r (3,kgq%nr))
  kgq%r(1:3,1:kgq%nr)=kpq%rk(1:3,1:kpq%nrk)
  allocate(indxkq (kgq%nr))
  do ii=1,kgq%nr
    indxkq(ii) = ii
  enddo
  ! Check that finite Q is the same as the difference between kgq and kg
  ! DYQ TODO: Fix symmetries before commiting. For now, no symmetries for WFNq_co
  if (.true.) then
    call fullbz(crysq,symsq,kgq,1,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  else
    call fullbz(crysq,symsq,kgq,symsq%ntran,skip_checkbz,wigner_seitz=.true.,paranoid=.true.)
  endif
  if (flagbz.eq.0.and.peinf%inode.eq.0) write(6,901)
  if (flagbz.eq.1.and.peinf%inode.eq.0) write(6,902)
901 format(1x,'Using symmetries to expand the shifted coarse-grid sampling')
902 format(1x,'No symmetries used in the shifted coarse-grid sampling')
!----------------------------------------------------------------------
! Find mapping between kgq and kg
!
  do ik=1,kgq%nf
    ikq = 0
    delta = 0.1d0
    do while ((delta .gt. TOL_Small) .and. (ikq.lt.kg%nf))
      ikq = ikq+1
      qq_temp(:) = kgq%f(:,ikq) - kg%f(:,ik) - xct%finiteq(:)
      do jj=1,3
        qq_temp(jj) = qq_temp(jj) - anint( qq_temp(jj) )
      enddo
      delta=sqrt((qq_temp(1))**2+(qq_temp(2))**2+(qq_temp(3))**2)
    enddo
    if (delta .gt. TOL_Small) then
      if(peinf%inode.eq.0) then
        write(0,*) '  Could not find point equivalent to ', (kg%f(ii,ik),ii=1,3)
      endif
      call die('Finite momentum not commensurate with kgrid of WFN_co',only_root_writes = .true.)
    endif
    xct%indexq(ik)=ikq
    ! kg%f and kgq%f may differ by a lattice vector near the zone edge
    do jj=1,3
      umk = nint(kg%f(jj,ik) - kgq%f(jj,ikq) + xct%finiteq(jj))
      kgq%f(jj,ikq) = kgq%f(jj,ikq) + dble(umk) !kgq(indexq(ik)) = kg(ikq)
      kgq%kg0(jj,ikq) = kgq%kg0(jj,ikq) + umk
    enddo
  enddo
!-----------------------------------------------------------------------
! If it exists, read eqp_co_q.dat for interpolation
  allocate(eqp%evshift_co_q (xct%nvb_co,kpq%nrk,kpq%nspin))
  eqp%evshift_co_q=0D0
  allocate(eqp%ecshift_co_q (xct%ncb_co,kpq%nrk,kpq%nspin))
  eqp%ecshift_co_q=0D0
  if (xct%eqp_co_q_corrections) then
    fncor = 'eqp_co_q.dat'
    allocate(kpq%elda (kpq%mnband, kpq%nrk, kpq%nspin))
    kpq%el(:,:,:) = kpq%el(:,:,:) - xct%avgpot / ryd
    kpq%elda(:,:,:) = kpq%el(:,:,:)
    minband = minval(kpq%ifmax(:,:)-xct%nvb_co+1)
    maxband = maxval(kpq%ifmax(:,:)+xct%ncb_co)
    call eqpcor(fncor,peinf%inode,peinf%npes,kpq, &
      minband,maxband,0,0,kpq%el,eqp%evshift_co_q,eqp%ecshift_co_q,1,0)
    !TODO: implement scissors schift
    ! now we call again to initialize the eqp arrays
    if(xct%eqp_co_corrections) then
      call eqpcor(fncor,peinf%inode,peinf%npes,kpq,0,0, &
        xct%nvb_co,xct%ncb_co,kp_co%el,eqp%evshift_co_q,eqp%ecshift_co_q,1,2,dont_write=.true.)
    endif
    if(any(kpq%ifmax(:,:) == 0)) &
      call die("BSE codes cannot handle a system where some k-points have no occupied bands.", only_root_writes = .true.)
    kpq%nvband=minval(kpq%ifmax(:,:)-kpq%ifmin(:,:))+1
    kpq%ncband=kpq%mnband-maxval(kpq%ifmax(:,:))
  endif !eqp_co_q_correction
!-----------------------------------------------------------------------
! Initialization of distributed wavefunctions
  distgwfcoq%nk=kgq%nr
  distgwfcoq%ngm=kpq%ngkmax
  distgwfcoq%ns=kpq%nspin
  distgwfcoq%nspinor=kpq%nspinor
  distgwfcoq%nv=xct%nvb_co
  ! FHJ: Use standard BLACS distribution for G-vectors
  distgwfcoq%block_sz = ((distgwfcoq%ngm+peinf%npes-1)/(peinf%npes))
  ! ngl = local number of G-vectors that I own.
  distgwfcoq%ngl = NUMROC(distgwfcoq%ngm, distgwfcoq%block_sz, peinf%inode, 0, peinf%npes)
  ! Local to global index translation: ig_g = ig_l + tgl
  distgwfcoq%tgl = distgwfcoq%block_sz * peinf%inode
  allocate(distgwfcoq%ng (distgwfcoq%nk))
  allocate(distgwfcoq%isort (distgwfcoq%ngl,distgwfcoq%nk))
  allocate(distgwfcoq%zv (distgwfcoq%ngl,distgwfcoq%nv,distgwfcoq%ns*distgwfcoq%nspinor,distgwfcoq%nk))
  !allocate(distgwfcoq%zc (distgwfcoq%ngl,distgwfcoq%nc,distgwfcoq%ns*distgwfcoq%nspinor,distgwfcoq%nk))
  distgwfcoq%ng(:)=0
  distgwfcoq%isort(:,:)=0
  distgwfcoq%zv(:,:,:,:)=0.0d0
  !distgwfcoq%zc(:,:,:,:)=0.0d0
!-----------------------------------------------------------------------
! Read the wavefunctions and distribute
  allocate(wfnv%isort (gvec%ng))
  wfnv%nband=xct%nvb_co
  wfnv%nspin=kpq%nspin
  wfnv%nspinor=kpq%nspinor
  if ( xct%use_wfn_hdf5 ) then
    if (peinf%inode==0) write(6,*)
    if (peinf%inode==0) write(6,'(a)') ' Reading HDF5 wavefuntion (WFNq_co.h5)'
    ngktot = SUM(kpq%ngk)
    allocate(gvec_kpt_components_all (3,ngktot))
    ! only valence
    iband_min = MINVAL(kpq%ifmax(:,:)) - xct%nvb_co + 1
    iband_min = MAX(iband_min, 1)
    iband_max = MAXVAL(kpq%ifmax(:,:))
    iband_max = MIN(iband_max, kpq%mnband)
    tot_nbands = iband_max - iband_min + 1
    band_block = (tot_nbands + peinf%npes - 1) / peinf%npes
    min_band_block = 128.0D+00 / ( dble(ngktot) * dble(kpq%nspin*kpq%nspinor) * 16.0D+00 /1024.0D+00/1024D+00 )
    min_band_block = MAX(min_band_block, 1)
    if ( xct%wfn_hdf5_min_band_block > 0 ) min_band_block = xct%wfn_hdf5_min_band_block
    if ( min_band_block > band_block ) then
      band_block = min_band_block
    end if
    ib_first = iband_min + band_block * peinf%inode
    ib_last = min(ib_first + band_block - 1, iband_max)
    ib_size = ib_last - ib_first + 1
    if ( ib_size < 1 ) then
      ! don`t read
      ib_first = -1
      ib_last = -1
      ib_size = 0
    end if
    allocate(ib_size_array (peinf%npes))
    ib_size_array = 0
    ib_size_array(peinf%inode+1) = ib_size
    ib_size_max = MAXVAL( ib_size_array )
    allocate(my_cg_all (ngktot, kpq%nspin*kpq%nspinor, MAX(ib_size,1)))
    allocate(ipe_cg_all (ngktot, kpq%nspin*kpq%nspinor, ib_size_max))
    my_cg_all = 0.0d0
    allocate(wfnv%cg (MAXVAL(kpq%ngk),wfnv%nband,wfnv%nspin*wfnv%nspinor))
    allocate(cgarray (MAXVAL(kpq%ngk)))
  end if
  if ( xct%use_wfn_hdf5 ) then
    call progress_init(prog_info, 'reading wavefunctions (WFNq_co.h5)', 'MPI task', peinf%npes)
  else
    call progress_init(prog_info, 'reading wavefunctions (WFNq_co)', 'k-point', kpq%nrk)
  end if
  npes_hdf5 = 0
  if ( xct%use_wfn_hdf5 ) npes_hdf5 = peinf%npes-1
do ipe = 0, npes_hdf5
  if ( xct%use_wfn_hdf5 ) then
    call progress_step(prog_info, ipe+1)
    if ( ib_size_array(ipe+1) > 0) then
     ipe_cg_all = 0.0d0
     if ( ipe == peinf%inode ) then
       ipe_cg_all(:, :, 1:ib_size) = my_cg_all(:, :, 1:ib_size)
     end if
    else
      ! no elements on this MPI task, cycle
      cycle
    end if
  end if
  istart = 1
  do irk=1,kpq%nrk
    if ( .not. xct%use_wfn_hdf5 ) call progress_step(prog_info, irk)
    irks = 0
    do ii=1,kgq%nr
      if (indxkq(ii) == irk) then
        irks=ii
        exit
      endif
    enddo
    allocate(gvec_kpt%components (3, kpq%ngk(irk)))
    if ( xct%use_wfn_hdf5 ) then
      gvec_kpt%components(:,:) = gvec_kpt_components_all(1:3, istart:istart+kpq%ngk(irk)-1)
    else
      call read_binary_gvectors(26, kpq%ngk(irk), kpq%ngk(irk), gvec_kpt%components)
    end if
    allocate(cg (kpq%ngk(irk),kpq%nspin*kpq%nspinor))
    if(irks > 0) then
      do ii = 1, kpq%ngk(irk)
        call findvector(wfnv%isort(ii), gvec_kpt%components(:, ii), gvec)
        if (wfnv%isort(ii) == 0) call die('Could not find g-vector.')
      enddo
      wfnv%ng=kpq%ngk(irk)
      if ( xct%use_wfn_hdf5 ) then
        cgarray = 0.0d0
      else
        if(peinf%inode == 0) then
          allocate(wfnv%cg (wfnv%ng,wfnv%nband,wfnv%nspin*wfnv%nspinor))
          allocate(cgarray (kpq%ngk(irk)))
        endif
      end if
    endif
    if(associated(gvec_kpt%components))then;deallocate(gvec_kpt%components);nullify(gvec_kpt%components);endif
! Loop over the bands
    !XXX do ii=1,kpq%mnband
    ii_start = 1
    ii_end = kpq%mnband
    if ( xct%use_wfn_hdf5 ) then
      ii_start = iband_min + band_block * ipe
      ii_end = min(ii_start+ band_block - 1, iband_max)
    end if
    ii_loc = 0
    do ii = ii_start, ii_end
      ii_loc = ii_loc + 1
      if ( xct%use_wfn_hdf5 ) then
        cg = 0.0d0
        cg (:,:) = ipe_cg_all(istart:istart+kpq%ngk(irk)-1, 1:kpq%nspin*kpq%nspinor, ii_loc)
      else
        call read_binary_data(26, kpq%ngk(irk), kpq%ngk(irk), kpq%nspin*kpq%nspinor, cg, bcast=.false.)
      end if
      if(irks == 0) cycle
      if(peinf%inode == 0 .or. xct%use_wfn_hdf5) then ! set up wfnv on root only
        do is=1, kpq%nspin
          ! Only need valence wave functions
          if (ii .gt. (kpq%ifmax(irk,is)-xct%nvb_co) .and. ii .le. (kpq%ifmax(irk,is))) then
            do isp=1, kpq%nspinor
              do kk=1, kpq%ngk(irk)
                cgarray(kk)=cg(kk, is*isp)
              enddo
              if (peinf%verb_debug) then
                write(6,'(a,3(1x,i0),2(1x,f18.13))') 'input_co_q', irks, ii, is*isp, cgarray(1)
              endif
              wfnv%cg(1:wfnv%ng,kpq%ifmax(irk,is)-ii+1,is*isp)=cgarray(1:wfnv%ng)
            enddo
            call checknorm('WFNq_co',ii,irks,kpq%ngk(irk),is,kpq%nspinor,cg(:,:))
          end if
        end do ! loop over spins
      endif ! peinf%inode=0
    enddo ! ii (loop over bands)
    if(allocated(cg))then;deallocate(cg);endif
    if ( .not. xct%use_wfn_hdf5 ) then
      if(peinf%inode == 0) then
        if(allocated(cgarray))then;deallocate(cgarray);endif
      endif
    end if ! not HDF5
    distgwfcoq%ng(irks)=wfnv%ng
    do ii=1,distgwfcoq%ngl
      if (ii+distgwfcoq%tgl.le.wfnv%ng) &
        distgwfcoq%isort(ii,irks)=wfnv%isort(ii+distgwfcoq%tgl)
    enddo
    ! copy wfnv%cg to distgwfcoq%zv for g-vectors owned by current proc
    if ( xct%use_wfn_hdf5 ) then
      do is=1, kp_co%nspin
        do isp=1, kp_co%nspinor
          kk = is * isp
          do jj = ii_start, MIN(ii_end, kpq%ifmax(irk,is))
            if (jj .gt. (kpq%ifmax(irk,is)-xct%nvb_co) .and. jj .le. (kpq%ifmax(irk,is))) then
              do ii=1,distgwfcoq%ngl
                if (ii+distgwfcoq%tgl.le.wfnv%ng) then
                  distgwfcoq%zv(ii, kpq%ifmax(irk,is)-jj+1, kk,irks)=wfnv%cg(ii+distgwfcoq%tgl,kpq%ifmax(irk,is)-jj+1,kk)
                endif
              enddo
            end if
          end do
        end do
      end do
    else
      do kk=1,distgwfcoq%ns*distgwfcoq%nspinor
        do jj=1,distgwfcoq%nv
          do ii=1,distgwfcoq%ngl
            if (ii+distgwfcoq%tgl.le.wfnv%ng) then
              distgwfcoq%zv(ii,jj,kk,irks)=wfnv%cg(ii+distgwfcoq%tgl,jj,kk)
            endif
          enddo
        enddo
      enddo
    end if
    if ( .not. xct%use_wfn_hdf5 ) then
      if(associated(wfnv%cg))then;deallocate(wfnv%cg);nullify(wfnv%cg);endif
    end if
    istart = istart + kpq%ngk(irk)
  enddo ! irk (loop over k-points)
end do ! ipe loop for HDF5 case
  call progress_free(prog_info)
  if ( xct%use_wfn_hdf5 ) then
    if(allocated(gvec_kpt_components_all))then;deallocate(gvec_kpt_components_all);endif
    if(allocated(ib_size_array))then;deallocate(ib_size_array);endif
    if(allocated(my_cg_all))then;deallocate(my_cg_all);endif
    if(allocated(ipe_cg_all))then;deallocate(ipe_cg_all);endif
    if(allocated(cgarray))then;deallocate(cgarray);endif
    if(associated(wfnv%cg))then;deallocate(wfnv%cg);nullify(wfnv%cg);endif
  end if
  if(associated(wfnv%isort))then;deallocate(wfnv%isort);nullify(wfnv%isort);endif
  if(allocated(indxkq))then;deallocate(indxkq);endif
  if (peinf%inode.eq.0) then
    write(6,'(/,1x,a)') 'Coarse-grid wavefunctions read from file WFNq_co:'
    write(6,'(1x,a,i0)') '- Number of k-points in irreducible BZ: ', kgq%nr
    write(6,'(1x,a,i0)') '- Number of k-points in full BZ: ', kgq%nf
    if (peinf%verb_high) then
      write(6,'(1x,a)') '- Listing all k-points:'
      write(6,'(1(2x,3(1x,f10.6)))') (kgq%r(:,jj), jj=1,kgq%nr)
    endif
    if ( .not. xct%use_wfn_hdf5 ) call close_file(26)
  endif ! node 0
  if(allocated(kpq%rk))then;deallocate(kpq%rk);endif
  if(allocated(kpq%ifmin))then;deallocate(kpq%ifmin);endif
  if(allocated(kpq%ifmax))then;deallocate(kpq%ifmax);endif
  if(allocated(kpq%el))then;deallocate(kpq%el);endif
  if(allocated(kp_co%rk))then;deallocate(kp_co%rk);endif
  if(allocated(kp_co%ifmin))then;deallocate(kp_co%ifmin);endif
  if(allocated(kp_co%ifmax))then;deallocate(kp_co%ifmax);endif
  if(allocated(kp_co%el))then;deallocate(kp_co%el);endif
 
  return
end subroutine input_co_q
end module input_co_q_m
