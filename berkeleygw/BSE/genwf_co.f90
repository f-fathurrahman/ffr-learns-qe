!==========================================================================
!
! Routines:
!
! (1) genwf_co() Originally By MLT Last Modified Apr/2016 (FHJ)
!
! input: crys, gvec, kg, kgq, syms, xct, distgwfco, distgwfcoq types
! ik label of k-point in FBZ
! ikq label of k-point in FBZ Q Shifted
!
! output: wfnc conduction wavefunctions at point k
! wfnv valence wavefunctions at point k
!
! FHJ: All the communication is performed with non-blocking MPI routines.
! TODO: Someone should change the parallelization scheme from
! splitting the G-vectors to a distribution of whole bands.
!
!==========================================================================

module genwf_co_m
  use global_m
  use gmap_m
  use input_utils_m
  use misc_m
  use sort_m
  use susymmetries_m
  implicit none
  private
  public :: &
    genwf_co
contains
subroutine genwf_co(crys,gvec,kg,kgq,syms,wfnc,wfnv,xct, &
  ik,ikq,distgwfco,distgwfcoq,workco,workcoq)
  type t_requests
    integer, allocatable :: isort_send(:), isort_recv(:)
    integer, allocatable :: cg_send(:,:,:), cg_recv(:,:,:)
  endtype t_requests
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (grid), intent(in) :: kg
  type (grid), intent(in) :: kgq
  type (symmetry), intent(in) :: syms
  type (wavefunction), intent(out) :: wfnc
  type (wavefunction), intent(out) :: wfnv
  type (xctinfo), intent(in) :: xct
  integer, intent(in) :: ik
  integer, intent(in) :: ikq
  type (tdistgwf), intent(in) :: distgwfco
  type (tdistgwf), intent(in) :: distgwfcoq
  type (work_genwf), intent(inout) :: workco
  type (work_genwf), intent(inout) :: workcoq
  integer :: irk, irkq, nss
  type(t_requests) :: reqs_v, reqs_c
 
! Prepare buffers for the (possibly shifted) val and cond wavefunctions
!-----------------------------------------------------------------------
  if (xct%qflag == 0) then
    irkq = kgq%indr(ikq)
  else
    irkq = kg%indr(ikq)
  endif
  call wfn_prepare_bufs(distgwfcoq, workcoq, wfnv, ikq, irkq, .true.)
  irk = kg%indr(ik)
  call wfn_prepare_bufs(distgwfco, workco, wfnc, ik, irk, .false.)
! Start non-blocking receives and sends
!-----------------------------------------------------------------------
  nss = workco%ns * workco%nspinor
  ! FHJ: We allocate even in the serial case, otherwise we would have to
  ! change all call signatures below.
  allocate(reqs_v%isort_send (peinf%npes))
  allocate(reqs_v%isort_recv (peinf%npes))
  allocate(reqs_c%isort_send (peinf%npes))
  allocate(reqs_c%isort_recv (peinf%npes))
  allocate(reqs_v%cg_send (workcoq%nb, nss, peinf%npes))
  allocate(reqs_v%cg_recv (workcoq%nb, nss, peinf%npes))
  allocate(reqs_c%cg_send (workco%nb, nss, peinf%npes))
  allocate(reqs_c%cg_recv (workco%nb, nss, peinf%npes))
  call wfn_start_recv(distgwfcoq, workcoq, ikq, irkq, reqs_v, .true.)
  call wfn_start_send(distgwfcoq, workcoq, ikq, irkq, reqs_v, .true.)
  call wfn_start_recv(distgwfco, workco, ik, irk, reqs_c, .false.)
  call wfn_start_send(distgwfco, workco, ik, irk, reqs_c, .false.)
! Call gmap and copy valence and conduction WFNs
!-----------------------------------------------------------------------
  if (xct%qflag == 0) then
    call wfn_rotate(workcoq, wfnv, kgq, ikq, reqs_v, .true.)
  else
    call wfn_rotate(workcoq, wfnv, kg, ikq, reqs_v, .true.)
  endif
  call wfn_rotate(workco, wfnc, kg, ik, reqs_c, .false.)
! Make sure the send buffers are ready to be deallocated and exit
!-----------------------------------------------------------------------
  if(allocated(reqs_v%isort_send))then;deallocate(reqs_v%isort_send);endif
  if(allocated(reqs_v%isort_recv))then;deallocate(reqs_v%isort_recv);endif
  if(allocated(reqs_c%isort_send))then;deallocate(reqs_c%isort_send);endif
  if(allocated(reqs_c%isort_recv))then;deallocate(reqs_c%isort_recv);endif
  if(allocated(reqs_v%cg_send))then;deallocate(reqs_v%cg_send);endif
  if(allocated(reqs_v%cg_recv))then;deallocate(reqs_v%cg_recv);endif
  if(allocated(reqs_c%cg_send))then;deallocate(reqs_c%cg_send);endif
  if(allocated(reqs_c%cg_recv))then;deallocate(reqs_c%cg_recv);endif
 
  return
contains
  !> Prepares and allocates work and output WFN buffers work_wfn and wfn.
  subroutine wfn_prepare_bufs(dist_wfn, work_wfn, wfn, ik_, irk_, is_val)
    type (tdistgwf), intent(in) :: dist_wfn
    type (work_genwf), intent(inout) :: work_wfn
    type (wavefunction), intent(out) :: wfn
    integer, intent(in) :: ik_
    integer, intent(in) :: irk_
    logical, intent(in) :: is_val
   
    if (ik_ /= work_wfn%ikold) then
      work_wfn%ng = dist_wfn%ng(irk_)
      if (is_val) then
        work_wfn%nb = dist_wfn%nv
      else
        work_wfn%nb = dist_wfn%nc
      endif
      work_wfn%ns = dist_wfn%ns
      work_wfn%nspinor = dist_wfn%nspinor
      if (work_wfn%ikold /= 0) then
        if(associated(work_wfn%cg))then;deallocate(work_wfn%cg);nullify(work_wfn%cg);endif
        if(associated(work_wfn%ph))then;deallocate(work_wfn%ph);nullify(work_wfn%ph);endif
        if(associated(work_wfn%ind))then;deallocate(work_wfn%ind);nullify(work_wfn%ind);endif
        if(associated(work_wfn%isort))then;deallocate(work_wfn%isort);nullify(work_wfn%isort);endif
      endif
      allocate(work_wfn%cg (work_wfn%ng,work_wfn%nb,work_wfn%ns*work_wfn%nspinor))
      allocate(work_wfn%ind (work_wfn%ng))
      allocate(work_wfn%ph (work_wfn%ng))
      allocate(work_wfn%isort (gvec%ng))
    endif
    wfn%ng = work_wfn%ng
    wfn%nband = work_wfn%nb
    wfn%nspin = work_wfn%ns
    wfn%nspinor = work_wfn%nspinor
    if (work_wfn%ns /= xct%nspin) then
      if (peinf%inode==0) then
        write(0,'(/a,i0,a,/,7x,a,i0,a)') 'ERROR: The given number of spins, ', xct%nspin, &
          ', does not match with', 'the number of spins in the coarse grid, ', work_wfn%ns, '.'
        if (is_val) then
          write(0,'(7x,a/)') 'Error was triggered when processing valence states.'
        else
          write(0,'(7x,a/)') 'Error was triggered when processing conduction states.'
        endif
      endif
      call die('Spin mismatch in genwf_co.wfn_prepare_bufs.', only_root_writes=.true.)
    endif
    allocate(wfn%cg (wfn%ng,wfn%nband,wfn%nspin*wfn%nspinor))
    allocate(wfn%isort (gvec%ng))
   
  end subroutine wfn_prepare_bufs
  !> Call gmap and rotate a wfn stored in work_wfn%cg. This routine will wait
  !! for the non-blocking MPI_IRecv`s to finish. The rotated WFNs will be
  !! stored both in work_wfn and wfn.
  subroutine wfn_rotate(work_wfn, wfn, kg_, ik_, reqs, is_val)
    type (work_genwf), intent(inout) :: work_wfn
    type (wavefunction), intent(inout) :: wfn
    type (grid), intent(in) :: kg_
    integer, intent(in) :: ik_
    type(t_requests), intent(inout) :: reqs
    logical, intent(in) :: is_val
    character(len=7) :: wfn_name
    real(DP) :: ekin(gvec%ng), qk(3)
    integer :: isorti(gvec%ng), ig, is, ispinor, ib
   
    if (is_val) then
      wfn_name = 'wfnv%cg'
    else
      wfn_name = 'wfnc%cg'
    endif
    if (ik_ /= work_wfn%ikold) then
! Compute inverse index array of Fourier components around rk-kpoint
      isorti(:) = 0
      do ig = 1, wfn%ng
        isorti(work_wfn%isort(ig)) = ig
      enddo
! Compute index array of Fourier components around fk-kpoint
      qk(1:3) = kg_%f(1:3,ik_)
      call kinetic_energies(gvec, crys%bdot, ekin, qvec=qk)
      call sortrx(gvec%ng, ekin, work_wfn%isort, gvec=gvec%components)
! Find ind and ph relating wavefunctions in fk to rk-kpoint
      work_wfn%ind(:) = 0
      work_wfn%ph(:) = 0.0d0
      call gmap(gvec, syms, wfn%ng, kg_%itran(ik_), kg_%kg0(:,ik_), &
        work_wfn%isort, isorti, work_wfn%ind, work_wfn%ph, .true.)
! Compute and renormalize valence wavefunctions
      do is = 1, wfn%nspin
        do ib = 1, wfn%nband
          do ispinor = 1, wfn%nspinor
            do ig = 1, wfn%ng
              if (work_wfn%ind(ig) > 0) then
                wfn%cg(ig,ib,is*ispinor) = work_wfn%ph(ig) * &
                  work_wfn%cg(work_wfn%ind(ig), ib, is*ispinor)
              endif
            enddo
          enddo
          call checknorm(wfn_name, ib, ik, wfn%ng, is, wfn%nspinor, wfn%cg(1:wfn%ng,ib,:))
        enddo
      enddo
      work_wfn%cg(:,:,:) = wfn%cg(:,:,:)
    endif ! ik_ /= work_wfn%ikold
    ! In spinor case, we must rotate spinors according to spinor rotation matrix umtrx
    wfn%cg(:,:,:) = work_wfn%cg(:,:,:)
    wfn%isort(:) = work_wfn%isort(:)
    work_wfn%ikold = ik_
   
  end subroutine wfn_rotate
  !> Prepares buffers work_wfn%cg to receive distributed WFNs through
  !! non-blocking MPI calls (MPI_IRecv).
  subroutine wfn_start_recv(dist_wfn, work_wfn, ik_, irk_, reqs, is_val)
    type (tdistgwf), intent(in) :: dist_wfn
    type (work_genwf), intent(inout) :: work_wfn
    integer, intent(in) :: ik_
    integer, intent(in) :: irk_
    type(t_requests), intent(inout) :: reqs
    logical, intent(in) :: is_val
    integer :: is, ib, ng_ik, ng_offset, ng_ipe, tag, nb, ipe
   
   
  end subroutine wfn_start_recv
  !> Prepares buffers dist_wfn%zv/zc to send distributed WFNs through
  !! non-blocking MPI calls (MPI_ISend).
  subroutine wfn_start_send(dist_wfn, work_wfn, ik_, irk_, reqs, is_val)
    type (tdistgwf), intent(in) :: dist_wfn
    type (work_genwf), intent(inout) :: work_wfn
    integer, intent(in) :: ik_
    integer, intent(in) :: irk_
    type(t_requests), intent(inout) :: reqs
    logical, intent(in) :: is_val
    integer, dimension(peinf%npes) :: dist_irk, dist_ik, dist_ikold
    integer :: is, ib, ng_ik, ng_offset, ng_ipe, tag, nb, ipe
    !integer :: ikr
   
    dist_irk(:) = 0
    dist_irk(peinf%inode+1) = irk_
    dist_ik(:) = 0
    dist_ik(peinf%inode+1) = ik_
    dist_ikold(:) = 0
    dist_ikold(peinf%inode+1) = work_wfn%ikold
    nb = work_wfn%nb
    if (ik_ /= work_wfn%ikold) then
      ng_ik = work_wfn%ng ! FHJ: Same as dist_wfn%ng(dist_irk(peinf%inode+1))
      work_wfn%cg(:,:,:) = 0.0d0
      work_wfn%isort(:) = 0
      work_wfn%isort(1:ng_ik) = dist_wfn%isort(1:ng_ik, irk_)
      if (is_val) then
        work_wfn%cg(1:ng_ik, 1:nb, 1:nss) = dist_wfn%zv(1:ng_ik, 1:nb, 1:nss, irk_)
      else
        work_wfn%cg(1:ng_ik, 1:nb, 1:nss) = dist_wfn%zc(1:ng_ik, 1:nb, 1:nss, irk_)
      endif
    endif
   
  end subroutine wfn_start_send
end subroutine genwf_co
end module genwf_co_m
