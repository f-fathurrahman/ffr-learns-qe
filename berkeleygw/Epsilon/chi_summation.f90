!==============================================================================
!
! Modules:
!
! chi_summation_m Last Modified: 04/19/2012 (FHJ)
!
! This module creates the (distributed) polarizability matrix chi by summing
! the (distributed) pol%gme matrices. There are routines that communicate the
! gme matrix using either matrix or element partitioning scheme.
!
!==============================================================================

module chi_summation_m
  use acc_linalg_m
  use algos_epsilon_m
  use global_m
  use blas_m
  use mtxelmultiply_m
  use scalapack_m
  use lin_denominator_m
  use io_utils_m
  use vcoul_generator_m
  use misc_m, only : procmem
  use timing_m, only: timing => epsilon_timing
  implicit none
  private
  public :: create_chi_summator, free_chi_summator,&
    chi_summation_comm_matrix, chi_summation_comm_elements
  !> FHJ: the chi_summator "object"
  type chi_summator_t
    real(DP) :: fact
    !> DVF: below are some temporary buffers needed for the chi summation. They are
    !! described in detail in this comment.
    !! gme = g-vector matrix element
    !! gmetempX where X = n,r,c are temporary arrays for static calculations
    !! n = normalized by the proper factor used in BGW
    !! r = row, meaning the matrix elements for all bands (nv*nc*nk) that the proc owns
    !! for the G-vectors in the rows of chi currently being computed
    !! c = column, the matrix elements for all bands (nv*nc*nk) that the proc owns
    !! for the G`-vectors in the rows of chi currently being computed
    !! the RDyn arrays are needed for full-frequency (FF) calculations, real and complex
    !! r2 is used in FF with matrix communication because you have a different denominator for
    !! each frequency. Only the r2 array (not the r) array is used for element communication
    !! the denominators are built into the gme`s for static calculations
    !! eden arrays hold the energy denominators for FF
    !! chilocal holds the contribution of a given processor to the GG` chunk of epsilon
    !! being computed
    !! chilocal2 holds the MPI reduced GG` chunk of epsilon being computed
    real(DP), allocatable :: chilocal(:,:)
    real(DP), allocatable :: chilocal2(:,:,:)
    complex(DPC), allocatable :: chilocalRDyn(:,:,:)
    complex(DPC), allocatable :: chilocal2RDyn(:,:,:,:)
    complex(DPC), allocatable :: chiRDyntmp(:)
    complex(DPC), allocatable :: chiRDynOld(:,:,:,:)
    real(DP), allocatable :: gmetempr(:,:),gmetempc(:,:)
    real(DP), allocatable :: gmetempn(:)
    complex(DPC), allocatable :: gmeRDyntempn(:)
    complex(DPC), allocatable :: gmeRDyntempr2(:,:,:)
    complex(DPC), allocatable :: gmeRDyntempr3(:,:)
    complex(DPC), allocatable :: gmeRDyntempc(:,:)
    complex(DPC), allocatable :: gmeRDyntempcs(:,:)
    integer, allocatable :: deltaCount(:,:)
    integer, allocatable :: deltaCountReduce(:,:)
  end type chi_summator_t
  public :: chi_summator_t
contains
  subroutine create_chi_summator(this, pol, scal, fact, nspin)
    type(chi_summator_t), intent(INOUT) :: this !<the chi_summator_t object
    type(polarizability), intent(IN) :: pol
    type(scalapack), intent(IN) :: scal
    real(DP), intent(IN) :: fact
    integer, intent(IN) :: nspin
   
    this%fact = fact
    if (pol%freq_dep .eq. 0) then
      allocate(this%gmetempn (pol%nmtx))
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
      allocate(this%gmeRDyntempn (pol%nmtx))
      allocate(this%chiRDyntmp (pol%nfreq_in_group))
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      allocate(this%gmeRDyntempn (pol%nmtx))
      allocate(this%chiRDyntmp (pol%nfreq_in_group))
      allocate(this%chiRDynOld (pol%nfreq_in_group,scal%npr,scal%npc,nspin))
    endif
    if (pol%freq_dep .eq. 3) then
      allocate(this%gmeRDyntempn (pol%nmtx))
      allocate(this%chiRDyntmp (pol%nfreq_in_group))
    endif
   
    return
  end subroutine create_chi_summator
  subroutine free_chi_summator(this, pol)
    type(chi_summator_t), intent(INOUT) :: this !<the chi_summator_t object
    type(polarizability), intent(IN) :: pol
   
    if (pol%freq_dep .eq. 0) then
      if(allocated(this%gmetempn))then;deallocate(this%gmetempn);endif
    endif
    if (pol%freq_dep .eq. 2) then
      if(allocated(this%gmeRDyntempn))then;deallocate(this%gmeRDyntempn);endif
      if(allocated(this%chiRDyntmp))then;deallocate(this%chiRDyntmp);endif
      if ((pol%freq_dep_method .eq. 1)) then
        if(allocated(this%chiRDynOld))then;deallocate(this%chiRDynOld);endif
      endif
    endif
    if (pol%freq_dep .eq. 3) then
      if(allocated(this%gmeRDyntempn))then;deallocate(this%gmeRDyntempn);endif
      if(allocated(this%chiRDyntmp))then;deallocate(this%chiRDyntmp);endif
    endif
   
    return
  end subroutine free_chi_summator
  !-----------------------------------------------------------------------------
  ! GCOMM_MATRIX
  !-----------------------------------------------------------------------------
  !> Create the pol%chi matrix using gcomm_matrix sceheme
  subroutine chi_summation_comm_matrix(this,pol,scal,kp,kpq,vwfn,&
    nst,nrk,indt,pht)
    type(chi_summator_t), intent(INOUT) :: this
    type(polarizability), intent(INOUT) :: pol
    type(scalapack), intent(in) :: scal
    type(kpoints), intent(IN) :: kp,kpq
    type(valence_wfns), intent(IN) :: vwfn
    integer, intent(IN) :: nst(:)
    integer, intent(IN) :: nrk
    integer, intent(INOUT) :: indt(:,:,:)
    real(DP), intent(INOUT) :: pht(:,:,:)
    integer :: ntot_members(pol%nfreq_group)
    integer :: icurr,ntot,ntot2,itot,ntotmax,ifreq_para
    integer :: ipe, ilimit, ii, jj, kk, ll, iv, irk, it, ispin,im,mytot
    integer :: i_myband,tag,irank,tmp_iv,im_proc
    complex(DPC) :: negfact
    real(DP) :: zvalue, cv_energy
    type(cvpair_info) :: cvpair_temp
    integer, allocatable :: tmprowindex(:),tmpcolindex(:)
    real(DP), allocatable :: tmprowph(:),tmpcolph(:), sendbuf(:,:)
    complex(DPC), allocatable :: gmeRDyntempr(:,:)
    complex(DPC), allocatable :: edenDRtemp(:,:)
    ! frequency points for the spectral functions of the polarizability
    integer :: isfreql, isfreqr, nwarn
    real(DP) :: sfreql, sfreqr, wr, wl
    ! Hilbert transform coefficients
    complex(DPC) :: htwR(pol%nfreq,pol%nsfreq), htwA(pol%nfreq,pol%nsfreq)
    complex(DPC) :: c11,c12,c13,c14,c21,c22,c23,c24
    complex(DPC) :: cA11,cA12,cA13,cA14,cA21,cA22,cA23,cA24
    integer :: isf,nf,nsf, request
    real(DP) :: sf1,sf2
    real(DP) :: step1,step2,fqt,eta
    complex(DPC) :: c1,c2,j_dpc,cA1,cA2
    logical :: first_reduce
    ! variables for non-blocking cyclic scheme
    integer :: isend_static, irec_static
    integer :: actual_send, actual_rec
    integer :: nsend_row, nsend_col, nrec_row, nrec_col
    integer :: req_send, tag_send, req_rec, tag_rec
    real(DP), allocatable :: buf_send(:,:,:), buf_rec(:,:,:), buf_temp(:,:,:)
    type(progress_info) :: prog_info
   
    !call alloc_summation_buffers(pol, this%fact)
    if (pol%freq_dep.eq.0) then
      allocate(this%chilocal2 (scal%npr,scal%npc,kp%nspin))
!disabled PARALLEL DO collapse(3)
      do ii = 1, kp%nspin
        do jj = 1, scal%npc
          do kk = 1, scal%npr
            this%chilocal2(kk,jj,ii)=0D0
          enddo
        enddo
      enddo
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
      allocate(this%chilocal2RDyn (scal%npr,scal%npc,pol%nfreq_in_group,kp%nspin))
!disabled PARALLEL DO collapse(4)
      do ii = 1, kp%nspin
        do jj = 1, pol%nfreq_in_group
          do kk = 1, scal%npc
            do ll = 1, scal%npr
              this%chilocal2RDyn(ll,kk,jj,ii)=0D0
            enddo
          enddo
        enddo
      enddo
    endif
! At this moment the spectral method only works for gcomm_elements.
! call die in inread.f90
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      allocate(this%chilocal2RDyn (scal%npr,scal%npc,pol%nfreq_in_group,kp%nspin))
!disabled PARALLEL DO collapse(4)
      do ii = 1, kp%nspin
        do jj = 1, pol%nfreq_in_group
          do kk = 1, scal%npc
            do ll = 1, scal%npr
              this%chilocal2RDyn(ll,kk,jj,ii)=0D0
            enddo
          enddo
        enddo
      enddo
    endif
    if (pol%freq_dep .eq. 3) then
      allocate(this%chilocal2RDyn (scal%npr,scal%npc,pol%nfreq_in_group,kp%nspin))
!disabled PARALLEL DO collapse(4)
      do ii = 1, kp%nspin
        do jj = 1, pol%nfreq_in_group
          do kk = 1, scal%npc
            do ll = 1, scal%npr
              this%chilocal2RDyn(ll,kk,jj,ii)=0D0
            enddo
          enddo
        enddo
      enddo
    endif
    ntot=0
    ntot2=0
    ntot = peinf%nvownactual*peinf%ncownactual
    do irk = 1, nrk
      ntot2=ntot2 + nst(irk)
    enddo
    ntot=ntot*ntot2
    !-------------------------------------------------------------------
    ! Static Be Here
    if (pol%freq_dep .eq. 0) then
      first_reduce = .true.
      call progress_init(prog_info, 'building polarizability matrix', 'processor', &
        peinf%npes)
      if(pol%nonblocking_cyclic) then
        ! initialize buffer for non-blocking cyclic scheme
        ! static process for the communication
        isend_static = MOD(peinf%inode + 1 + peinf%npes, peinf%npes)
        irec_static = MOD(peinf%inode - 1 + peinf%npes, peinf%npes)
        ! allocate my size for the first send
        nsend_row = scal%nprd(peinf%inode+1)
        nsend_col = scal%npcd(peinf%inode+1)
        allocate(buf_send (nsend_row, nsend_col, kp%nspin))
        do ispin = 1 , kp%nspin
          !disabled PARALLEL DO collapse(2)
          do ii = 1, nsend_col
            do jj = 1, nsend_row
              buf_send(jj,ii,ispin) = 0.0d0
            enddo
          enddo
        end do
      end if
      do ipe = 1, peinf%npes
        call progress_step(prog_info)
        if(pol%nonblocking_cyclic) then
          ! calculate the actual process we have to send and we are receiving
          actual_send = MOD(peinf%inode + ipe + peinf%npes, peinf%npes)
          actual_rec = MOD(peinf%inode - ipe + peinf%npes, peinf%npes)
          nrec_row = scal%nprd(actual_rec+1)
          nrec_col = scal%npcd(actual_rec+1)
          ! allocate reciving buffer
          allocate(buf_rec (nrec_row,nrec_col,kp%nspin))
          do ispin = 1 , kp%nspin
            !disabled PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_rec(jj,ii,ispin) = 0.0d0
              enddo
            enddo
          end do
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_comm)
          endif
          ! allocate working array
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_array_alloc)
          endif
          allocate(this%chilocal (nrec_row, nrec_col))
!disabled PARALLEL DO collapse(2)
          do ii = 1, nrec_col
            do jj = 1, nrec_row
              this%chilocal(jj,ii)=0D0
            enddo
          enddo
          allocate(buf_temp (nrec_row, nrec_col, kp%nspin))
          do ispin = 1 , kp%nspin
!disabled PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_temp(jj,ii,ispin) = 0.0d0
              enddo
            enddo
          end do
          allocate(this%gmetempr (nrec_row, ntot))
          allocate(this%gmetempc (nrec_col, ntot))
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_array_alloc)
          endif
        else
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_array_alloc)
          endif
          allocate(this%chilocal (scal%nprd(ipe),scal%npcd(ipe)))
!disabled PARALLEL DO collapse(2)
          do ii = 1, scal%npcd(ipe)
            do jj = 1, scal%nprd(ipe)
              this%chilocal(jj,ii)=0D0
            enddo
          enddo
          allocate(this%gmetempr (scal%nprd(ipe),ntot))
! JRD XXX Should change order here like do in FF case. I think I did...
          allocate(this%gmetempc (scal%npcd(ipe),ntot))
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_array_alloc)
          endif
        end if ! pol%nonblocking_cyclic
        do ispin = 1 , kp%nspin
          if(pol%nonblocking_cyclic) then
            call mtxelmultiply(scal,ntot,nrk,nst,this%fact,vwfn, &
              this%gmetempr,this%gmetempc,this%chilocal,pol%gme,pol,indt,pht,actual_rec+1,ispin)
!disabled PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_temp(jj,ii,ispin) = this%chilocal(jj,ii)
              enddo
            enddo
          else
            call mtxelmultiply(scal,ntot,nrk,nst,this%fact,vwfn, &
              this%gmetempr,this%gmetempc,this%chilocal,pol%gme,pol,indt,pht,ipe,ispin)
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_comm)
          endif
          this%chilocal2(:,:,ispin)=this%chilocal(:,:)
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_comm)
          endif
          ! in case nonblocking_cyclic communication will be finalize outside the spin-loop
          end if ! pol%nonblocking_cyclic
        enddo ! ispin
        if(allocated(this%chilocal))then;deallocate(this%chilocal);endif
        if(allocated(this%gmetempr))then;deallocate(this%gmetempr);endif
        if(allocated(this%gmetempc))then;deallocate(this%gmetempc);endif
        if(pol%nonblocking_cyclic) then
          ! accumulate contribution into receiving buffer
          ! buf_rec(:,:,:) = buf_rec(:,:,:) + buf_temp
          do ispin = 1 , kp%nspin
!disabled PARALLEL DO collapse(2)
            do ii = 1, nrec_col
              do jj = 1, nrec_row
                buf_rec(jj,ii,ispin) = buf_rec(jj,ii,ispin) + buf_temp(jj,ii,ispin)
              enddo
            enddo
          end do
          if(allocated(buf_temp))then;deallocate(buf_temp);endif
          ! copy the messega to the sending buffer for the next cycle
          if(allocated(buf_send))then;deallocate(buf_send);endif
          allocate(buf_send (nrec_row, nrec_col, kp%nspin))
          buf_send = buf_rec
          nsend_row = nrec_row
          nsend_col = nrec_col
          ! deallocate receiving buffer
          if(allocated(buf_rec))then;deallocate(buf_rec);endif
        end if
      enddo ! ipe
      if(pol%nonblocking_cyclic) then
        ! done
        this%chilocal2(:,:,:) = buf_send(:,:,:)
        if(allocated(buf_send))then;deallocate(buf_send);endif
      else
      end if
      call progress_free(prog_info)
      do ispin =1, kp%nspin
        pol%chi(:,:,ispin) = this%chilocal2(:,:,ispin)
      enddo
      if(allocated(this%chilocal2))then;deallocate(this%chilocal2);endif
      !XXXXXXXXXXx
      ! call diagonalize_scalapack(pol%nmtx, scal, pol%chi(:,:,1), 1.0D-3, icurr, this%chilocal)
      !XXXXXXXXXXX
    endif ! pol%freq_dep .eq. 0
    !-------------------------------------
    ! Full Frequency Be Here
    negfact = -1D0*this%fact
    if ( ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) .or. pol%freq_dep .eq. 3 ) then
      if(pol%nfreq_group .gt. 1 .and. peinf%inode .lt. peinf%npes) then
      endif
      call progress_init(prog_info, 'building polarizability matrix', 'processor', &
        peinf%npes_freqgrp)
      do ipe = 1, peinf%npes_freqgrp
        call progress_step(prog_info)
        if (peinf%verb_debug .and. peinf%inode==0) then
          write(6,'(A,i8,6x,A,i8,A)') '### ipe=',ipe,'(npes=',peinf%npes,')'
        endif
        allocate(this%gmeRDyntempr3 (scal%nprd(ipe),ntot))
        allocate(this%gmeRDyntempc (scal%npcd(ipe),ntot))
        allocate(this%chilocalRDyn (scal%nprd(ipe),scal%npcd(ipe),pol%nfreq_in_group))
        allocate(gmeRDyntempr (scal%nprd(ipe),ntot))
        allocate(edenDRtemp (ntot, pol%nfreq_in_group))
!disabled PARALLEL DO collapse(2)
        do ii = 1, pol%nfreq_in_group
          do jj = 1, scal%npcd(ipe)
            do kk = 1, scal%nprd(ipe)
              this%chilocalRDyn(kk,jj,ii)=0D0
            enddo
          enddo
        enddo
        do ispin = 1 , kp%nspin
          itot = 0
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_prep)
          endif
          allocate(tmprowindex (scal%nprd(ipe)))
          allocate(tmpcolindex (scal%npcd(ipe)))
          allocate(tmprowph (scal%nprd(ipe)))
          allocate(tmpcolph (scal%npcd(ipe)))
          do im=1,pol%nfreq_group ! im labels which member of the mtxel comm group are you
            im_proc=peinf%rank_f+1+(im-1)*peinf%npes_freqgrp ! im_proc gives this mtxel comm group member`s global
            do irk = 1, nrk ! proc number
              do it = 1, nst(irk)
                do icurr=1,scal%nprd(ipe)
                  tmprowindex(icurr) = indt(scal%imyrowd(icurr,ipe),it,irk)
                  tmprowph(icurr) = pht(scal%imyrowd(icurr,ipe),it,irk)
                enddo
                do icurr=1,scal%npcd(ipe)
                  tmpcolindex(icurr) = indt(scal%imycold(icurr,ipe),it,irk)
                  tmpcolph(icurr) = pht(scal%imycold(icurr,ipe),it,irk)
                enddo
! JRD XXX - Would rather put OMP here with collapse 2 - but unbalanced... need to loop over tmp_iv = 1, nvown instead
                do iv = 1,vwfn%nband+pol%ncrit
                  tmp_iv = peinf%global_indexv(iv,im_proc)
                  if (peinf%does_it_ownv(iv,im_proc)) then
                    ilimit = peinf%global_ncown(im_proc)
                  else
                    ilimit = 0
                  endif
                  !disabled PARALLEL DO private (mytot,zvalue,jj,icurr,ifreq_para,cvpair_temp,cv_energy)
                  do i_myband = 1, ilimit
                    mytot = itot + i_myband
                    zvalue = pol%edenDyn(peinf%global_indexv(iv,im_proc),i_myband,ispin,irk,im)
                    if(pol%lin_denominator<TOL_Zero) then
                      ! this is when the lin_denominator mode is not active.
                      ! DVF: These factors of (peinf%inode)/peinf%npes*100000 are so that you don`t access the
                      ! array edenDRtemp and other similar arrays for the processors that are not used
                      ! when pol%nfreq_group does not divide the original number of processors in the
                      ! calculation (peinf%npes_orig). This factor makes the loop start at some huge
                      ! number, so the loop will in fact not execute at all
                      do jj=1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
                        ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
                        if (abs(zvalue) .gt. Tol_Zero) then
                          ! DYQ: To get a polarizability consistent with the Tamm-Dancoff Approximation,
                          ! we only include positive frequency terms in the full frequency calculation
                          if(pol%tda) then
                            edenDRtemp(mytot,ifreq_para)= -0.5d0*( &
                              1d0/(zvalue+(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd))
                          else
                            ! JRD XXX - INDEX ORDER here worries me
                            edenDRtemp(mytot,ifreq_para)= -0.5d0*( &
                              1d0/(zvalue-(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd)+ &
                              1d0/(zvalue+(pol%dFreqBrd(jj)+pol%dFreqGrid(jj))/ryd))
                          endif
                        else
                          edenDRtemp(mytot,ifreq_para)= 0D0
                        endif
                      enddo
                    endif
                    do icurr=1,scal%nprd(ipe)
                      gmeRDyntempr(icurr,mytot)=pol%gme(tmprowindex(icurr), &
                        i_myband,tmp_iv,ispin,irk,im) * tmprowph(icurr)
                    enddo
                    !do jj = 1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
                    ! ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
                    ! this%gmeRDyntempr2(:,mytot,ifreq_para)=gmeRDyntempr(:)*edenDRtemp(ifreq_para)
                    !enddo
                    do icurr=1,scal%npcd(ipe)
                      this%gmeRDyntempc(icurr,mytot) = &
                        (pol%gme(tmpcolindex(icurr),i_myband,tmp_iv,ispin,irk,im) * tmpcolph(icurr))
                    enddo
                  enddo ! i_myband
                  !disabled END PARALLEL DO
                  itot = itot+ilimit
                enddo ! iv
              enddo ! it
            enddo ! irk
          enddo ! im
          if(allocated(tmprowindex))then;deallocate(tmprowindex);endif
          if(allocated(tmpcolindex))then;deallocate(tmpcolindex);endif
          if(allocated(tmprowph))then;deallocate(tmprowph);endif
          if(allocated(tmpcolph))then;deallocate(tmpcolph);endif
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_prep)
          endif
          !Do the zgemm`s
          if(ntot > 0) then
            do jj =1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
              ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_gemm)
              endif
              !disabled parallel do private(icurr)
              do mytot = 1, ntot
                do icurr=1,scal%nprd(ipe)
                  this%gmeRDyntempr3(icurr,mytot)=gmeRDyntempr(icurr,mytot)*edenDRtemp(mytot,ifreq_para)
                enddo
              enddo
              call zgemm('n','t',scal%nprd(ipe),scal%npcd(ipe),ntot, &
                negfact,this%gmeRDyntempr3(:,:),scal%nprd(ipe),this%gmeRDyntempc(:,:),scal%npcd(ipe),&
                (0D0,0D0),this%chilocalRDyn(:,:,ifreq_para),scal%nprd(ipe))
              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_gemm)
              endif
            enddo
          endif
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_comm)
          endif
          if(pol%nfreq_group .eq. 1) then
          elseif(pol%nfreq_group .gt. 1 .and. peinf%inode .lt. peinf%npes) then
          endif
          this%chilocal2RDyn(:,:,:,ispin)=this%chilocalRDyn(:,:,:)
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_comm)
          endif
        enddo ! ispin
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_array_alloc)
        endif
        if(allocated(this%chilocalRDyn))then;deallocate(this%chilocalRDyn);endif
        if(allocated(edenDRtemp))then;deallocate(edenDRtemp);endif
        if(allocated(gmeRDyntempr))then;deallocate(gmeRDyntempr);endif
        if(allocated(this%gmeRDyntempr3))then;deallocate(this%gmeRDyntempr3);endif
        if(allocated(this%gmeRDyntempc))then;deallocate(this%gmeRDyntempc);endif
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_array_alloc)
        endif
      enddo ! ipe
      call progress_free(prog_info)
      do ispin =1, kp%nspin
        do jj=1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
          ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
! JRD XXX This assignment is now a waste of time and memory. Should just set pol%chiR(A)Dyn
! directly above
          pol%chiRDyn(:,:,ifreq_para,ispin) = this%chilocal2RDyn(:,:,ifreq_para,ispin)
        enddo ! jj
      enddo ! ispin
      if(allocated(this%chilocal2RDyn))then;deallocate(this%chilocal2RDyn);endif
      !XXXXXXXXXXX
      ! DO jj = 1, pol%nfreq
      ! IF(peinf%inode .eq. 0) WRITE(2000,*) jj
      ! call diagonalize_scalapack(pol%nmtx, scal, pol%chiRDyn(:,:,jj,1), 1.0D-3, icurr, this%chilocal, eigenval)
      ! DEALLOCATE(eigenval)
      ! DEALLOCATE(this%chilocal)
      ! END DO
      !XXXXXXXXXXX
    endif ! (pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0)
    !-------------------------------------
    ! Full Frequency Be Here.
    ! M. Shishkin and G. Kresse, Implementation and performance of the
    ! frequency-dependent GW method within the PAW framework,
    ! PHYSICAL REVIEW B 74, 035101, 2006.
    negfact = -1D0*this%fact
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      if(pol%nfreq_group .gt. 1 .and. peinf%inode .lt. peinf%npes) then
      endif
      ! -------------------------------------------
      ! compute the Hilbert transform coefficients
      ! -------------------------------------------
      j_dpc=(0.0,1.0)
      htwR(:,:)=(0.0,0.0)
      nf=pol%nfreq
      nsf=pol%nsfreq
      do jj=1,nf
        eta=pol%dBrdning/ryd
        fqt=pol%dFreqGrid(jj)/ryd
        do isf=1,nsf
          if (isf==1) then
            c1=(0.0,0.0)
            step1=1.0
          else
            sf1=pol%dSFreqGrid(isf-1)/ryd
            sf2=pol%dSFreqGrid(isf)/ryd
            step1=sf2-sf1
            c11=((sf1-fqt)*(-1.0*j_dpc)-eta)*(atan((fqt-sf2)/(eta))&
              & -atan((fqt-sf1)/(eta)))
            c12=0.50*(sf1-fqt+(-1.0*j_dpc)*eta)*log(((fqt-sf2)**2.0+eta*eta)&
              & /(((fqt-sf1)**2.0+eta*eta)))
            c13=-((sf1+fqt)*j_dpc-eta)*(atan((fqt+sf2)/(eta))&
              & -atan((fqt+sf1)/(eta)))
            c14=0.50*(sf1+fqt+j_dpc*eta)*log(((fqt+sf2)**2.0+eta*eta)&
              & /(((fqt+sf1)**2.0+eta*eta)))
            c1=c11+c12+c13+c14
            c1=c1/step1
          endif
          if (isf==nsf) then
            c2=(0.0,0.0)
            step2=1.0
          else
            sf1=pol%dSFreqGrid(isf)/ryd
            sf2=pol%dSFreqGrid(isf+1)/ryd
            step2=sf2-sf1
            c21=((sf2-fqt)*(-1.0*j_dpc)-eta)*(atan((fqt-sf1)/(eta))&
              & -atan((fqt-sf2)/(eta)))
            c22=0.50*(sf2-fqt+(-1.0*j_dpc)*eta)*log(((fqt-sf1)**2.0+eta*eta)&
              & /(((fqt-sf2)**2.0+eta*eta)))
            c23=-((sf2+fqt)*j_dpc-eta)*(atan((fqt+sf1)/(eta))&
              & -atan((fqt+sf2)/(eta)))
            c24=0.50*(sf2+fqt+j_dpc*eta)*log(((fqt+sf1)**2.0+eta*eta)&
              & /(((fqt+sf2)**2.0+eta*eta)))
            c2=c21+c22+c23+c24
            c2=c2/step2
          endif
          if (isf==1.or.isf==nsf) then
            htwR(jj,isf)=0.5d0*(c1/step1+c2/step2)
          else
            htwR(jj,isf)=1.0d0*(c1+c2)/(step1+step2)
          endif
        enddo
      enddo
      ! ----------------------------------------------------
      ! compute the spectral functions of the polarizability
      ! ----------------------------------------------------
      call progress_init(prog_info, 'building polarizability matrix', 'processor', &
        peinf%npes_freqgrp)
      nwarn=0
      do ipe = 1, peinf%npes_freqgrp
        call progress_step(prog_info)
        if (peinf%verb_debug .and. peinf%inode==0) then
          write(6,'(A,i8,6x,A,i8,A)') '### ipe=',ipe,'(npes=',peinf%npes,')'
        endif
        allocate(gmeRDyntempr (scal%nprd(ipe),1))
        allocate(this%gmeRDyntempr2 (scal%nprd(ipe),ntot,pol%os_nsfreq_para))
        allocate(this%gmeRDyntempc (ntot,scal%npcd(ipe)))
        allocate(this%chilocalRDyn (scal%nprd(ipe),scal%npcd(ipe),pol%os_nsfreq_para))
! JRD XXX should thread if keep COMM elements
        this%chilocalRDyn=0
! JRD XXX if we keep COMM elements we should thread this
        gmeRDyntempr=(0.0,0.0)
        this%gmeRDyntempr2=(0.0,0.0)
        this%gmeRDyntempc=(0.0,0.0)
        do ispin = 1 , kp%nspin
          itot = 0
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_prep)
          endif
          allocate(tmprowindex (scal%nprd(ipe)))
          allocate(tmpcolindex (scal%npcd(ipe)))
          allocate(tmprowph (scal%nprd(ipe)))
          allocate(tmpcolph (scal%npcd(ipe)))
          do im=1,pol%nfreq_group ! im labels which member of the mtxel comm group are you
            im_proc=peinf%rank_f+1+(im-1)*peinf%npes_freqgrp ! im_proc gives this mtxel comm group member`s global
            do irk = 1, nrk ! proc number
              do it = 1, nst(irk)
                do icurr=1,scal%nprd(ipe)
                  tmprowindex(icurr) = indt(scal%imyrowd(icurr,ipe),it,irk)
                  tmprowph(icurr) = pht(scal%imyrowd(icurr,ipe),it,irk)
                enddo
                do icurr=1,scal%npcd(ipe)
                  tmpcolindex(icurr) = indt(scal%imycold(icurr,ipe),it,irk)
                  tmpcolph(icurr) = pht(scal%imycold(icurr,ipe),it,irk)
                enddo
                do iv = 1,vwfn%nband+pol%ncrit
                  tmp_iv = peinf%global_indexv(iv,im_proc)
                  if (peinf%does_it_ownv(iv,im_proc)) then
                    ilimit = peinf%global_ncown(im_proc)
                  else
                    ilimit = 0
                  endif
                  !disabled PARALLEL private (mytot,zvalue,gmeRDyntempr,icurr, &
                  !disabled isfreql,isfreqr,sfreql,sfreqr, &
                  !disabled wl,wr,i_myband,jj)
                  !disabled DO
                  do i_myband = 1, ilimit
                    zvalue = -pol%edenDyn(peinf%global_indexv(iv,im_proc),i_myband,ispin,irk,im)
                    if (abs(zvalue) .gt. Tol_Zero) then
                      mytot = itot + i_myband
                      isfreql=-1
                      do jj=pol%nsfreq,1,-1
                        if ((pol%dSFreqGrid(jj)/ryd)<zvalue) then
                          isfreql=jj
                          EXIT
                        endif
                      enddo
                      if (isfreql.eq.pol%nsfreq) then
                        nwarn=nwarn+1
                        if (nwarn==1.and.peinf%inode.eq.0) then
                          write(0,*) 'WARNING: for accuracy, sfrequency_high_cutoff should be '
                          write(0,*) 'larger than energy difference between highest unoccupied '
                          write(0,*) 'state and lowest occupied state.'
                        endif
                        cycle
                      endif
                      sfreql=pol%dSFreqGrid(isfreql)/ryd
                      isfreqr=isfreql+1
                      sfreqr=pol%dSFreqGrid(isfreqr)/ryd
                      wr= (zvalue-sfreql)/(sfreqr-sfreql)
                      wl= -(zvalue-sfreqr)/(sfreqr-sfreql)
                      do icurr=1,scal%nprd(ipe)
                        gmeRDyntempr(icurr,1)=pol%gme(tmprowindex(icurr), &
                          i_myband,tmp_iv,ispin,irk,im) * tmprowph(icurr)
                      enddo
                      this%gmeRDyntempr2(:,mytot,isfreql)=gmeRDyntempr(:,1)*wl
                      this%gmeRDyntempr2(:,mytot,isfreqr)=gmeRDyntempr(:,1)*wr
                      do icurr=1,scal%npcd(ipe)
                        this%gmeRDyntempc(mytot,icurr) = &
                          (pol%gme(tmpcolindex(icurr),i_myband,tmp_iv,ispin,irk,im) * tmpcolph(icurr))
                      enddo
                    endif
                  enddo ! i_myband
                  !disabled END DO
                  !disabled END PARALLEL
                  itot = itot+ilimit
                enddo ! iv
              enddo ! it
            enddo ! irk
          enddo ! im
          if(allocated(tmprowindex))then;deallocate(tmprowindex);endif
          if(allocated(tmpcolindex))then;deallocate(tmpcolindex);endif
          if(allocated(tmprowph))then;deallocate(tmprowph);endif
          if(allocated(tmpcolph))then;deallocate(tmpcolph);endif
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_prep)
          endif
          !Do the zgemm`s
          if(ntot > 0) then
            do jj =1+peinf%rank_mtxel, pol%nsfreq,pol%nfreq_group
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_gemm)
              endif
              call zgemm('n','n',scal%nprd(ipe),scal%npcd(ipe),ntot, &
                negfact,this%gmeRDyntempr2(:,:,jj),scal%nprd(ipe),this%gmeRDyntempc(:,:),ntot,&
                (0D0,0D0),this%chilocalRDyn(:,:,jj),scal%nprd(ipe))
              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_gemm)
              endif
            enddo
          endif
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_comm)
          endif
          if(pol%nfreq_group .eq. 1) then
          elseif(pol%nfreq_group .gt. 1 .and. peinf%inode .lt. peinf%npes) then
          endif
          this%chilocal2RDyn(:,:,:,ispin)=this%chilocalRDyn(:,:,:)
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_comm)
          endif
        enddo ! ispin
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_array_alloc)
        endif
        if(allocated(this%chilocalRDyn))then;deallocate(this%chilocalRDyn);endif
        if(allocated(gmeRDyntempr))then;deallocate(gmeRDyntempr);endif
        if(allocated(this%gmeRDyntempr2))then;deallocate(this%gmeRDyntempr2);endif
        if(allocated(this%gmeRDyntempc))then;deallocate(this%gmeRDyntempc);endif
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_array_alloc)
        endif
      enddo ! ipe
      call progress_free(prog_info)
      do ispin =1, kp%nspin
        do jj=1+peinf%rank_mtxel,pol%nsfreq,pol%nfreq_group
          pol%chiTDyn(jj,:,:,ispin) = this%chilocal2RDyn(:,:,jj,ispin)
        enddo ! jj
      enddo ! ispin
      if(allocated(this%chilocal2RDyn))then;deallocate(this%chilocal2RDyn);endif
      ! -----------------------------
      ! Hilbert transform
      ! ------------------------------
! JRD XXX - pol%chi... arrays need to be reordered here.
! I think we specifying LDC wrong here - should pol%nfreq_in_group not pol%nfreq right?
      call zgemm('n','n',pol%nfreq,scal%npr*scal%npc*kp%nspin,pol%os_nsfreq_para, &
        (-1D0,0D0),htwR(:,:),pol%nfreq,pol%chiTDyn(:,:,:,:),pol%os_nsfreq_para, &
        (0D0,0D0),this%chiRDynOld(:,:,:,:),pol%nfreq)
      do ispin =1, kp%nspin
        do jj=1,pol%nfreq_in_group
          pol%chiRDyn(:,:,jj,ispin) = this%chiRDynOld(jj,:,:,ispin)
        enddo ! jj
      enddo ! ispin
    endif ! pol%freq_dep.eq.2.and.pol%freq_dep_method.eq.1
    !call free_summation_buffers(pol)
   
    return
  end subroutine chi_summation_comm_matrix
  !-----------------------------------------------------------------------------
  ! GCOMM_ELEMENTS
  !-----------------------------------------------------------------------------
  !> Create the pol%chi matrix using gcomm_elements sceheme
  subroutine chi_summation_comm_elements(this,pol,scal,kp,vwfn,cwfn,&
    nst,nrk,indt,pht)
    type(chi_summator_t), intent(INOUT) :: this
    type(polarizability), intent(INOUT) :: pol
    type(scalapack), intent(in) :: scal
    type(kpoints), intent(IN) :: kp
    type(valence_wfns), intent(IN) :: vwfn
    type(conduction_wfns), intent(IN) :: cwfn
    integer, intent(IN) :: nst(:)
    integer, intent(IN) :: nrk
    integer, intent(INOUT) :: indt(:,:,:)
    real(DP), intent(INOUT) :: pht(:,:,:)
    integer :: icurr,ntot,itot
    integer :: iv, ic, irk, it, ispin
    real(DP) :: temp_gme
    integer, allocatable :: iowna(:)
    integer :: isend
    complex(DPC), allocatable :: edenDRtemp(:)
    real(DP) :: zvalue
    ! frequency points for the spectral functions of the polarizability
    integer :: isfreql, isfreqr
    real(DP) :: sfreql, sfreqr, wr, wl
    ! Hilbert tranform coefficients
    complex(DPC) :: htwR(pol%nfreq,pol%nsfreq), htwA(pol%nfreq,pol%nsfreq)
    complex(DPC) :: c11,c12,c13,c14,c21,c22,c23,c24
    complex(DPC) :: cA11,cA12,cA13,cA14,cA21,cA22,cA23,cA24
    integer :: isf,nf,nsf,ifreq_para
    real(DP) :: sf1,sf2
    real(DP) :: step1,step2,fqt,eta
    complex(DPC) :: c1,c2,j_dpc,cA1,cA2
    integer :: ii, jj
    type(progress_info) :: prog_info
    integer :: nsftot, il, ir, n1, n2, n3, max_nv, nwarn
    integer, allocatable :: count_v(:), ind_v(:,:), ind_sf(:)
   
    !call alloc_summation_buffers(pol, this%fact)
    if(peinf%inode .eq. 0) then
      call timing%start(timing%chi_sum_array_alloc)
    endif
    allocate(iowna (vwfn%nband+pol%ncrit))
    if (pol%freq_dep .eq. 0) then
      allocate(this%chilocal (scal%npr,scal%npc))
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
      allocate(this%chilocalRDyn (scal%npr,scal%npc,pol%nfreq_in_group))
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      allocate(this%chilocalRDyn (scal%npr,scal%npc,pol%os_nsfreq_para))
      allocate(count_v (pol%os_nsfreq_para))
      allocate(ind_sf (pol%os_nsfreq_para))
    endif
    if(peinf%inode .eq. 0) then
      call timing%stop(timing%chi_sum_array_alloc)
    endif
    call logit("Starting chi Sum")
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      ! -------------------------------------------
      ! compute the Hilbert transform coefficients
      ! -------------------------------------------
      if(peinf%inode .eq. 0) then
        call timing%start(timing%chi_sum_ht_nb)
      endif
      j_dpc=(0.0,1.0)
      htwR(:,:)=(0.0,0.0)
      nf=pol%nfreq
      nsf=pol%nsfreq
      do jj=1,nf
        eta=pol%dBrdning/ryd
        fqt=pol%dFreqGrid(jj)/ryd
        do isf=1,nsf
          if (isf==1) then
            c1=(0.0,0.0)
            step1=1.0
          else
            sf1=pol%dSFreqGrid(isf-1)/ryd
            sf2=pol%dSFreqGrid(isf)/ryd
            step1=sf2-sf1
            c11=((sf1-fqt)*(-1.0*j_dpc)-eta)*(atan((fqt-sf2)/(eta))&
              & -atan((fqt-sf1)/(eta)))
            c12=0.50*(sf1-fqt+(-1.0*j_dpc)*eta)*log(((fqt-sf2)**2.0+eta*eta)&
              & /(((fqt-sf1)**2.0+eta*eta)))
            c13=-((sf1+fqt)*j_dpc-eta)*(atan((fqt+sf2)/(eta))&
              & -atan((fqt+sf1)/(eta)))
            c14=0.50*(sf1+fqt+j_dpc*eta)*log(((fqt+sf2)**2.0+eta*eta)&
              & /(((fqt+sf1)**2.0+eta*eta)))
            c1=c11+c12+c13+c14
            c1=c1/step1
          endif
          if (isf==nsf) then
            c2=(0.0,0.0)
            step2=1.0
          else
            sf1=pol%dSFreqGrid(isf)/ryd
            sf2=pol%dSFreqGrid(isf+1)/ryd
            step2=sf2-sf1
            c21=((sf2-fqt)*(-1.0*j_dpc)-eta)*(atan((fqt-sf1)/(eta))&
              & -atan((fqt-sf2)/(eta)))
            c22=0.50*(sf2-fqt+(-1.0*j_dpc)*eta)*log(((fqt-sf1)**2.0+eta*eta)&
              & /(((fqt-sf2)**2.0+eta*eta)))
            c23=-((sf2+fqt)*j_dpc-eta)*(atan((fqt+sf1)/(eta))&
              & -atan((fqt+sf2)/(eta)))
            c24=0.50*(sf2+fqt+j_dpc*eta)*log(((fqt+sf1)**2.0+eta*eta)&
              & /(((fqt+sf2)**2.0+eta*eta)))
            c2=c21+c22+c23+c24
            c2=c2/step2
          endif
          if (isf==1.or.isf==nsf) then
            htwR(jj,isf)=0.5d0*(c1/step1+c2/step2)
          else
            htwR(jj,isf)=1.0d0*(c1+c2)/(step1+step2)
          endif
        enddo
      enddo
      if(peinf%inode .eq. 0) then
        call timing%stop(timing%chi_sum_ht_nb)
      endif
    endif!(pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      pol%chiTDyn=(0.0,0.0)
    endif
    call progress_init(prog_info, 'building polarizability matrix', 'blocks', &
      nrk*kp%nspin*(cwfn%nband-vwfn%nband))
    nwarn=0
    do irk=1,nrk
      if (peinf%verb_debug .and. peinf%inode==0) then
        write(6,'(A,i8,6x,A,i8,A)') '### irk=',irk,'(nrk=',nrk,')'
      endif
      do ispin=1,kp%nspin
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_array_alloc)
        endif
        iowna(:)=1
        ntot=(vwfn%nband+pol%ncrit)*nst(irk)
        if (pol%freq_dep .eq. 0) then
!JRD XXX Should thread
          this%chilocal=0
          allocate(this%gmetempr (scal%npr,ntot))
          allocate(this%gmetempc (ntot,scal%npc))
        endif
        if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
!JRD XXX Should thread
          this%chilocalRDyn=0
          allocate(this%gmeRDyntempr2 (scal%npr,ntot,pol%nfreq_in_group))
          allocate(this%gmeRDyntempc (ntot,scal%npc))
        endif
        if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
!JRD XXX Should thread
          this%chilocalRDyn=(0.0,0.0)
          max_nv=0
          do ic=1,cwfn%nband-vwfn%nband
            count_v=0
            do iv=1,(vwfn%nband+pol%ncrit)
              isend=peinf%global_pairowner(iv,ic)-1
              if (isend .lt. 0) then
                write(0,*) 'Illegal value for mpi proc, isend: ',iv,ic
                call die("internal error in chi_summation")
              endif
              if (isend .eq. peinf%inode) then
                zvalue=-pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1)
                iowna(iv)=iowna(iv)+1
              endif
              if (scal%npr*scal%npc .ne. 0) then
                do it=1, nst(irk)
                  if (abs(zvalue) .gt. Tol_Zero) then
                    isfreql=-1
                    do jj=pol%nsfreq,1,-1
                      if ((pol%dSFreqGrid(jj)/ryd)<zvalue) then
                        isfreql=jj
                        EXIT
                      endif
                    enddo
                    if (isfreql.eq.pol%nsfreq) then
                      nwarn=nwarn+1
                      if (nwarn==1.and.peinf%inode.eq.0) then
                        write(0,*) 'WARNING: for accuracy, sfrequency_high_cutoff should be '
                        write(0,*) 'larger than energy difference between highest unoccupied '
                        write(0,*) 'state and lowest occupied state.'
                      endif
                      cycle
                    endif
                    isfreqr=isfreql+1
                    count_v(isfreql)=count_v(isfreql)+1
                    count_v(isfreqr)=count_v(isfreqr)+1
                  endif
                enddo !it
              endif
            enddo !iv
            if (max_nv<maxval(count_v(:))) then
              max_nv=maxval(count_v(:))
            endif
          enddo !ic
          allocate(this%gmeRDyntempr2 (scal%npr,max_nv,pol%os_nsfreq_para))
          allocate(this%gmeRDyntempc (ntot,scal%npc))
          allocate(this%gmeRDyntempcs (max_nv,scal%npc))
          this%gmeRDyntempr2=(0.0,0.0)
          this%gmeRDyntempc=(0.0,0.0)
          this%gmeRDyntempcs=(0.0,0.0)
          allocate(ind_v (pol%os_nsfreq_para, max_nv))
          this%gmeRDyntempr2=(0.0,0.0)
          iowna(:)=1
        endif!(pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_array_alloc)
        endif
        do ic=1,cwfn%nband-vwfn%nband
          call progress_step(prog_info)
          ! We do two giant loops here for freq_dep cases
          if (pol%freq_dep .eq. 0) then
            itot=0
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_comm)
              endif
              do iv=1,(vwfn%nband+pol%ncrit)
              isend=peinf%global_pairowner(iv,ic)-1
              if (isend .lt. 0) then
                write(0,*) 'Illegal value for mpi proc, isend:',&
                  peinf%inode,iv,ic,peinf%global_pairowner(iv,ic)
                call die("internal error in chi_summation")
              endif
              if (isend .eq. peinf%inode) then
                if (iowna(iv) .gt. peinf%ncownactual) call die('iowna(iv) bigger than ncownactual')
                this%gmetempn(:) = pol%gme(:,iowna(iv),peinf%indexv(iv), &
                  ispin,irk,1) * sqrt(this%fact)
                iowna(iv)=iowna(iv)+1
              endif
              if (scal%npr*scal%npc .ne. 0) then
                do it =1, nst(irk)
                  itot = itot + 1
                  do icurr=1,scal%npr
                    this%gmetempr(icurr,itot)=this%gmetempn(indt(scal%imyrow(icurr),it,irk)) * &
                      pht(scal%imyrow(icurr),it,irk)
                  enddo
                  do icurr=1,scal%npc
                    temp_gme = this%gmetempn(indt(scal%imycol(icurr),it,irk))
                    this%gmetempc(itot,icurr)=(temp_gme * pht(scal%imycol(icurr),it,irk))
                  enddo
                enddo ! it
              endif
            enddo ! iv
            if(peinf%inode .eq. 0) then
              call timing%stop(timing%chi_sum_comm)
            endif
            ! JRD: Using Level3 BLAS here for better performance
            if (scal%npr*scal%npc .ne. 0 .and. ntot > 0) then
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_gemm)
              endif
              call dgemm('n','n',scal%npr,scal%npc,ntot, &
                -1.0d0,this%gmetempr,scal%npr,this%gmetempc,ntot,1.0d0,this%chilocal,scal%npr)
              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_gemm)
              endif
            endif
          endif ! pol%freq_dep .eq. 0
          !---------------------
          ! JRD: Full Frequency Be Here
          if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
            if(peinf%inode .eq. 0) then
              call timing%start(timing%chi_sum_array_alloc)
            endif
            allocate(edenDRtemp (pol%nfreq))
            if(peinf%inode .eq. 0) then
              call timing%stop(timing%chi_sum_array_alloc)
            endif
            itot=0
!JRD XXX Should thread
            do iv=1,(vwfn%nband+pol%ncrit)
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_comm)
              endif
              isend=peinf%global_pairowner(iv,ic)-1
              if (isend .lt. 0) then
                write(0,*) 'Illegal value for mpi proc, isend: ',iv,ic
                call die("internal error in chi_summation")
              endif
              if (isend .eq. peinf%inode) then
                this%gmeRDyntempn(:) = pol%gme(:,iowna(iv),peinf%indexv(iv), &
                  ispin,irk,1) * sqrt(this%fact)
                if (abs(pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1)) .gt. Tol_Zero) then
                  do jj=1,pol%nfreq
                    if (pol%tda) then
                      edenDRtemp(jj)= -0.50d0 * ( 1d0/(pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1)+&
                        (pol%dFreqGrid(jj)+pol%dFreqBrd(jj))/ryd))
                    else
                      edenDRtemp(jj)= -0.50d0 * ( 1d0/(pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1) &
                        -(pol%dFreqGrid(jj)+pol%dFreqBrd(jj))/ryd)+&
                        1d0/(pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1)+&
                        (pol%dFreqGrid(jj)+pol%dFreqBrd(jj))/ryd))
                    endif
                  enddo
                else
                  edenDRtemp(:)=0D0
                endif
                iowna(iv)=iowna(iv)+1
              endif
              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_comm)
              endif
              if (scal%npr*scal%npc .ne. 0) then
                do it =1, nst(irk)
                  if(peinf%inode .eq. 0) then
                    call timing%start(timing%chi_sum_row)
                  endif
                  itot = itot + 1
                  do icurr=1,scal%npr
                    do jj =1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
                      ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
                      this%gmeRDyntempr2(icurr,itot,ifreq_para)= &
                        (this%gmeRDyntempn(indt(scal%imyrow(icurr),it,irk))*pht(scal%imyrow(icurr),it,irk))*edenDRtemp(jj)
                    enddo
                  enddo
                  if(peinf%inode .eq. 0) then
                    call timing%stop(timing%chi_sum_row)
                  endif
                  if(peinf%inode .eq. 0) then
                    call timing%start(timing%chi_sum_column)
                  endif
                  do icurr=1,scal%npc
                    this%gmeRDyntempc(itot,icurr) = &
                      conjg(this%gmeRDyntempn(indt(scal%imycol(icurr),it,irk))*pht(scal%imycol(icurr),it,irk))
                  enddo
                  if(peinf%inode .eq. 0) then
                    call timing%stop(timing%chi_sum_column)
                  endif
                enddo ! it
              endif
            enddo ! iv
            ! JRD: Using Level3 BLAS here for better performance
            if (scal%npr*scal%npc .ne. 0 .and. ntot > 0) then
              do jj =1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
                if(peinf%inode .eq. 0) then
                  call timing%start(timing%chi_sum_gemm)
                endif
                ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
                call zgemm('n','n',scal%npr,scal%npc,ntot,(-1D0,0D0),this%gmeRDyntempr2(:,:,ifreq_para),scal%npr, &
                  this%gmeRDyntempc(:,:),ntot,(1D0,0D0),this%chilocalRDyn(:,:,ifreq_para),scal%npr)
                if(peinf%inode .eq. 0) then
                  call timing%stop(timing%chi_sum_gemm)
                endif
              enddo
            endif
            if(peinf%inode .eq. 0) then
              call timing%start(timing%chi_sum_array_alloc)
            endif
            if(allocated(edenDRtemp))then;deallocate(edenDRtemp);endif
            if(peinf%inode .eq. 0) then
              call timing%stop(timing%chi_sum_array_alloc)
            endif
          endif ! (pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0)
          !---------------------
          ! Full Frequency Be Here(shishkin and Kresse 2006)
          if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
            count_v=0
            ind_v=0
            ind_sf=0
            nsftot=0
            itot=0
            do iv=1,(vwfn%nband+pol%ncrit)
              if(peinf%inode .eq. 0) then
                call timing%start(timing%chi_sum_comm)
              endif
              isend=peinf%global_pairowner(iv,ic)-1
              if (isend .lt. 0) then
                write(0,*) 'Illegal value for mpi proc, isend: ',iv,ic
                call die("internal error in chi_summation")
              endif
              if (isend .eq. peinf%inode) then
                this%gmeRDyntempn(:) = pol%gme(:,iowna(iv),peinf%indexv(iv), &
                  ispin,irk,1) * sqrt(this%fact)
                zvalue=-pol%edenDyn(peinf%indexv(iv),iowna(iv),ispin,irk,1)
                iowna(iv)=iowna(iv)+1
              endif
              if(peinf%inode .eq. 0) then
                call timing%stop(timing%chi_sum_comm)
              endif
              ! compute spectral functions of the polarizability
              if (scal%npr*scal%npc .ne. 0) then
                do it=1, nst(irk)
                  if (abs(zvalue) .gt. Tol_Zero) then
                    if(peinf%inode .eq. 0) then
                      call timing%start(timing%chi_sum_row)
                    endif
                    itot=itot+1
                    isfreql=-1
                    do jj=pol%nsfreq,1,-1
                      if ((pol%dSFreqGrid(jj)/ryd)<zvalue) then
                        isfreql=jj
                        EXIT
                      endif
                    enddo
                    if (isfreql.eq.pol%nsfreq) then
                      cycle
                    endif
                    isfreqr=isfreql+1
                    count_v(isfreql)=count_v(isfreql)+1
                    count_v(isfreqr)=count_v(isfreqr)+1
                    il=count_v(isfreql)
                    ir=count_v(isfreqr)
                    ind_v(isfreql,il)=itot
                    ind_v(isfreqr,ir)=itot
                    sfreql=pol%dSFreqGrid(isfreql)/ryd
                    sfreqr=pol%dSFreqGrid(isfreqr)/ryd
                    wl=-(zvalue-sfreqr)/(sfreqr-sfreql)
                    wr=(zvalue-sfreql)/(sfreqr-sfreql)
                    do icurr=1,scal%npr
                      this%gmeRDyntempr2(icurr,il,isfreql)=this%gmeRDyntempn( &
                        indt(scal%imyrow(icurr),it,irk))*pht(scal%imyrow(icurr),it,irk)*wl
                      this%gmeRDyntempr2(icurr,ir,isfreqr)=this%gmeRDyntempn( &
                        indt(scal%imyrow(icurr),it,irk))*pht(scal%imyrow(icurr),it,irk)*wr
                    enddo
                    if(peinf%inode .eq. 0) then
                      call timing%stop(timing%chi_sum_row)
                    endif
                    if(peinf%inode .eq. 0) then
                      call timing%start(timing%chi_sum_column)
                    endif
                    do icurr=1,scal%npc
                      this%gmeRDyntempc(itot,icurr) = &
                        conjg(this%gmeRDyntempn(indt(scal%imycol(icurr),it,irk))*pht(scal%imycol(icurr),it,irk))
                    enddo
                    if(peinf%inode .eq. 0) then
                      call timing%stop(timing%chi_sum_column)
                    endif
                  endif
                enddo ! it
              endif
            enddo ! iv
            if(peinf%inode .eq. 0) then
              call timing%start(timing%chi_sum_array_alloc)
            endif
            jj=0
            do ii=1+peinf%rank_mtxel,pol%nsfreq,pol%nfreq_group
              if (count_v(ii)>0) then
                jj=jj+1
                ind_sf(jj)=ii
              endif
            enddo
            nsftot=jj
            if(peinf%inode .eq. 0) then
              call timing%stop(timing%chi_sum_array_alloc)
            endif
            if (scal%npr*scal%npc .ne. 0 .and. ntot > 0) then
              do ii=1, nsftot
                if(peinf%inode .eq. 0) then
                  call timing%start(timing%chi_sum_gemm)
                endif
                n1=ind_sf(ii)
                n2=count_v(n1)
                do jj=1,n2
                  n3=ind_v(n1,jj)
                  this%gmeRDyntempcs(jj,:)=this%gmeRDyntempc(n3,:)
                enddo
                call zgemm('n','n',scal%npr,scal%npc,n2,(-1D0,0D0),this%gmeRDyntempr2(:,:,n1),scal%npr, &
                  this%gmeRDyntempcs(:,:),max_nv,(1D0,0D0),this%chilocalRDyn(:,:,n1),scal%npr)
                if(peinf%inode .eq. 0) then
                  call timing%stop(timing%chi_sum_gemm)
                endif
              enddo
            endif
          endif ! (pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)
        enddo ! ic (loop over conduction bands)
        if(peinf%inode .eq. 0) then
          call timing%start(timing%chi_sum_array_alloc)
        endif
        if (pol%freq_dep .eq. 0) then
          if (scal%npr*scal%npc .ne. 0) then
            pol%chi(:,:,ispin) = pol%chi(:,:,ispin) + this%chilocal(:,:)
          endif
          if(allocated(this%gmetempr))then;deallocate(this%gmetempr);endif
          if(allocated(this%gmetempc))then;deallocate(this%gmetempc);endif
        endif
        if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
          if (scal%npr*scal%npc .ne. 0) then
            do jj = 1+peinf%rank_mtxel+(peinf%inode)/peinf%npes*100000,pol%nfreq,pol%nfreq_group
              ifreq_para=(jj+pol%nfreq_group-1)/pol%nfreq_group
! JRD XXX This copy is now a waste of time. Should set pol%chi directly above in the zgemm
              pol%chiRDyn(:,:,ifreq_para,ispin) = pol%chiRDyn(:,:,ifreq_para,ispin) + this%chilocalRDyn(:,:,ifreq_para)
            enddo
          endif
          if(allocated(this%gmeRDyntempr2))then;deallocate(this%gmeRDyntempr2);endif
          if(allocated(this%gmeRDyntempc))then;deallocate(this%gmeRDyntempc);endif
        endif
        if(peinf%inode .eq. 0) then
          call timing%stop(timing%chi_sum_array_alloc)
        endif
        if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
          if(peinf%inode .eq. 0) then
            call timing%start(timing%chi_sum_array_alloc)
          endif
          if (scal%npr*scal%npc .ne. 0) then
            do jj = 1+peinf%rank_mtxel, pol%nsfreq,pol%nfreq_group
              pol%chiTDyn(jj,:,:,ispin) = pol%chiTDyn(jj,:,:,ispin) + this%chilocalRDyn(:,:,jj)
            enddo
          endif
          if(allocated(this%gmeRDyntempr2))then;deallocate(this%gmeRDyntempr2);endif
          if(allocated(this%gmeRDyntempc))then;deallocate(this%gmeRDyntempc);endif
          if(allocated(this%gmeRDyntempcs))then;deallocate(this%gmeRDyntempcs);endif
          if(allocated(ind_v))then;deallocate(ind_v);endif
          if(peinf%inode .eq. 0) then
            call timing%stop(timing%chi_sum_array_alloc)
          endif
        endif
      enddo ! ispin (loop over spins)
    enddo ! irk (loop over k-points in set rk)
    call progress_free(prog_info)
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      ! -------------------------
      ! Hilbert transform
      ! -------------------------
      if(peinf%inode .eq. 0) then
        call timing%start(timing%chi_sum_ht_nb)
      endif
! JRD XXX pol%chi arrays out of order below
      call zgemm('n','n',pol%nfreq,scal%npr*scal%npc*kp%nspin,pol%os_nsfreq_para, &
        (-1D0,0D0),htwR(:,:),pol%nfreq,pol%chiTDyn(:,:,:,:),pol%os_nsfreq_para, &
        (0D0,0D0),this%chiRDynOld(:,:,:,:),pol%nfreq)
      do ispin =1, kp%nspin
        do jj=1,pol%nfreq_in_group
          pol%chiRDyn(:,:,jj,ispin) = this%chiRDynOld(jj,:,:,ispin)
        enddo ! jj
      enddo ! ispin
      if(peinf%inode .eq. 0) then
        call timing%stop(timing%chi_sum_ht_nb)
      endif
    endif !(pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)
    if(peinf%inode .eq. 0) then
      call timing%start(timing%chi_sum_array_alloc)
    endif
    if (pol%freq_dep .eq. 0) then
      if(allocated(this%chilocal))then;deallocate(this%chilocal);endif
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 0 .or. pol%freq_dep_method .eq. 2)) then
      if(allocated(this%chilocalRDyn))then;deallocate(this%chilocalRDyn);endif
    endif
    if ((pol%freq_dep .eq. 2).and.(pol%freq_dep_method .eq. 1)) then
      if(allocated(this%chilocalRDyn))then;deallocate(this%chilocalRDyn);endif
      if(allocated(count_v))then;deallocate(count_v);endif
      if(allocated(ind_sf))then;deallocate(ind_sf);endif
    endif
    if(allocated(iowna))then;deallocate(iowna);endif
    !call free_summation_buffers(pol)
    if(peinf%inode .eq. 0) then
      call timing%stop(timing%chi_sum_array_alloc)
    endif
   
    return
  end subroutine chi_summation_comm_elements
! close the condition ( if defined MPI && defined USESCALAPACK) for subspace method
end module chi_summation_m
