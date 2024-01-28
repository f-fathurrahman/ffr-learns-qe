!==========================================================================
!
! Routines:
!
! (1) mtxel_vxc() Last Modified: 5/12/2008 (JRD)
!
! Calculates matrix elements of the DFT exchange-correlation potential
! sig%vxc and puts the results into alda:
!
! alda(in,1) = <nk|Vxc(r)|nk> with n = sig%diag(in)
! or = <nk+phonq|dVxc(r)|mk> for EP
!
!==========================================================================
module mtxel_vxc_m
  use fftw_m
  use global_m
  use misc_m
  implicit none
  private
  public :: &
    mtxel_vxc
contains
subroutine mtxel_vxc(kp,gvec,sig,wfnk,wfnkoff,alda,vxc_type,wfnk_phonq,wfnk_phonq_off)
  type (kpoints), intent(in) :: kp
  type (gspace), intent(in) :: gvec
  type (siginfo), intent(in) :: sig
  type (wfnkstates), intent(inout) :: wfnk, wfnkoff
  real(DP), intent(out) :: alda(sig%ndiag+sig%noffdiag,sig%nspin)
  integer, intent(in) :: vxc_type
  type (wfnkstates), intent(inout), optional :: wfnk_phonq, wfnk_phonq_off
  integer :: in,ioff,ispin
  integer :: jsp,jspmin,jspmax
  integer :: ig,iso,jso,ik,ii,jj,ij,js,source,dest,tag
  character (len=100) :: tmpstr
  real(DP), allocatable :: alda2(:,:)
  complex(DPC), allocatable :: &
    fftbox1(:,:,:,:),fftbox2(:,:,:,:),fftbox3(:,:,:,:)
  integer, dimension(:), allocatable :: gvecindex
  integer, dimension(3) :: Nfft
  integer :: ix, iy, iz
  real(DP) :: scale
  real(DP), allocatable :: vxctemp(:,:)
  real(DP) :: temp
!-------------------- Begin Routine --------------------------------------------
 
! Consistency check
  if ((.not.sig%elph) .and. (vxc_type .eq. 3)) then
    call die('Not an electron-phonon calculation, but calculating dVXC!')
  endif
  if (vxc_type .eq. 3) then
    if ((.not.present(wfnk_phonq)) .or. (.not.present(wfnk_phonq_off))) then
      call die('Calculation of dVXC mtxel needs wfnk_phonq and wfnk_phonq_off.')
    endif
  endif
! Nullify output
  alda=0.0d0
! Using FFT to compute matrix elements
! Allocate temporary wavefunction array
  if (sig%noffdiag.gt.0) then
    if (mod(peinf%inode,peinf%npes/peinf%npools).eq.0) then
      allocate(wfnkoff%zk (2*wfnk%nkpt,kp%nspinor))
      if (vxc_type .eq. 3) then
        allocate(wfnk_phonq_off%zk (2*wfnk_phonq%nkpt,kp%nspinor)) ! nkpt is nkg, 2 is for left and right two wfns
      endif
    endif
  endif
! Get FFT box sizes and scale factor
  call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)
! Get FFT box sizes and scale factor
! Allocate FFT boxes
  allocate(fftbox1 (Nfft(1),Nfft(2),Nfft(3),sig%nspin))
  allocate(fftbox2 (Nfft(1),Nfft(2),Nfft(3),kp%nspinor))
  allocate(fftbox3 (Nfft(1),Nfft(2),Nfft(3),kp%nspinor))
  allocate(vxctemp (gvec%ng,sig%nspin))
! Put the vxc data into fftbox 1 and FFT it to real space
  allocate(gvecindex (gvec%ng))
  do ig=1,gvec%ng
    gvecindex(ig)=ig
  enddo
  do ispin=1,sig%nspin
    if (vxc_type .eq. 1) then
      vxctemp(:,ispin) = sig%vxc(:,sig%spin_index(ispin))
    else if (vxc_type .eq. 2) then
      vxctemp(:,ispin) = sig%vxc2(:,sig%spin_index(ispin))
    else if (vxc_type .eq. 3) then
      vxctemp(:,ispin) = sig%dvxc(:,sig%spin_index(ispin))
    else
      call die('vxc_type can only be 1, 2, 3, for VXC, VXC2, dVXC, respectively')
    endif
    call put_into_fftbox(gvec%ng,vxctemp(:,ispin),gvec%components,gvecindex,fftbox1(:,:,:,ispin),Nfft)
    call do_FFT(fftbox1(:,:,:,ispin),Nfft,1)
    ! ZL: check dVXC
! if (vxc_type .eq. 3) then
! if (peinf%inode.eq.0) write(*,*) "check VXC after FFT to real space", vxc_type
! if (peinf%inode.eq.0) write(*,*) fftbox1(:,:,:,ispin)
! endif
! Loop over the bands for which we need the matrix elements
! For each one, put the band into fftbox2, FFT to real space,
! and then integrate in real space vxc(r)*|psi(r)|^2
! Store result into alda(in,ispin).
    do in=1,peinf%ndiag_max
      temp=0.0d0
      write(tmpstr,'("Computing <n|Vxc|n> for n=",i4)') &
        sig%diag(peinf%index_diag(in))
      call logit(tmpstr)
        if (mod(peinf%inode,peinf%npes/peinf%npools).eq.0) then
          do jsp=1,kp%nspinor
            ij = ispin*jsp ! for shortening lengthy call to put_into_fftbox
            call put_into_fftbox(wfnk%nkpt,wfnk%zk((in-1)*wfnk%nkpt+1:,ij),gvec%components,wfnk%isrtk,fftbox2(:,:,:,jsp),Nfft)
            call do_FFT(fftbox2(:,:,:,jsp),Nfft,1)
            if (vxc_type .eq. 3) then
              call put_into_fftbox(wfnk_phonq%nkpt,wfnk_phonq%zk((in-1)*wfnk_phonq%nkpt+1:,ij), &
                                   gvec%components,wfnk_phonq%isrtk,fftbox3(:,:,:,jsp),Nfft)
              call do_FFT(fftbox3(:,:,:,jsp),Nfft,1)
              call conjg_fftbox(fftbox3(:,:,:,jsp),Nfft)
            endif
          enddo
          do iz=1,Nfft(3)
            do iy=1,Nfft(2)
              do ix=1,Nfft(1)
                do jsp=1,kp%nspinor
                  if (vxc_type .eq. 3) then
                    temp = temp + fftbox1(ix,iy,iz,ispin)*fftbox2(ix,iy,iz,jsp)*fftbox3(ix,iy,iz,jsp)
                  else
                    temp = temp + fftbox1(ix,iy,iz,ispin) * abs(fftbox2(ix,iy,iz,jsp))**2
                  endif
                enddo
              enddo
            enddo
          enddo
        endif
      if (peinf%flag_diag(in)) then
        alda(peinf%index_diag(in),ispin)=temp*scale
      endif
    enddo ! in
    do ioff=1,peinf%noffdiag_max
      write(tmpstr,'("Computing offdiag <n|Vxc|m> for n,m=",2i4)') &
        sig%off1(peinf%index_offdiag(ioff)), &
        sig%off2(peinf%index_offdiag(ioff))
      call logit(tmpstr)
      temp=0.0d0
! (gsm) begin gathering wavefunctions over pools
! $$$ inefficient communication, this should be rewritten $$$
      do jj=1,peinf%npools
        dest=(jj-1)*(peinf%npes/peinf%npools)
        do ii=1,2
          iso=sig%offmap(peinf%index_offdiag(ioff),ii)
          jso=(iso-1)/peinf%npools+1
          source=mod(iso-1,peinf%npools)*(peinf%npes/peinf%npools)
! ZL: for wfnk and wfnkoff
          if (peinf%inode.eq.source.and.peinf%inode.eq.dest) then
            do jsp=1,kp%nspinor
              do ik=1,wfnk%nkpt
                wfnkoff%zk((ii-1)*wfnk%nkpt+ik,jsp)= &
                  wfnk%zk((jso-1)*wfnk%nkpt+ik,ispin*jsp)
              enddo
            enddo
          else
            do jsp=1,kp%nspinor
              do ik=1,wfnk%nkpt
                wfnkoff%zk((ii-1)*wfnk%nkpt+ik,jsp)= &
                  wfnk%zk((jso-1)*wfnk%nkpt+ik,ispin*jsp)
              enddo
            enddo
          endif
! ZL: for wfnk_phonq and wfnk_phonq_off
          if (vxc_type .eq. 3) then
            if (peinf%inode.eq.source.and.peinf%inode.eq.dest) then
              do jsp=1,kp%nspinor
                do ik=1,wfnk_phonq%nkpt
                  wfnk_phonq_off%zk((ii-1)*wfnk_phonq%nkpt+ik,jsp)= &
                    wfnk_phonq%zk((jso-1)*wfnk_phonq%nkpt+ik,ispin*jsp)
                enddo
              enddo
            else
              do jsp=1,kp%nspinor
                do ik=1,wfnk_phonq%nkpt
                  wfnk_phonq_off%zk((ii-1)*wfnk_phonq%nkpt+ik,jsp)= &
                    wfnk_phonq%zk((jso-1)*wfnk_phonq%nkpt+ik,ispin*jsp)
                enddo
              enddo
            endif
          endif ! (vxc_type .eq. 3)
! ZL: done for wfnk, wfnkoff, wfnk_phonq, wfnk_phonq_off
        enddo
      enddo
! (gsm) end gathering wavefunctions over pools
      if (mod(peinf%inode,peinf%npes/peinf%npools).eq.0) then
        do jsp=1,kp%nspinor
          ! ZL: the index here are important, as they represent n, m
          if (vxc_type .eq. 3) then
            call put_into_fftbox(wfnk_phonq%nkpt,wfnk_phonq_off%zk(1:,jsp), &
                                 gvec%components,wfnk_phonq%isrtk,fftbox2(:,:,:,jsp),Nfft)
          else
            call put_into_fftbox(wfnk%nkpt,wfnkoff%zk(1:,jsp),gvec%components,wfnk%isrtk,fftbox2(:,:,:,jsp),Nfft)
          endif
          call do_FFT(fftbox2(:,:,:,jsp),Nfft,1)
          call conjg_fftbox(fftbox2(:,:,:,jsp),Nfft)
          call put_into_fftbox(wfnk%nkpt,wfnkoff%zk(wfnk%nkpt+1:,jsp),gvec%components,wfnk%isrtk,fftbox3(:,:,:,jsp),Nfft)
          call do_FFT(fftbox3(:,:,:,jsp),Nfft,1)
        enddo
        do jsp=1,kp%nspinor
          do iz=1,Nfft(3)
            do iy=1,Nfft(2)
              do ix=1,Nfft(1)
                temp = temp + fftbox1(ix,iy,iz,ispin)*fftbox2(ix,iy,iz,jsp)*fftbox3(ix,iy,iz,jsp)
              enddo
            enddo
          enddo
        enddo
      endif
      if (peinf%flag_offdiag(ioff)) then
        alda(peinf%index_offdiag(ioff)+sig%ndiag,ispin)=temp*scale*ryd
      endif
    enddo ! ioff
  enddo ! ispin
! Deallocate temporary wavefunction array
  if (sig%noffdiag.gt.0) then
    if (mod(peinf%inode,peinf%npes/peinf%npools).eq.0) then
      if(associated(wfnkoff%zk))then;deallocate(wfnkoff%zk);nullify(wfnkoff%zk);endif
      if (vxc_type .eq. 3) then
        if(associated(wfnk_phonq_off%zk))then;deallocate(wfnk_phonq_off%zk);nullify(wfnk_phonq_off%zk);endif
      endif
    end if
  endif
! Deallocate FFT boxes
  if(allocated(gvecindex))then;deallocate(gvecindex);endif
  if(allocated(fftbox1))then;deallocate(fftbox1);endif
  if(allocated(fftbox2))then;deallocate(fftbox2);endif
  if(allocated(fftbox3))then;deallocate(fftbox3);endif
! If MPI add up all the work done in parallel
 
  return
end subroutine mtxel_vxc
end module mtxel_vxc_m
