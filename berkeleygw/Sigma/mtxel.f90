!============================================================================
!
! Routines:
!
! (1) mtxel() Originally By ? Last Modified 7/8/2008 (JRD)
!
! Subroutine computes required matrix elements
! of the form <nn,k|exp{i(q+G).r}|n1,k-q> = M_{nn,n1}(k-q,q,G)
!
! FHJ: Note that eqns. (20-24) of the BerkleyGW arxiv paper use the
! quantity [M_{n``,n}(k,-q,-G)]^*, which is the same as what we compute
! here, M_{nn,n1}(k-q,q,G) = aqs(G,n1). The indices in the arxiv paper and
! in the code are related, respectively, by:
! - n and n` <-> nn
! - n`` <-> n1
!
! input nn band index for "outer" band
! input ncoul number of matrix elements required
! input isrtrq index array for g-vectors in
! <nk|exp(i(q+g).r)|n1k-q>
! output aqs matrix elements required
!
!============================================================================

module mtxel_m
  use algos_sigma_m
  use fftw_m
  use global_m
  use misc_m
  use acc_mtxel_kernels_m, only: acc_mtxel_sig, &
                                 acc_zero_box, acc_put_into_fftbox, acc_box_multiply, acc_get_from_fftbox, &
                                 acc_run_fft
  implicit none
  private
  public :: &
    mtxel
contains
subroutine mtxel(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: nn
  type (gspace), intent(in) :: gvec
  type (wfnkqstates), intent(in) :: wfnkq
  type (wfnkstates), intent(in) :: wfnk
  integer, intent(in) :: ncoul
  integer, intent(in) :: isrtrq(gvec%ng)
  real(DP), intent(out) :: aqs(ncoul,peinf%ntband_max)
  integer, intent(in) :: ispin
  type (kpoints), intent(inout) :: kp
 
  select case (mtxel_algo)
  case (CPU_ALGO)
    call mtxel_cpu(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  case (OPENACC_ALGO)
    call mtxel_openacc(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  case (OMP_TARGET_ALGO)
    call mtxel_omp_target(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  case default
    call die("Invald algorithm for mtxel", only_root_writes = .true.)
  end select
 
  return
end subroutine mtxel
subroutine mtxel_cpu(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: nn
  type (gspace), intent(in), target :: gvec
  type (wfnkqstates), intent(in), target :: wfnkq
  type (wfnkstates), intent(in), target :: wfnk
  integer, intent(in), target :: isrtrq(gvec%ng)
  type (kpoints), intent(inout) :: kp
  integer, intent(in) :: ncoul
  real(DP), intent(out) :: aqs(ncoul,peinf%ntband_max)
  integer, intent(in) :: ispin
  integer :: n1,jsp
! We use FFT to compute <u_nn,k|e^(iG.r)|u_n1,k-q> elements where
! u_nk is the periodic part of the wave function.
! The calculation is done in real space, and integration over
! the grid is replaced by the sum over the grid points p:
!
! <u_nn,k|e^(iG.r)|u_n1,k-q> =
! Volume/Np * sum_p { conj(u_nn,k(p))*e^(iG.p)*u_n1k-q(p) }
!
! Since u_nk(p) = Volume^-0.5 * sum_G { cnk(G)*e^(iG.p) },
! and FFT is defined as FFT(cnk,+,p) = sum_G { cnk(G)*e^{+iG.p} },
! we must compute
!
! <u_nn,k|e^(iG.r)|u_n1,k-q>
! = 1/Np * sum_p { conj(FFT(c_nn k,+,p))*e^(iG.p)*FFT(c_n1 k-q,+,p) }
! = 1/Np * FFT(conj(FFT(c_nn k,+,:)).*FFT(c_n1 k-q,+,:),+,G)
!
! where .* is a point by point multiplication on the grid
  complex(DPC), dimension(:,:,:), allocatable :: fftbox1,fftbox2
  integer, dimension(3) :: Nfft
  real(DP) :: scale
  real(DP), dimension(:), allocatable :: tmparray
! Compute size of FFT box we need and scale factor
 
  call setup_FFT_sizes(gvec%FFTgrid,Nfft,scale)
! Allocate FFT boxes
  allocate(fftbox1 (Nfft(1),Nfft(2),Nfft(3)))
  allocate(fftbox2 (Nfft(1),Nfft(2),Nfft(3)))
! Put the data for band nn into FFT box 1 and do the FFT,zk(:,1)
  allocate(tmparray (ncoul))
  do jsp = ispin,ispin*kp%nspinor
    call put_into_fftbox(wfnk%nkpt,wfnk%zk((nn-1)*wfnk%nkpt+1:,jsp),gvec%components,wfnk%isrtk,fftbox1,Nfft)
    call do_FFT(fftbox1,Nfft,1)
    ! We need the complex conjugate of u_{nn,k)(r) for the cross correlation
    call conjg_fftbox(fftbox1,Nfft)
! Now we loop over the n1 states and get the matrix elements:
! Get n1 wave function and put it into box 2,
! do FFT, get u_{n1k-q}(r)
! multiply by box1 contents, get
! do FFT again,
! and extract the resulting matrix elements
! Now we loop over the n1 states and get the matrix elements:
! 1. Get conduction wave function and put it into box 2,
! 2. do FFT, get u_{n1,k-q}(r)
! 3. multiply by box1 contents, get F(r) = [u_{nn,k)(r)]^* u_{n1,k-q}(r)
! 4. do FFT again, and extract the resulting matrix elements,
! <nn,k|exp{i(q+G).r}|n1,k-q>
    do n1=1,peinf%ntband_node
      call put_into_fftbox(wfnkq%nkpt,wfnkq%zkq((n1-1)*wfnkq%nkpt+1:,jsp),gvec%components,wfnkq%isrtkq,fftbox2,Nfft)
      call do_FFT(fftbox2,Nfft,1)
      call multiply_fftboxes(fftbox1,fftbox2,Nfft)
      call do_FFT(fftbox2,Nfft,1)
      call get_from_fftbox(ncoul,tmparray,gvec%components,isrtrq,fftbox2,Nfft,scale)
      if (kp%nspinor.eq.1 .or. jsp.eq. 1) then
        aqs(:,n1) = tmparray(:)
      else
        aqs(:,n1) = aqs(:,n1) + tmparray(:)
      endif
    enddo
  enddo
  if(allocated(tmparray))then;deallocate(tmparray);endif
! We are done, so deallocate FFT boxes
  if(allocated(fftbox1))then;deallocate(fftbox1);endif
  if(allocated(fftbox2))then;deallocate(fftbox2);endif
 
  return
end subroutine mtxel_cpu
subroutine mtxel_openacc(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: nn
  type (gspace), intent(in):: gvec
  type (wfnkqstates), intent(in) :: wfnkq
  type (wfnkstates), intent(in) :: wfnk
  integer, intent(in) :: isrtrq(gvec%ng)
  type (kpoints), intent(inout) :: kp
  integer, intent(in) :: ncoul
  real(DP), intent(out) :: aqs(ncoul,peinf%ntband_max)
  integer, intent(in) :: ispin
 
  call die("OpenACC version of mtxel requested, but OpenACC not compiled&
           & into this executable", only_root_writes = .true.)
 
  return
end subroutine mtxel_openacc
subroutine mtxel_omp_target(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: nn
  type (gspace), intent(in):: gvec
  type (wfnkqstates), intent(in) :: wfnkq
  type (wfnkstates), intent(in) :: wfnk
  integer, intent(in) :: isrtrq(gvec%ng)
  type (kpoints), intent(inout) :: kp
  integer, intent(in) :: ncoul
  real(DP), intent(out) :: aqs(ncoul,peinf%ntband_max)
  integer, intent(in) :: ispin
 
  call die("OpenMP Target version of mtxel requested, but OpenMP Target&
           & not compiled into this executable", only_root_writes = .true.)
 
  return
end subroutine mtxel_omp_target
subroutine mtxel_openacc_2(nn,gvec,wfnkq,wfnk,ncoul,isrtrq,aqs,ispin,kp)
  integer, intent(in) :: nn
  type (gspace), intent(in):: gvec
  type (wfnkqstates), intent(in) :: wfnkq
  type (wfnkstates), intent(in) :: wfnk
  integer, intent(in) :: isrtrq(gvec%ng)
  type (kpoints), intent(inout) :: kp
  integer, intent(in) :: ncoul
  real(DP), intent(out) :: aqs(ncoul,peinf%ntband_max)
  integer, intent(in) :: ispin
 
  call die("OpenACC version of mtxel (2) requested, but OpenACC not compiled&
          & into this executable", only_root_writes = .true.)
 
  return
end subroutine mtxel_openacc_2
end module mtxel_m
