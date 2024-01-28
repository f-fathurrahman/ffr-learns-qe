!=================================================================================
!
! Routines:
!
! (1) g_sum() Originally By JRD Last Modified 2/11/2012 (FHJ)
!
! Performs the sum over G for each (v,c,v',c') set.
!
! This routine scales as N^5. This routine represents the worst
! scaling with number of atoms in Kernel code.
!
! Bow down before it.
!
!=================================================================================
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
module g_sum_m
  use global_m
  use algos_kernel_m
  use blas_m
  implicit none
  public :: g_sum_TDA, g_sum_extended
  private
contains
  !> FHJ: Performs the summation over Gp to get the direct kernel term.
  !! This routine is specialized for TDA kernels.
  subroutine g_sum_TDA(xct,invband,incband,temph,tempb,mccp, &
    bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik,&
    tempw,bsedwing)
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: invband,incband,leading_dim
    integer, intent(in) :: ivp_in,icp_in,ikp,iv_in,ic_in,ik
    real(DP), intent(in), target :: tempb(:,:,:,:), temph(:,:,:), mccp(:,:,:,:)
    real(DP), intent(inout), target :: bsedhead(:,:,:), bsedbody(:,:,:)
    real(DP), intent(in), optional, target :: tempw(:,:,:,:)
    real(DP), intent(inout), optional, target :: bsedwing(:,:,:)
   
    if (present(bsedwing) .and. .not. present(tempw)) then
      call die("Internal error in g_sum_TDA: bsedwing present but tempw not present", &
               only_root_writes = .true.)
    end if
    select case (g_sum_algo)
    case (CPU_ALGO)
      if (present(bsedwing)) then
        call g_sum_TDA_cpu(xct, invband, incband, temph, tempb, mccp, &
                           bsedhead, bsedbody, leading_dim, ivp_in, icp_in, &
                           ikp, iv_in, ic_in, ik, tempw, bsedwing)
      else
        call g_sum_TDA_cpu(xct, invband, incband, temph, tempb, mccp, &
                           bsedhead, bsedbody, leading_dim, ivp_in, icp_in, &
                           ikp, iv_in, ic_in, ik)
      end if
    case (OPENACC_ALGO)
      if (present(bsedwing)) then
        call g_sum_TDA_openacc(xct, invband, incband, temph, tempb, mccp, &
                           bsedhead, bsedbody, leading_dim, ivp_in, icp_in, &
                           ikp, iv_in, ic_in, ik, tempw, bsedwing)
      else
        call g_sum_TDA_openacc(xct, invband, incband, temph, tempb, mccp, &
                           bsedhead, bsedbody, leading_dim, ivp_in, icp_in, &
                           ikp, iv_in, ic_in, ik)
      end if
    case default
      call die("Invald algorithm for g_sum_TDA", &
               only_root_writes = .true.)
    end select
   
  end subroutine g_sum_TDA
  subroutine g_sum_TDA_cpu(xct,invband,incband,temph,tempb,mccp, &
    bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik,&
    tempw,bsedwing)
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: invband,incband,leading_dim
    integer, intent(in) :: ivp_in,icp_in,ikp,iv_in,ic_in,ik
    real(DP), intent(in) :: tempb(:,:,:,:), temph(:,:,:), mccp(:,:,:,:)
    real(DP), intent(inout) :: bsedhead(:,:,:), bsedbody(:,:,:)
    real(DP), intent(in), optional :: tempw(:,:,:,:)
    real(DP), intent(inout), optional :: bsedwing(:,:,:)
    integer :: isv, isc, iv, ivp, ic, icp, iit, iitend, iitbeG
   
    do isv=1,xct%nspin
      do isc=1,xct%nspin
        if (xct%icpar .eq. 0) then
          iit = peinf%wown(1,1,ikp,1,1,ik) - 1
        else if (xct%ivpar .eq. 0) then
          iit = peinf%wown(1,icp_in,ikp,1,ic_in,ik) - 1
        else
          iit = peinf%wown(ivp_in,icp_in,ikp,iv_in,ic_in,ik) - 1
        endif
        iitbeg = iit + 1
        iitend = iit + invband*invband*incband*incband
        if (incband .gt. 1) then
          call dgemm('t','n',invband*invband,incband*incband,xct%ng, &
            1.0d0,tempb(:,:,:,isv),xct%ng,mccp(:,:,:,isc),xct%ng, &
            1.0d0,bsedbody(iitbeg:iitend,isc,isv),invband*invband)
          if (present(bsedwing)) &
            call dgemm('t','n',invband*invband,incband*incband,xct%ng, &
            1.0d0,tempw(:,:,:,isv),xct%ng,mccp(:,:,:,isc),xct%ng, &
            1.0d0,bsedwing(iitbeg:iitend,isc,isv),invband*invband)
        else
          call dgemv('t',xct%ng,invband*invband,1.0d0,tempb(:,:,:,isv),xct%ng,mccp(:,:,:,isc), &
            1,0.0d0,bsedbody(iitbeg:iitend,isc,isv),1)
          if (present(bsedwing)) &
            call dgemv('t',xct%ng,invband*invband,1.0d0,tempw(:,:,:,isv),xct%ng,mccp(:,:,:,isc), &
            1,0.0d0,bsedwing(iitbeg:iitend,isc,isv),1)
        endif
        if (iit + (incband*invband)**2 > leading_dim) then
          write(0,*) peinf%inode, ik, ikp
          write(0,*) peinf%wown(1,1,ikp,1,1,ik), iit
          write(0,*) peinf%nckpe, peinf%myown, leading_dim
          call die("Internal error in g_sum with array dimensions.")
        endif
        ! Head
        do icp = 1, incband
          do ic = 1, incband
            do ivp = 1, invband
              do iv = 1, invband
                iit = iit + 1
                bsedhead(iit,isc,isv) = bsedhead(iit,isc,isv) + &
                  mccp(1,ic,icp,isc) * temph(iv,ivp,isv)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
   
  end subroutine g_sum_TDA_cpu
  subroutine g_sum_TDA_openacc(xct,invband,incband,temph,tempb,mccp, &
    bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik,&
    tempw,bsedwing)
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: invband,incband,leading_dim
    integer, intent(in) :: ivp_in,icp_in,ikp,iv_in,ic_in,ik
    real(DP), intent(in), target :: tempb(:,:,:,:), temph(:,:,:), mccp(:,:,:,:)
    real(DP), intent(inout), target :: bsedhead(:,:,:), bsedbody(:,:,:)
    real(DP), intent(in), optional, target :: tempw(:,:,:,:)
    real(DP), intent(inout), optional, target :: bsedwing(:,:,:)
   
    call die("OpenACC version of g_sum requested, but OpenACC not compiled "//&
             "into this executable", only_root_writes = .true.)
   
  end subroutine g_sum_TDA_openacc
  !> FHJ: Performs the summation over Gp to get the direct kernel term.
  !! This routine is generalized for extended kernels calculations.
  subroutine g_sum_extended(xct,ofs2,ofs2p,n2,n2p,temph,tempb,m22p, &
    bsedhead,bsedbody,leading_dim,ivp_in,icp_in,ikp,iv_in,ic_in,ik, &
    tempw, bsedwing)
    type (xctinfo), intent(in) :: xct
    integer, intent(in) :: ofs2, ofs2p
    integer, intent(in) :: n2, n2p
    integer, intent(in) :: leading_dim
    integer, intent(in) :: ivp_in,icp_in,ikp,iv_in,ic_in,ik
    real(DP), intent(in) :: tempb(:,:,:,:), temph(:,:,:), m22p(:,:,:,:)
    real(DP), intent(inout) :: bsedhead(:,:,:), bsedbody(:,:,:)
    real(DP), intent(in), optional :: tempw(:,:,:,:)
    real(DP), intent(inout), optional :: bsedwing(:,:,:)
    integer :: isv, isc, i1, i1p, i2, i2p, it, itbeg, it_buf
    integer :: n_left, n_right, buf_sz
    real(DP), pointer :: buf_b(:), buf_w(:)
   
    ! FHJ: n_left/n_right is the total number of val/cond states that this PE
    ! deals with. The number will be nv+nc for extended kernel calculations.
    ! However, whenever we deal with m22p, we should actually use n2/n2p +
    ! offsets, b/c the matrix m22p is decomposed in blocks.
    n_left = xct%n1b_co
    if (xct%ivpar==1) n_left = 1
    n_right = xct%n2b_co
    if (xct%icpar==1) n_right = 1
    buf_sz = n_left**2 * n2 * n2p
    allocate(buf_b (buf_sz))
    if (present(bsedwing)) then
      allocate(buf_w (buf_sz))
    endif
    do isv=1,xct%nspin
      do isc=1,xct%nspin
        if (xct%icpar .eq. 0) then
          itbeg = peinf%wown(1,1,ikp,1,1,ik)
        else if (xct%ivpar .eq. 0) then
          itbeg = peinf%wown(1,icp_in,ikp,1,ic_in,ik)
        else
          itbeg = peinf%wown(ivp_in,icp_in,ikp,iv_in,ic_in,ik)
        endif
        if (n2>1 .or. n2p>1) then
          call dgemm('t','n',n_left**2,n2*n2p,xct%ng, &
            1.0d0,tempb(:,:,:,isv),xct%ng,m22p(:,:,:,isc),xct%ng, &
            0.0d0,buf_b(:),n_left**2)
          if (present(bsedwing)) &
            call dgemm('t','n',n_left**2,n2*n2p,xct%ng, &
            1.0d0,tempw(:,:,:,isv),xct%ng,m22p(:,:,:,isc),xct%ng, &
            0.0d0,buf_w(:),n_left**2)
        else
          call dgemv('t',xct%ng,n_left**2,1.0d0,tempb(:,:,:,isv),xct%ng,m22p(:,:,:,isc), &
            1,0.0d0,buf_b(:),1)
          if (present(bsedwing)) &
            call dgemv('t',xct%ng,n_left**2,1.0d0,tempw(:,:,:,isv),xct%ng,m22p(:,:,:,isc), &
            1,0.0d0,buf_w(:),1)
        endif
        if (itbeg-1 + (n_left*n_right)**2 > leading_dim) then
          write(0,*) peinf%inode, ik, ikp
          write(0,*) peinf%wown(1,1,ikp,1,1,ik), itbeg, n_left, n_right
          write(0,*) peinf%nckpe, peinf%myown, leading_dim
          call die("Internal error in g_sum with array dimensions.")
        endif
        do i2p = 1, n2p
          do i2 = 1, n2
            ! FHJ: max(it_buf) = n_left**2 * n2*n2p
            it_buf = 1 + n_left*n_left*((i2-1) + n2*(i2p-1))
            ! FHJ: max(it) = n_left**2 * n_right**2
            it = itbeg + n_left*n_left*((ofs2+i2-1) + n_right*(ofs2p+i2p-1))
            do i1p = 1, n_left
              do i1 = 1, n_left
                bsedhead(it,isc,isv) = temph(i1,i1p,isv) * m22p(1,i2,i2p,isc)
                if (present(bsedwing)) bsedwing(it,isc,isv) = buf_w(it_buf)
                bsedbody(it,isc,isv) = buf_b(it_buf)
                it = it + 1
                it_buf = it_buf + 1
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    if(associated(buf_b))then;deallocate(buf_b);nullify(buf_b);endif
    if (present(bsedwing)) then
      if(associated(buf_w))then;deallocate(buf_w);nullify(buf_w);endif
    endif
   
  end subroutine g_sum_extended
end module g_sum_m
