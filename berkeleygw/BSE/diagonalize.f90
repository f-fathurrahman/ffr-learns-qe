!==============================================================================
!
! Routines:
!
! (1) diagonalize() Originally AC Last Modified 01/2015 (FHJ)
!
! Generalised parallel Hermitian eigenvalue solver,
! Returns eigenvectors and eigenvalues of the BSE Hamiltonian hbse_a.
!
! neig = number of eigenvectors
! nmat = lda of hbse_a and sis (=nv in code)
!
! interface routine for the pzheevx/pdsyevx, generalised parallel eigenvalue
! solver. Starts with the distributed S and H matrices. Redistributes
! them to blacs layout, then calls pzheevx/pdsyevx and
! then collects all the eigenvectors onto
! all pes. For more details of the scalapack routines and data layouts
! see http://www.netlib.org/scalapack/scalapack_home.html
!
! based on pssyevx_inter/pcheevx_inter, written by Andrew Canning
! NERSC/LBL 1998
!
! distributed solver (block cyclic blacs layout):
!
! nbl = blocksize
! nprow = processor grid row
! npcol = processor grid column
!
! double-precision eigenvalue solvers, pzheevx/pdsyevx
!
!===================================================================================

module diagonalize_m
  use global_m
  use scalapack_m
  implicit none
  private
  public :: &
    diagonalize
contains
subroutine diagonalize(xct, neig, nmat, hbse_a, evals, evecs_r, hbse_b, evecs_l)
  type(xctinfo), intent(in) :: xct
  integer, intent(in) :: neig
  integer, intent(in) :: nmat
  real(DP), intent(inout) :: hbse_a(nmat,peinf%nblocks*peinf%block_sz)
  real(DP), intent(out) :: evals(neig)
  real(DP), intent(out) :: evecs_r(:,:) !<({nmat,2*nmat},pblock)
  real(DP), intent(inout), optional :: hbse_b(nmat,peinf%nblocks*peinf%block_sz)
  real(DP), intent(out), optional :: evecs_l(:,:) !<({nmat,2*nmat},pblock)
  integer :: ii, info
  character :: range
  character*100 :: tmpstr
 
  if (.not.xct%tda.and.(.not.present(hbse_b).or..not.present(evecs_l))) then
    call die('Internal parameter error in diagonalize.', only_root_writes=.true.)
  endif
! If neig > nmat, only the first neig eigenvectors/eigenvalues will
! be computed. Otherwise, calculate all of them.
! Set range
  range = 'A'
  if (xct%tda) then
    if (neig==nmat) then
      range='A'
    else if (neig<nmat) then
      range='I'
    else
      write(tmpstr,'(a,i0,a,i0)') 'diagonalize: ', neig, ' eigenvalues requested > matrix size ', nmat
      call die(tmpstr)
    endif
  endif
  ! FHJ: Full BSE
  if (.not.xct%tda) then
    ! FHJ: Even if we have only 1 processor, we use the ScaLAPACK routines from
    ! SSEIG because it`s faster than our pure LAPACK algorithm.
    call serial_full()
    !call serial_full_gvx()
   
    return
  endif
  ! FHJ: Serial TDA
  if (peinf%npes==1) then
    call serial_tda()
   
    return
  endif
  ! FHJ: Parallel TDA below
! End of USESCALAPACK region
  ! FHJ: If we got here, we had a weird combination of npes>0 but no ScaLAPACK
  call die('Running with more than one processor without ScaLAPACK!')
 
contains
  !============================================================================
  !> Serial code to diagonalize the TDA Hamiltonian.
  !============================================================================
  subroutine serial_tda()
    real(DP) :: evals_t(nmat)
    real(DP) :: abstol
    integer, allocatable :: iwork(:), ifail(:)
    real(DP), allocatable :: work(:)
    integer :: iloop, jloop
    integer :: nfound, ilow, iup, info
    real(DP) :: ellow, elup
   
    if (nmat.ne.peinf%nblocks*peinf%block_sz) then
      call die('Hamiltonian matrix does not seem to be square!')
    endif
    if (peinf%inode==0) write(6,*)
    call calc_broken_herm(hbse_a, nmat, 'H')
    if (peinf%inode==0) write(6,*)
    abstol=0.0
    allocate(ifail (nmat))
    allocate(work (10*nmat))
    allocate(iwork (5*nmat))
    ilow = 1
    iup = neig
    call dsyevx('V',range,'U',nmat,hbse_a,nmat,ellow,elup,ilow,iup, &
      abstol,nfound,evals_t,evecs_r,nmat,work,10*nmat,iwork,ifail,info)
    if(nfound .lt. neig) then
      write(tmpstr,'(a, i10, a, i10, a)') 'Diagonalization with dsyevx failed: only ', &
        nfound, ' of ', neig, ' eigenvectors found.'
      call die(tmpstr)
    endif
    if(info.lt.0) then
      write(tmpstr,*) "Problem in input parameters for zheevx/dsyevx: info = ",info
      call die(tmpstr)
    endif
    if(info.gt.0) then
      write(0,*) "Convergence problems in zheevx/dsyevx: info = ",info
      write(tmpstr,*) 'The following eigenvector failed to converge: ifail = ',ifail
      call die(tmpstr)
    endif
    evals(1:neig) = evals_t(1:neig)
    if(allocated(iwork))then;deallocate(iwork);endif
    if(allocated(work))then;deallocate(work);endif
    if(allocated(ifail))then;deallocate(ifail);endif
   
  end subroutine serial_tda
  !============================================================================
  !> Serial code to diagonalize the full BSE Hamiltonian (Originally by FHJ).
  !============================================================================
  subroutine serial_full()
    real(DP), allocatable :: work(:)
    real(DP) :: evals_r(2*nmat), evals_i(2*nmat)
    real(DP) :: hbse(2*nmat,2*nmat), tmp_norm
    integer :: lwork, info, jj
   
    if (nmat.ne.peinf%nblocks*peinf%block_sz) then
      call die('Hamiltonian matrix does not seem to be square!')
    endif
    if (peinf%inode==0) write(6,*)
    call calc_broken_herm(hbse_a, nmat, 'A')
    call calc_broken_herm(hbse_b, nmat, 'B')
    if (peinf%inode==0) write(6,*)
    ! FHJ: symmetrize hbse_a and hbse_b using the lower triangular portions.
    ! We need to to this to obtain the same result as SSEIG solver.
    do ii=1,nmat
      hbse_a(ii, ii+1:nmat) = (hbse_a(ii+1:nmat, ii))
      hbse_b(ii, ii+1:nmat) = hbse_b(ii+1:nmat, ii)
    enddo
    hbse(1:nmat, 1:nmat) = hbse_a(1:nmat, 1:nmat)
    hbse(nmat+1:2*nmat, nmat+1:2*nmat) = -(hbse_a(1:nmat, 1:nmat))
    hbse(1:nmat, nmat+1:2*nmat) = hbse_b(1:nmat, 1:nmat)
    hbse(nmat+1:2*nmat, 1:nmat) = -(hbse_b(1:nmat, 1:nmat))
    allocate(work (10))
    write(6,'(1x,a,i0)') 'Beginning LAPACK diagonalization. Size: ', 2*nmat
    call dgeev('V', 'V', 2*nmat, hbse, 2*nmat, &
      evals_r, evals_i, &
      evecs_l, 2*nmat, evecs_r, 2*nmat, work, -1, &
      info)
    if(info/=0) then
      write(tmpstr,*) "problem in xgeev/query mode: got info = ", info
      call die(tmpstr)
    endif
    lwork = max(1, int(work(1)))
    if(allocated(work))then;deallocate(work);endif
    allocate(work (lwork))
    call dgeev('V', 'V', 2*nmat, hbse, 2*nmat, &
      evals_r, evals_i, &
      evecs_l, 2*nmat, evecs_r, 2*nmat, work, lwork, &
      info)
    if(info/=0) then
      write(tmpstr,*) "problem in xgeev: got info = ", info
      call die(tmpstr)
    endif
    write(6,*) 'Done with LAPACK diagonalization.'
    evals(1:neig) = evals_r(1:neig)
    ! FHJ: No, the left/right eigenvectors that come out of LAPACK are *not*
    ! correctly normalized against each other!
    do jj=1,2*nmat
      tmp_norm = sum(evecs_l(1:2*nmat,jj)*evecs_r(1:2*nmat,jj))
      evecs_l(1:2*nmat,jj) = evecs_l(1:2*nmat,jj) / tmp_norm
    enddo
   
  end subroutine serial_full
  !============================================================================
  !> Serial code to diagonalize the full BSE Hamiltonian (Orig. by FHJ).
  !! This version calls a generalized eigensolver from LAPACK.
  !============================================================================
  subroutine serial_full_gvx()
    real(DP) :: hbse(2*nmat,2*nmat)
    integer :: lwork, liwork, lrwork, nfound, ii
    integer :: ifail(2*nmat), iwork(10*nmat)
    real(DP) :: b_mat(2*nmat,2*nmat), s_mat(2*nmat,2*nmat)
    real(DP), allocatable :: work(:)
    real(DP) :: abstol, evals_t(2*nmat)
   
    ! FHJ: zero coupling block if the user requests. Useful to test the solver.
    if (xct%zero_coupling_block) hbse_b = 0.0d0
    if (nmat.ne.peinf%nblocks*peinf%block_sz) then
      call die('Hamiltonian matrix does not seem to be square!')
    endif
    ! FHJ: symmetrize hbse_a and hbse_b using the lower triangular portions.
    ! We need to to this to obtain the same result as SSEIG solver.
    do ii=1,nmat
      hbse_a(ii, ii+1:nmat) = (hbse_a(ii+1:nmat, ii))
      hbse_b(ii, ii+1:nmat) = hbse_b(ii+1:nmat, ii)
    enddo
    hbse(1:nmat, 1:nmat) = hbse_a(1:nmat, 1:nmat)
    hbse(nmat+1:2*nmat, nmat+1:2*nmat) = (hbse_a(1:nmat, 1:nmat))
    hbse(1:nmat, nmat+1:2*nmat) = hbse_b(1:nmat, 1:nmat)
    hbse(nmat+1:2*nmat, 1:nmat) = (hbse_b(1:nmat, 1:nmat))
    ! Create matrix B
    ! b_mat = [ 1 0 ]
    ! [ 0 -1 ]
    call dlaset('A', nmat, nmat, 0.0d0, 1.0d0, b_mat(1,1), 2*nmat)
    call dlaset('A', nmat, nmat, 0.0d0, -1.0d0, b_mat(nmat+1,nmat+1), 2*nmat)
    ! FHJ: Call *evx in query mode
    abstol = 0.0d0
    allocate(work (10))
    call dsygvx &
      (2, 'V', range, 'U', 2*nmat, b_mat, 2*nmat, hbse, 2*nmat, &
      0d0, 0d0, 1, neig, abstol, nfound, &
      evals_t, evecs_r, 2*nmat, work, -1, &
      iwork, ifail, info)
    if (info/=0) then
      !if(peinf%inode==0) then
        write(0,'(/,a,i0,/)') 'ERROR: Query mode for ???gvx failed with info=', info
      !endif
      call die("???gvx query failed", only_root_writes=.true.)
    endif
    lwork = max(1,int(work(1))) + 10*nmat
    if(allocated(work))then;deallocate(work);endif
    allocate(work (lwork))
    ! FHJ: Call p*evx for realz
    if (peinf%inode==0) write(6,'(1x,a,i0)') 'Beginning LAPACK diagonalization. Size: ', 2*nmat
    call dsygvx &
      (2, 'V', range, 'U', 2*nmat, b_mat, 2*nmat, hbse, 2*nmat, &
      0d0, 0d0, 1, neig, abstol, nfound, &
      evals_t, evecs_r, 2*nmat, work, lwork, &
      iwork, ifail, info)
    if (peinf%inode==0) write(6,*) 'Done LAPACK diagonalization'
    if(allocated(work))then;deallocate(work);endif
    if(nfound<neig) then
      if (peinf%inode==0) then
        write(0,'(/,a)') 'ERROR: Diagonalization with zhegvx/dsygvx failed:'
        write(0,'(3(a,i0),/)') 'only ', nfound, ' out of ', neig, ' eigenvalues found.'
      endif
        write(0,'(/,a,i0,/)') 'ERROR: Diagonalization with ???evx failed with info=',info
        write(0,'(3(a,i0),/)') 'only ', nfound, ' out of ', neig, ' eigenvalues found.'
      call die("LAPACK found wrong number of eigenvalues", only_root_writes=.true.)
    endif
    if (info/=0.and.info/=2) then
      !if(peinf%inode==0) then
        write(0,'(/,a,i0,/)') 'ERROR: Diagonalization with ???evx failed with info=',info
      !endif
      call die("p???evx diagonalization failed", only_root_writes=.true.)
    endif
    ! FHJ: Copy eigenvalues/eigenvectors
    evals(1:neig) = evals_t(1:neig)
    if (peinf%inode==0) write(6,*) 'Calculating overlap matrix'
    ! FHJ: Calculate overlap S = evecs_r evecs_r^H.
    ! Note: other codes/works define S by evecs_r^H evecs_r
    call dsyrk&
      ('L', 'N', 2*nmat, 2*nmat, 1.0d0, evecs_r, 2*nmat, 0.0d0, s_mat, 2*nmat)
    if (peinf%inode==0) write(6,*) 'Performing Cholesky decomposition'
    ! FHJ: Invert overlap matrix. First, do Cholesky decomposition.
    call dpotrf('L', 2*nmat, s_mat, 2*nmat, info)
    if (info/=0) then
      if (peinf%inode==0) write(0,*) 'ERROR: got info=', info
      call die('Cholesky decomposition failed', only_root_writes=.true.)
    endif
    if (peinf%inode==0) write(6,*) 'Inverting overlap matrix'
    ! FHJ: Now, invert the matrix
    call dpotri('L', 2*nmat, s_mat, 2*nmat, info)
    if (info/=0) then
      if (peinf%inode==0) write(0,*) 'ERROR: got info=', info
      call die('matrix inversion failed', only_root_writes=.true.)
    endif
    if (peinf%inode==0) write(6,*) 'Calculating left evecs'
    ! FHJ: Multiply S^{-1} by evecs_r to get evecs_l^H
    call dsymm&
      ('L', 'L', 2*nmat, 2*nmat, 1.0d0, s_mat, 2*nmat, &
      evecs_r, 2*nmat, 0.0d0, evecs_l, 2*nmat)
    ! FHJ: Apply complex conjugation to get evecs_l^T
    evecs_l = (evecs_l)
    if (peinf%inode==0) write(6,*) 'Diagonalization done'
   
  end subroutine serial_full_gvx
end subroutine diagonalize
!> FHJ: Internal subroutine that actually prints ||M - M^H||_F / ||M + M^H||_F.
subroutine print_broken_herm(norm, block)
  real(DP), intent(in) :: norm
  character, intent(in) :: block
  logical, save :: warned=.false.
 
  if (peinf%inode==0) then
    write(6,'(1x,a)',advance='no') 'Degree of broken Hermiticity of the '
    if (block=='A'.or.block=='B') then
      write(6,'(a)',advance='no') 'subblock '//block//' of the BSE Hamiltonian: '
    else
      write(6,'(a)',advance='no') 'BSE Hamiltonian: '
    endif
    if (norm>=1.and..not.warned) then
      write(6,'(f0.4,a)') norm, ' %'
      write(0,*)
      write(0,*) 'WARNING: non-Hermiticity of the BSE Hamiltonian is large!'
      write(0,*) 'Some possible reasons include:'
      write(0,*) '- Coarse grid that is too coarse'
      write(0,*) '- Large q0 shift (>0.001)'
      write(0,*) '- Small dielectric cutoff (<8 Ry)'
      write(0,*)
      warned = .true.
    else
      write(6,'(f6.4,a)') norm, ' %'
    endif
  endif
 
end subroutine print_broken_herm
!> FHJ: Calculate the degree of broken Hermiticity: ||M - M^H||_F / ||M + M^H||_F.
subroutine calc_broken_herm(mat, nmat, block)
  integer, intent(in) :: nmat
  real(DP), intent(in) :: mat(nmat,nmat)
  character, intent(in) :: block
  real(DP) :: mat2(nmat,nmat)
  real(DP) :: norm, norm2
  real(DP), external :: dlange
 
  ! Calc ||A - A^H||_F
  mat2(:,:) = transpose(mat(:,:))
  if (block=='B') then
    mat2(:,:) = mat(:,:) - mat2(:,:)
  else
    mat2(:,:) = mat(:,:) - (mat2(:,:))
  endif
  norm = dlange('F', nmat, nmat, mat2, nmat, mat2)
  ! Calc ||A + A^H||_F
  mat2(:,:) = transpose(mat(:,:))
  if (block=='B') then
    mat2(:,:) = mat(:,:) + mat2(:,:)
  else
    mat2(:,:) = mat(:,:) + (mat2(:,:))
  endif
  norm2 = dlange('F', nmat, nmat, mat2, nmat, mat2)
  norm = norm / norm2 * 1d2
  call print_broken_herm(norm, block)
 
end subroutine calc_broken_herm
end module diagonalize_m
