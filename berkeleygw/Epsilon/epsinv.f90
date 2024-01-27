!=================================================================================
!
! Routines:
!
! (1) epsinv() Originally By (?) Last Modified 5/1/2008 (JRD)
!
! This routine:
!
! 1. Calculates epsilon based on chi.
! 2. Inverts epsilon.
! 3. Writes the result to unit=12 if q0="zero" and unit 13 otherwise.
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
module epsinv_m
  use global_m
  use inversion_m
  use io_utils_m
  use misc_m
  use scalapack_m
  use timing_m, only: timing => epsilon_timing
  use vcoul_generator_m
  use write_matrix_m
  implicit none
  private
  public :: &
    epsinv
contains
subroutine epsinv(gvec,pol,ekin,q0,is_q0,crys,scal,kp,omega_plasma,iq,E_rpa)
  type (gspace), intent(in) :: gvec
  type (polarizability), intent(in) :: pol
  real(DP), intent(in) :: ekin(gvec%ng)
  real(DP), intent(in) :: q0(3)
  logical, intent(in) :: is_q0
  type (crystal), intent(in) :: crys
  type (scalapack), intent(in) :: scal
  type (kpoints), intent(in) :: kp
  real(DP), intent(in) :: omega_plasma
  integer, intent(in) :: iq
  real(DP), optional, intent(out) :: E_rpa
  integer :: qgrid(3)
  real(DP) :: q0vec(3)
  type (twork_scell) :: work_scell
  integer :: i,j,iw,jj,ii,is,js,iunit,my_iq,ifreq
  integer :: irow, icol, icurr, irowm, icolm, icountr, icountc
  integer :: ig_l, ig_g, igp_l, igp_g
  integer, allocatable :: isorti(:)
  integer :: iscreen, nfq, iparallel, isize,ifreq_para,freq_grp_ind
  real(DP) :: vc, oneoverq, avgcut
  real(DP) :: epssum1R, epssum2R, epssum1R_rel,epssum2R_rel
  real(DP) :: chitmp
  real(DP), allocatable :: eps(:,:),ewng(:)
  real(DP), allocatable :: epsTemplate(:,:)
  real(DP), allocatable :: epsdiag(:,:,:),epsdiagt(:,:,:), vcoultemp(:)
  real(DP), allocatable :: vcoul(:)
  complex(DPC), allocatable :: chiRDyntmp(:)
  complex(DPC), allocatable :: epsRDyn(:,:,:), epsADyn(:,:,:)
  complex(DPC), allocatable :: epsRDyn_head(:), epsRDyn_head_temp(:)
  ! Auxiliary matrix for inversion
  complex(DPC), allocatable :: eps1Aux(:,:)
  ! Auxiliary matrices for subspace truncation method
  complex(DPC), allocatable :: epsRDyn_sub(:,:,:), eps1Aux_sub(:,:)
  logical :: subspace, keep_full_eps_static
  integer :: nrow_loc_sub, ncol_loc_sub, neig_sub
  real(DP) :: epsheaddummy, wcoul0
  character*80 :: filename
  type(progress_info) :: prog_info
  type (scalapack) :: scal_sub, scal_aux
  real(DP), allocatable :: integ_rpa_pol(:)
  integer :: freq_offset
  ! for timing (temporary)
  real(DP) :: t1, t2
  ! to avoid large freq allocation in the subspace case
  integer :: nfreq_group_alloc
 
  allocate(vcoul (pol%nmtx))
  subspace = .FALSE.
  if(pol%freq_dep==0) then
    allocate(eps (scal%npr,scal%npc))
    allocate(ewng (pol%nmtx))
  else
    allocate(chiRDyntmp (pol%nfreq_in_group))
    !XXXX
    nfreq_group_alloc = pol%nfreq_in_group
    if (subspace) then
      if ( pol%matrix_in_subspace_basis ) then
         nfreq_group_alloc = 1
      end if
    end if
    allocate(epsRDyn (scal%npr,scal%npc,nfreq_group_alloc))
    !XXXX
    allocate(eps1Aux (scal%npr,scal%npc))
    ! subspace truncation specific stuff
    IF(subspace) THEN
      nrow_loc_sub = pol%nrow_local_sub
      ncol_loc_sub = pol%ncol_local_sub
      neig_sub = pol%neig_sub
      allocate(epsRDyn_sub (MAX(1,nrow_loc_sub),MAX(1,ncol_loc_sub),pol%nfreq_in_group))
      allocate(eps1Aux_sub (MAX(1,nrow_loc_sub),MAX(1,ncol_loc_sub)))
      epsRDyn_sub = (0.d0,0.d0)
    END IF
  endif
  allocate(isorti (gvec%ng))
!------------------------------
! Invert isrtx
!
! SIB: isorti is the inverse sort order for pol%isrtx.
! pol%isrtx has the sort indices for |q0+gvec%components|^2
!
! JRD XXX is the initialization necessary. probably very slow
  if (pol%freq_dep==0) then
    eps(:,:)=0.0d0
  else
    epsRDyn(:,:,:)=(0.d0,0.d0)
  endif
  vcoul(:)=0.0d0
  do i=1,gvec%ng
    isorti(pol%isrtx(i)) = i
  end do
!-------------- Construct Dielectric Matrix ---------------------------
!
! e(q+g,q+g`) = del(g,g`) - (8pi/(q+g)**2) chi(q+g,q+g`). For spin-polarized
! calc., e(G+q,G`+q)=del(G,G`)- (8PI/(G+q)^2) SUM_spin chi(G+q,G`+q,ispin)
! Which is pol%chi(j,1) as compiled in epsilon_main.f90. If q--> 0 , we have to treat
! the wings separately
!
! SIB: using the Rydberg as our unit of energy
! if pol%icutv is on (coulomb cutoff) then we multiply the coulomb
! interaction by the appropriate factor (1-cos(vcut*|q0+g|))
!
! if (peinf%inode .eq. 0) then
! write(6,*) ' '
! write(6,*) 'Calculating Coulomb Potential'
! endif
  icurr=0
! Generator Coulomb Interaction Array Vcoul
! For Epsilon, We want to treat all types of screening the same for vcoul. Because
! we calculate it exactly
  avgcut=TOL_ZERO
  iscreen=0
  ! FHJ: The following are not used when peinf%jobtypeeval==1
  nfq = 0
  q0vec = 0d0
  iparallel = 1
  ! FHJ: FIXME - supercell truncation actually cares about qgrid.
  ! But is it doing the correct thing?
  qgrid(:) = 1
  call timing%start(timing%epsinv_vcoul)
  epsheaddummy=0.0d0
  wcoul0=0.0d0
  IF(subspace) THEN
    ! we already have vcoul
    vcoul(:) = pol%vcoul_sub(:)
  ELSE
    if(pol%nfreq_group .eq. 1) then
      call vcoul_generator(pol%icutv,pol%truncval,gvec,crys%bdot,crys%celvol, &
        nfq,pol%nmtx,pol%isrtx,iscreen,q0,q0vec,vcoul, &
        pol%iwritecoul,iparallel,avgcut,oneoverq,qgrid,epsheaddummy, &
        work_scell,.false.,wcoul0)
    else
      call vcoul_generator(pol%icutv,pol%truncval,gvec,crys%bdot,crys%celvol, &
        nfq,pol%nmtx,pol%isrtx,iscreen,q0,q0vec,vcoul, &
        pol%iwritecoul,iparallel,avgcut,oneoverq,qgrid,epsheaddummy, &
        work_scell,.false.,wcoul0,nfreq_group=pol%nfreq_group)
    endif
  END IF
  call timing%stop(timing%epsinv_vcoul)
! write(6,*) 'Done VCoul'
  if (pol%freq_dep==0) then
    do i=1,pol%nmtx
      irow=MOD(INT(((i-1)/scal%nbl)+TOL_SMALL),scal%nprow)
      if(irow.ne.scal%myprow) cycle
      vc = vcoul(i)
!-------------------------
! Construct non-interacting epsilon
! Static case
      do j=1,pol%nmtx
        icol=MOD(INT(((j-1)/scal%nbl)+TOL_SMALL),scal%npcol)
        if(icol .eq. scal%mypcol) then
          icurr=icurr+1
          irowm=INT((icurr-1)/scal%npc+TOL_SMALL)+1
          icolm=MOD((icurr-1),scal%npc)+1
          eps(irowm,icolm) = 0.0d0
          if (i.eq.j) eps(irowm,icolm) = 1.0d0
          chitmp = pol%chi(irowm,icolm,1)
          eps(irowm,icolm) = eps(irowm,icolm) - vc*chitmp
        endif
      end do
    end do
  else !freq_dep>0 below
    IF(subspace) THEN
      ! at this point we have the symmetrized polarizability
      ! (already scaled by the SQRT of the coulomb potential)
      ! for all frequency defined within the subspace of omega=0,
      ! additionally for omega=0 we have the full (original basis) polarizability
      ! we calculate epsilon just changing sing of chiRDyn and adding one on the diagonal
      epsRDyn_sub(:,:,:) = - pol%chiRDyn(:,:,:,1)
      icountr = 0
      icountc = 0
      DO i = 1, neig_sub
        icol = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL),scal%npcol)
        irow = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL),scal%nprow)
        IF(irow == scal%myprow) icountr = icountr + 1
        IF(icol == scal%mypcol) icountc = icountc + 1
        IF(irow /= scal%myprow) CYCLE
        IF(icol /= scal%mypcol) CYCLE
        epsRDyn_sub(icountr,icountc,:) = epsRDyn_sub(icountr,icountc,:) + 1.0D+00
      END DO
      ! same for omega=0 full epsilon
      epsRDyn(:,:,1) = - pol%chiRDyn_sym_omega0(:,:)
      icountr = 0
      icountc = 0
      DO i=1, pol%nmtx
        icol = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL),scal%npcol)
        irow = MOD(INT(((i-1)/scal%nbl)+TOL_SMALL),scal%nprow)
        IF(irow == scal%myprow) icountr = icountr + 1
        IF(icol == scal%mypcol) icountc = icountc + 1
        IF(irow /= scal%myprow) CYCLE
        IF(icol /= scal%mypcol) CYCLE
        epsRDyn(icountr,icountc,1) = epsRDyn(icountr,icountc,1) + 1.0D+00
      END DO
    ELSE
      allocate(epsTemplate (scal%npr,scal%npc))
      allocate(vcoultemp (scal%npr))
      epsTemplate = 0D0
! JRD Set diagonal elements to 1 and define vcoultemp
      icountc = 0
      icountr = 0
      do i=1,pol%nmtx
        icol=MOD(INT(((i-1)/scal%nbl)+TOL_SMALL),scal%npcol)
        irow=MOD(INT(((i-1)/scal%nbl)+TOL_SMALL),scal%nprow)
        if(icol.eq.scal%mypcol) icountc = icountc+1
        if(irow.eq.scal%myprow) icountr = icountr+1
        if(irow.ne.scal%myprow) cycle
        vcoultemp(icountr) = vcoul(i)
        if(icol.ne.scal%mypcol) cycle
        epsTemplate(icountr,icountc) = 1D0
      enddo
! JRD XXX if we block this we can probably hold epsTemplate in cache
      do ifreq = 1, pol%nfreq_in_group
        do i=1,scal%npc
            epsRDyn(:,i,ifreq) = epsTemplate(:,i)- &
              vcoultemp(:)*pol%chiRDyn(:,i,ifreq,1)
        end do
      end do
      if(allocated(epsTemplate))then;deallocate(epsTemplate);endif
      if(allocated(vcoultemp))then;deallocate(vcoultemp);endif
    END IF
  endif
3999 format(1x,a,i6,a,2es25.15e3)
  if (peinf%inode==0 .and. pol%freq_dep==0) then
    write(6,3999) 'q-pt ', iq, ': Head of Epsilon         = ', eps(1,1)
    write(6,3999) 'q-pt ', iq, ': Epsilon(2,2)            = ', eps(2,2)
  endif
  if (peinf%inode==0 .and. pol%freq_dep>0) then
    write(6,3999) 'q-pt ', iq, ': Head of Epsilon         = ', epsRDyn(1,1,1)
    write(6,3999) 'q-pt ', iq, ': Epsilon(2,2)            = ', epsRDyn(2,2,1)
  endif
!-------------------------------------------------------------
! Print head versus frequency
  if(pol%nfreq_group .eq. 1 .and. peinf%rank_f .eq. 0 .and. pol%freq_dep .eq. 2) then
    allocate(epsRDyn_head (pol%nfreq))
    epsRDyn_head=cmplx(0d0,0d0,kind=DPC)
    IF(subspace) THEN
      epsRDyn_head(:)=epsRDyn_sub(1,1,:)
    ELSE
      epsRDyn_head(:)=epsRDyn(1,1,:)
    END IF
  elseif(pol%nfreq_group > 1 .and. peinf%inode .lt. peinf%npes .and. peinf%rank_f .eq. 0 .and. pol%freq_dep .eq. 2) then
    allocate(epsRDyn_head (pol%nfreq))
    allocate(epsRDyn_head_temp (pol%nfreq))
    epsRDyn_head=cmplx(0d0,0d0,kind=DPC)
    epsRDyn_head_temp=cmplx(0d0,0d0,kind=DPC)
! JRD XXX this nonsense may be slow
    do ifreq=1,pol%nfreq
      freq_grp_ind=mod(ifreq-1,pol%nfreq_group)
      ifreq_para=(ifreq+pol%nfreq_group-1)/pol%nfreq_group
      if(freq_grp_ind .eq. peinf%rank_mtxel) then
        epsRDyn_head_temp(ifreq)=epsRDyn(1,1,ifreq_para)
      endif
    enddo
  endif
  if (peinf%rank_mtxel .eq. 0 .and. peinf%rank_f .eq. 0 .and. pol%freq_dep .eq. 2) then
    write(52,'("# q= ",3f12.6," nmtx= ",i0)') q0(:), pol%nmtx
    write(52,*)
    do iw = 1, pol%nfreq
      write(52,'(f12.6,4(1x,es13.6))') pol%dFreqGrid(iw), &
        dble(epsRDyn_head(iw)),aimag(epsRDyn_head(iw))
    enddo
  endif
!------------ Here we invert the epsilon matrix -----------------------------
!
  call timing%start(timing%epsinv_invert)
! Serial Version
  if (pol%freq_dep==0) then
    call dinvert_serial(pol%nmtx,eps)
  else
    call progress_init(prog_info, 'inversion of dielectric matrix', 'frequency', pol%nfreq_in_group)
    do iw = 1, pol%nfreq_in_group
      call progress_step(prog_info)
! JRD XXX Copy no longer necessary
      eps1Aux(:,:) = epsRDyn(:,:,iw)
      call zinvert_serial(pol%nmtx,eps1Aux)
      epsRDyn(:,:,iw) = eps1Aux(:,:)
    enddo
    call progress_free(prog_info)
  endif
  call timing%stop(timing%epsinv_invert)
! Done inverting
!-----------------------------------------------------------------------------
  if (peinf%inode==0 .and. pol%freq_dep==0) then
    write(6,3999) 'q-pt ', iq, ': Head of Epsilon Inverse = ', eps(1,1)
    write(6,3999) 'q-pt ', iq, ': Epsilon Inverse(2,2)    = ', eps(2,2)
  endif
  if (peinf%inode==0 .and. pol%freq_dep>0) then
    if(subspace) then
      if((.not. pol%keep_full_eps_static) .and. pol%matrix_in_subspace_basis) then
        write(6,'(1x,A)') 'Within the subspace basis:'
      end if
    end if
    write(6,3999) 'q-pt ', iq, ': Head of Epsilon Inverse = ', epsRDyn(1,1,1)
    write(6,3999) 'q-pt ', iq, ': Epsilon Inverse(2,2)    = ', epsRDyn(2,2,1)
  endif
!----------- Print out independent matrix elements ---------------------------
! JRD XXX this nonsense may be slow
  if(pol%nfreq_group .eq. 1 .and. peinf%rank_f .eq. 0 .and. pol%freq_dep .eq. 2) then
    if(nfreq_group_alloc == pol%nfreq_in_group) then
      epsRDyn_head(:)=epsRDyn(1,1,:)
    endif ! (nfreq_group_alloc == pol%nfreq_in_group)
  elseif(pol%nfreq_group > 1 .and. peinf%inode .lt. peinf%npes .and. peinf%rank_f .eq. 0 .and. pol%freq_dep .eq. 2) then
    if(nfreq_group_alloc == pol%nfreq_in_group) then
      epsRDyn_head=cmplx(0d0,0d0,kind=DPC)
      epsRDyn_head_temp=cmplx(0d0,0d0,kind=DPC)
      do ifreq=1,pol%nfreq
        freq_grp_ind=mod(ifreq-1,pol%nfreq_group)
        ifreq_para=(ifreq+pol%nfreq_group-1)/pol%nfreq_group
        if(freq_grp_ind .eq. peinf%rank_mtxel) then
          epsRDyn_head_temp(ifreq)=epsRDyn(1,1,ifreq_para)
        endif
      enddo
    endif ! (nfreq_group_alloc == pol%nfreq_in_group)
  endif
  if (peinf%inode==0 .and. pol%freq_dep==2) then
    write(51,'("# q= ",3f12.6," nmtx= ",i0)') q0(:), pol%nmtx
    write(51,*)
    do iw = 1, pol%nfreq
      write(51,'(f12.6,4(1x,es15.6e3))') pol%dFreqGrid(iw), &
        dble(epsRDyn_head(iw)),aimag(epsRDyn_head(iw))
    enddo
  endif
  if (peinf%inode.eq.0 .and. pol%freq_dep .eq. 0) then
    ! JRD Warn User about possible lack of symmetry
    if (is_q0) then
      write(7,*)
      write(7,*) 'For q0 points, you should check the symmetry (eps(G,G'') = eps*(-G,-G'')) by'
      write(7,*) 'using the eps0sym code. Wavefunction convergence, as well as a finite q-shift'
      write(7,*) 'may cause this property of eps(G,G'') to be broken.'
      write(6,*)
      write(6,*) 'For q0 points, you should check the symmetry (eps(G,G'') = eps*(-G,-G'')) by'
      write(6,*) 'using the eps0sym code. Wavefunction convergence, as well as a finite q-shift'
      write(6,*) 'may cause this property of eps(G,G'') to be broken.'
    endif
    write(7,4000) kp%nspin
    do i=1,scal%npr
      is=scal%isrtxrow(i)
      do j=1,scal%npc
        js=scal%isrtxcol(j)
        if (i .eq. j .or. i .eq. j+1) then
          write(7,4200) gvec%components(1:3,is), gvec%components(1:3,js), eps(i,j)
        endif
      end do
    end do
  end if
! JRD the i and j loops are out of order below
  if (peinf%inode==0 .and. pol%freq_dep>0) then
    write(7,4001) kp%nspin
    do iw = 1, min(pol%nfreq_in_group, nfreq_group_alloc)
      do i=1,scal%npr
        is=scal%isrtxrow(i)
        do j=1,scal%npc
          js=scal%isrtxcol(j)
          if (i .eq. j .or. i .eq. j+1) then
            write(7,4300) gvec%components(1:3,is), gvec%components(1:3,js), epsRDyn(i,j,iw)
          endif
        end do
      end do
    end do
  end if
4000 format(/ /,13x,'g',19x,'gp',9x, &
       'inverse epsilon           nspin= ',1i1)
4001 format(/ /,13x,'g',19x,'gp',9x, &
  'inverse epsilon RDyn      nspin= ',1i1)
4200 format(5x,3i5,5x,3i5,5x,2f13.8)
4300 format(5x,3i5,5x,3i5,5x,4f13.8)
!---------- Full-Frequency Sum-Rule -----------------------------------------
  epssum1R=0D0
  epssum2R=0D0
  ! FHJ: These sum rules are *wrong* for systems without inversion symmetry.
  ! The sum rules apply to the Hermitian/Antihermitian components of the
  ! dielectric matrix. See D. J. Johnson, Phys. Rev. B 9, 4475 (1974),
  ! eqns. 2.1 to 2.7. Note that the derivation of the complex GPP is right,
  ! even if the notation is misleading (see comment in LITERATURE.html).
  ! Anyone care to fix this??
  if (peinf%inode==0 .and. pol%freq_dep==2 .and. pol%freq_dep_method/=2) then
    do iw = 2, pol%nFreq-pol%nfreq_imag
      epssum1R=epssum1R+(1D0*Ryd/pol%dFreqGrid(iw))* &
        aimag(epsRDyn_head(iw))*(pol%dFreqGrid(iw)-pol%dFreqGrid(iw-1))/Ryd
      epssum2R=epssum2R+(pol%dFreqGrid(iw)/Ryd)* &
        aimag(epsRDyn_head(iw))*(pol%dFreqGrid(iw)-pol%dFreqGrid(iw-1))/Ryd
    enddo
    epssum1R=(2D0*epssum1R/Pi_D)+1D0
    epssum1R_rel=(epssum1R)/dble(epsRDyn(1,1,1))
    epssum2R_rel=(-1D0*epssum2R)/((Pi_D/2D0)*omega_plasma**2)
! Ref: Hybertsen & Louie PRB 34, 5390 (1986), eq. 29 and Appendix A
    write(6,*) ' '
    write(6,*) 'Full Frequency: Sum rules for head:'
    write(6,*) 'Int((1/w)*Im(eps^-1(w))) =', epssum1R_rel*100D0, ' % of exact'
    write(6,*) 'Int((w)*Im(eps^-1(w))) =', epssum2R_rel*100D0, ' % of exact'
  endif
!---------- Write inverse dielectric matrices to file -----------------------
  call timing%start(timing%epsinv_i_o)
  if (peinf%inode==0) write(6,'(/1x,a)') 'Writing dielectric matrix to file'
  !XXXXXXXX
  if (subspace) then
    if (keep_full_eps_static) then
      if (peinf%inode==0) write(6,'(1x,a/)') 'Subspace: Full Epsilon will be retained for omega=0'
    else
      if (peinf%inode==0) write(6,'(1x,a/)') 'Subspace: Full Epsilon will NOT be retained for omega=0'
    end if
  end if
  !XXXXXXXX
  if (is_q0) then
    filename = 'eps0mat.h5'
  else
    filename = 'epsmat.h5'
  endif
    iunit=13
    if (is_q0) iunit=12
    if (peinf%inode .eq. 0) then
      write(iunit) gvec%ng,pol%nmtx, &
        (pol%isrtx(i),isorti(i),i=1,gvec%ng)
      write(iunit) (ekin(i),i=1,gvec%ng)
      write(iunit) (q0(i),i=1,3)
    endif
  if (pol%freq_dep==0) then
      call write_matrix_d(scal, eps, pol%nmtx, iunit)
  else !freq_dep/=0 below
    IF(subspace .AND. pol%matrix_in_subspace_basis) THEN
      ! here we write out the subspace epsinv matrices
        ! here no HDF5
        ! write the coulomb potential
        if (peinf%inode .eq. 0) then
           ! already written:
           ! gvec%ng, pol%nmtx, pol%isrtx, isorti, ekin, q0
           ! Write coulomb potential
           write(iunit) (vcoul(i), i=1, pol%nmtx)
        end if
        ! if we want to use the full epsinv at omega zero write it out here
        IF(keep_full_eps_static) THEN
          if (peinf%inode==0) write(6,'(1x,a)') 'Writing Full Inv Epsilon for omega=0'
          eps1Aux(:,:) = epsRDyn(:,:,1)
          call write_matrix_d_sub(scal, eps1Aux, pol%nmtx, iunit)
        END IF
        ! write eigenvectors (just consider it as a square matrix)
        if (peinf%inode==0) write(6,'(1x,a)') 'Writing Eigenvectors of Epsilon omega=0'
        ! write the basis size here
        if (peinf%inode==0) then
          write(iunit) neig_sub
        end if
        eps1Aux(:,:) = pol%eigenvect_omega0(:,:)
        call write_matrix_d_sub(scal, eps1Aux, pol%nmtx, iunit, neig_sub)
        ! and all subspace matrices (including that for omega=zero),
        ! we don`t write the advace, just recalculate it in sigma as done in the case of use_hdf5
        ! create the sub scalapack environment
        if (peinf%inode==0) write(6,'(1x,a)') 'Writing Subspace Inv Epsilon for all other frequencies'
        if (peinf%inode==0) write(6,'(1x,a,i5)') 'Number of frequencies:', pol%nFreq
        scal_sub%nprow = scal%nprow
        scal_sub%npcol = scal%npcol
        scal_sub%myprow = scal%myprow
        scal_sub%mypcol = scal%mypcol
        scal_sub%nbl = scal%nbl
        scal_sub%icntxt = scal%icntxt
        scal_sub%npr = nrow_loc_sub
        scal_sub%npc = ncol_loc_sub
        call write_matrix_f(scal_sub, pol%nFreq, epsRDyn_sub, neig_sub, iunit, pol%nfreq_group)
    ELSE ! else (write subspace matrices)
        call write_matrix_f(scal, pol%nFreq, epsRDyn, pol%nmtx, iunit, pol%nfreq_group &
          )
    END IF ! write subspace matrices
  endif
  if (peinf%inode==0) write(6,'(1x,a/)') 'Ok'
  call timing%stop(timing%epsinv_i_o)
! Finished writing eps
!------------------------------------------------------------------------
  if(allocated(vcoul))then;deallocate(vcoul);endif
  if(allocated(isorti))then;deallocate(isorti);endif
  if(pol%freq_dep==0) then
    if(allocated(eps))then;deallocate(eps);endif
    if(allocated(ewng))then;deallocate(ewng);endif
  else
    if(allocated(epsRDyn))then;deallocate(epsRDyn);endif
    if(allocated(eps1Aux))then;deallocate(eps1Aux);endif
    if(allocated(chiRDyntmp))then;deallocate(chiRDyntmp);endif
  endif
  if (pol%freq_dep==2) then
    if (pol%nfreq_group==1.and.peinf%rank_f==0) then
      if(allocated(epsRDyn_head))then;deallocate(epsRDyn_head);endif
    elseif (pol%nfreq_group>1.and.peinf%inode<peinf%npes.and.peinf%rank_f==0) then
      if(allocated(epsRDyn_head))then;deallocate(epsRDyn_head);endif
      if(allocated(epsRDyn_head_temp))then;deallocate(epsRDyn_head_temp);endif
    endif
  endif
 
  return
end subroutine epsinv
end module epsinv_m
