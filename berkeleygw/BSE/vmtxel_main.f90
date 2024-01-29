!===============================================================================
!
! Routines:
!
! (1) vmtxel Originally by GKA (2018)
!
! Read in the wavefunctions and compute the velocity (or momentum) operator
! matrix elements for vertical transitions, unrestricted by the occupations.
!
! The program writes the matrix elements in the file 'vmtxel'
! or 'vmtxel_b1', 'vmtxel_b2', 'vmtxel_b3'
! Note however that, in these files, the ordering of the bands are always
! bottom up, for both initial and final states, at variance with the
! absorption executable, which counts valence bands from top to bottom.
!
!===============================================================================

program vmtxel
  use fullbz_m, only: dealloc_grid
  use genwf_m
  use global_m
  use input_fi_m
  use input_q_m
  use timing_m, only: timing => bse_timing
  use vmtxel_m
  use wfn_rho_vxc_io_m
  use write_program_header_m
  implicit none
  type (xctinfo) :: xct
  type (eqpinfo) :: eqp
  type (flags) :: flag
  type (crystal) :: crys
  type (gspace) :: gvec
  type (grid) :: kg_fi, kgq_fi
  type (kpoints) :: kp, kp_fi, kpq_fi
  type (symmetry) :: syms
  type (wavefunction) :: wfnc_fi
  type (wavefunction) :: wfnvq_fi
  type (work_genwf) :: work, workq
  type (int_wavefunction) :: intwfnc
  type (int_wavefunction) :: intwfnv
  type (vmtxel_t) :: dip
  character(len=3) :: sheader
  real(DP) :: omega_plasma
  integer :: iflavor
  integer :: is,ik,ikq,ic,iv,jdim,ipol,ikt
  integer :: error
  logical :: found_wfnq
  integer, allocatable :: indexq_fi(:)
  call peinfo_init()
  !----------------------------------------------------------------------------
  ! Initialization
  ! Write header
  call write_program_header('BSE/vmtxel', .false.)
  ! Initialize HDF5
  ! Initialize timer
  call timing%init()
  call timing%start(timing%total)
  ! Initialize some flags by hand
  flag%bz0 = 0 ! Use symmetryes in unshifted grid
  flag%bzq = 1 ! Do not use symmetries in shifted grid
  flag%vm = 0 ! Calculate velocity/momentum matrix elements
  flag%read_dtmat = .false.
  xct%rfermi = .true. ! Do not correct the fermi level
  xct%efermi = 0.0_dp !!
  xct%eqp_corrections = .false.
  xct%unrestricted_transf = .true.
  xct%vmin = 1 ! do not specify occupations
  xct%vmax = 0 !!
  ! avoid using HDF5 for now
  xct%use_wfn_hdf5 = .false.
  xct%wfn_hdf5_min_band_block = -1
  ! Check whether the WFNq_fi file exists
  inquire(file='WFNq_fi', exist=found_wfnq)
  if (found_wfnq) then
    ! Use velocity operator
    flag%opr = 0 ! Use velocity operator
    xct%npol = 1 ! GKA: Should check whether WFNq_fi is avalable
  else
    ! Use momentum operator
    flag%opr = 1 ! Use momentum operator
    xct%npol = 3 ! GKA: Should check whether WFNq_fi is avalable
  endif
  !----------------------------------------------------------------------------
  ! Read the header of WFN and extract some dimensions
  ! The following section is copied from bse_init
  ! Read the header of WFN
  if (peinf%inode == 0) &
    call open_file(25, file='WFN_fi', form='unformatted', status='old')
  sheader = 'WFN'
  iflavor = 0
  call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys, &
                               dont_warn_kgrid=.true.)
  if (peinf%inode == 0) call close_file(25)
  ! Use dimensions read from the header of WFN_fi
  xct%nspin = kp%nspin
  xct%nvb_fi = minval(kp%ifmax(:,:) - kp%ifmin(:,:)) + 1
  xct%ncb_fi = kp%mnband - xct%nvb_fi
  !----------------------------------------------------------------------------
  ! Read in wavefunctions
  call logit('Calling input')
  call timing%start(timing%input)
  call input_fi(crys,gvec,kg_fi,kp_fi,syms,eqp,xct,flag, &
                omega_plasma,.true.,intwfnc,read_all_bands=.true.)
  ! Print out some info
  if (peinf%inode.eq.0) then
    write(6,'(/1x,a)') 'Info on the wavefunctions:'
    write(6,'(1x,a,i0)') '- Number of valence bands: ', xct%nvb_fi
    write(6,'(1x,a,i0)') '- Number of cond. bands: ', xct%ncb_fi
    write(6,'(1x,a,i0)') '- Number of spins: ', xct%nspin
    write(6,'()')
  endif
  call timing%stop(timing%input)
  allocate(indexq_fi (xct%nkpt_fi))
  allocate(xct%indexq_fi (xct%nkpt_fi))
  if (flag%vm.ne.1.or. .not. flag%read_dtmat) then
    call timing%start(timing%input_q)
    call logit('Calling input_q')
    call input_q(kp_fi,crys,gvec,kg_fi,kgq_fi,kpq_fi,syms,xct,indexq_fi,eqp,&
                 flag,intwfnv,read_all_bands=.true.)
    call timing%stop(timing%input_q)
  endif
  !----------------------------------------------------------------------------
  ! Calculate the velocity (or momentum) matrix elements
  call logit('Calculating v/p matrix elememts')
  ! Initialize
  call dip%init(xct%nspin, xct%nkpt_fi, kp_fi%mnband, kpq_fi%mnband, &
                opr=flag%opr, npol=xct%npol, band_ordering=1, &
                with_velocity=.true.)
  ! Allocate memory
  call dip%alloc()
  ! Copy list of k-points
  do ik=1, dip%nk
    dip%kpt(:,ik) = kg_fi%f(:,ik)
  end do
  ! Set polarization vector
  dip%pol(:,1) = xct%pol
  call timing%start(timing%vmtxel)
  do ikt=1, peinf%ikt(peinf%inode+1)
    ik = peinf%ik(peinf%inode+1,ikt)
    ikq = indexq_fi(ik)
    call genwf(crys,gvec,kg_fi,syms,wfnc_fi,ik,ik,kp%nspin,kp%mnband,&
               work,intwfnc,1,is_cond=.true.)
    call genwf(crys,gvec,kgq_fi,syms,wfnvq_fi,ik,ikq,kpq_fi%nspin,kpq_fi%mnband,&
               workq,intwfnv,1,is_cond=.false.)
    call dip%compute_ik_vmtxel(ik, wfnc_fi, wfnvq_fi, gvec, xct%qshift, crys, eqp)
    if(associated(wfnc_fi%cg))then;deallocate(wfnc_fi%cg);nullify(wfnc_fi%cg);endif
    if(associated(wfnc_fi%isort))then;deallocate(wfnc_fi%isort);nullify(wfnc_fi%isort);endif
    if(associated(wfnvq_fi%cg))then;deallocate(wfnvq_fi%cg);nullify(wfnvq_fi%cg);endif
    if(associated(wfnvq_fi%isort))then;deallocate(wfnvq_fi%isort);nullify(wfnvq_fi%isort);endif
  enddo
  ! typedefs initializes all of these ikolds to 0
  if(work%ikold.ne.0) then
    if(associated(work%cg))then;deallocate(work%cg);nullify(work%cg);endif
    if(associated(work%ph))then;deallocate(work%ph);nullify(work%ph);endif
    if(associated(work%ind))then;deallocate(work%ind);nullify(work%ind);endif
    if(associated(work%isort))then;deallocate(work%isort);nullify(work%isort);endif
  endif
  if(workq%ikold.ne.0) then
    if(associated(workq%cg))then;deallocate(workq%cg);nullify(workq%cg);endif
    if(associated(workq%ph))then;deallocate(workq%ph);nullify(workq%ph);endif
    if(associated(workq%ind))then;deallocate(workq%ind);nullify(workq%ind);endif
    if(associated(workq%isort))then;deallocate(workq%isort);nullify(workq%isort);endif
  endif
  ! Share matrix elements
  call dip%reduce()
  ! Write to file
  call dip%write_vmtxel()
  call timing%stop(timing%vmtxel)
  !----------------------------------------------------------------------------
  ! Free memory
  call dip%free()
  call kp%free()
  call kp_fi%free()
  call kpq_fi%free()
  if (flag%vm.ne.1.or. .not. flag%read_dtmat) then
    call dealloc_grid(kgq_fi)
  endif
  if (flag%vm == 0 .and. .not. flag%read_dtmat) then
    if(associated(intwfnc%cgk))then;deallocate(intwfnc%cgk);nullify(intwfnc%cgk);endif
    if(associated(intwfnv%cgk))then;deallocate(intwfnv%cgk);nullify(intwfnv%cgk);endif
    if(associated(intwfnc%isort))then;deallocate(intwfnc%isort);nullify(intwfnc%isort);endif
    if(associated(intwfnv%isort))then;deallocate(intwfnv%isort);nullify(intwfnv%isort);endif
  endif
  if(allocated(indexq_fi))then;deallocate(indexq_fi);endif
  if(associated(xct%indexq_fi))then;deallocate(xct%indexq_fi);nullify(xct%indexq_fi);endif
  !----------------------------------------------------------------------------
  ! Finalize
  call write_memory_usage()
  call timing%stop(timing%total)
  call timing%print()
end program vmtxel
