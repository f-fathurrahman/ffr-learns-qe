!=================================================================================
!
! Routines:
!
! (1) epscopy() Last Modified 2014 (FHJ)
!
! Copy dielectric matrices from eps0mat(.h5)/epsmat(.h5) files to memory.
!
! This routine reads the eps0mat and epsmat files (units 10 and 11)
! and copies the relevant information about those dielectric matrices
! to memory (epsmpi). Processor 0 will actually read data and redistribute
! them to other processors. The eps0mat can have more than one q-point,
! which is useful when performing voronoi averages for the q->0 point.
!
!==================================================================================
module epscopy_m
  use global_m
  use misc_m
  use io_utils_m
  use scalapack_m
  use timing_m, only: timing => sigma_timing
  implicit none
  private
  public :: epscopy_init, epscopy
  type :: comm_buffer
    complex(DPC), allocatable :: msg(:,:)
    integer, allocatable :: col_global_indx(:)
    integer :: nrow, ncol
    integer :: proc
  end type
contains
!> Determines the number of q-points, the maximum rank of the dielectric matrix,
!! number of frequencies
subroutine epscopy_init(sig, neps)
  type (siginfo), intent(inout) :: sig
  integer, intent(out) :: neps !< Max. rank of dielectric matrix
  integer :: ngarbage1, nmtx, itrash, ig, igp, iq
  integer :: nFreq, nfreq_imag, freq_dep_flag
  real(DP), allocatable :: dFreqGrid(:)
  complex(DPC), allocatable :: dFreqBrd(:)
  logical, allocatable :: is_imag(:)
  logical :: file_exists
  logical :: subtrunc_flags(6)
  integer :: neig_sub, ierr
  real(DP) :: chi_eigenvalue_cutoff
 
  sig%subspace_q0 = .FALSE.
  sig%matrix_in_subspace_basis_q0 = .FALSE.
  sig%keep_full_eps_static_q0 = .FALSE.
  sig%subspace = .FALSE.
  sig%matrix_in_subspace_basis = .FALSE.
  sig%keep_full_eps_static = .FALSE.
  subtrunc_flags(:) = .FALSE.
  sig%neig_sub_max = 0
  if (sig%freq_dep>=0.and.sig%freq_dep<=3) then
      if (peinf%inode==0) then
        call open_file(unit=10,file='eps0mat',form='unformatted',status='old')
        read(10)
        ierr = 0
        read(10,IOSTAT=ierr) freq_dep_flag, sig%nFreq, &
                 sig%subspace_q0, sig%matrix_in_subspace_basis_q0, sig%keep_full_eps_static_q0
        IF(ierr .NE. 0) THEN
          ! probably this is due because of the absence of the subspace flags
          sig%subspace_q0 = .FALSE.
          sig%matrix_in_subspace_basis_q0 = .FALSE.
          sig%keep_full_eps_static_q0 = .FALSE.
        END IF
! make sure that in the case of subspace MPI and scalapack are linked
        IF(sig%subspace_q0) THEN
          call die('Subspace method only works with MPI and SCALAPACK')
        END IF
        ! Consistency check, before continuing
        if (freq_dep_flag==2.and.sig%freq_dep/=2) &
          call die('eps0mat is frequency-dependent, but this Sigma calculation is not.')
        if (freq_dep_flag==0.and.sig%freq_dep==2) then
          call die('This Sigma calculation is frequency-dependent, but eps0mat is not.')
        endif
        if (sig%freq_dep/=2 .and. sig%subspace_q0) then
          call die('Subspace approximation implemented only for freq_dep = 2')
        end if
        read(10)
        if (sig%freq_dep==2.or.sig%freq_dep==3) then
          allocate(dFreqGrid (sig%nFreq))
          allocate(dFreqBrd (sig%nFreq))
          allocate(is_imag (sig%nFreq))
          read(10) dFreqGrid, dFreqBrd
          is_imag = (dabs(dFreqGrid)<TOL_ZERO)
          is_imag(1) = .false.
          sig%nfreq_imag = count(is_imag)
          if(allocated(dFreqGrid))then;deallocate(dFreqGrid);endif
          if(allocated(dFreqBrd))then;deallocate(dFreqBrd);endif
          if(allocated(is_imag))then;deallocate(is_imag);endif
        else
          read(10)
        endif
        read(10)
        read(10)
        read(10)
        read(10) sig%nq0
        read(10)
        !XXXXXXXX
        IF(sig%subspace_q0 .AND. sig%matrix_in_subspace_basis_q0) THEN
          do iq=1,sig%nq0
            read(10) itrash, neps
            read(10) ! ekin
            read(10) ! q
            read(10) ! vcoul
            ! full epsilon omega zero
            IF(sig%keep_full_eps_static_q0) THEN
              do igp=1, neps
                read(10) ! For epsRDyn only
              end do
            END IF
            ! eigenvectors
            read(10) neig_sub
            IF(sig%neig_sub_max < neig_sub) sig%neig_sub_max = neig_sub
            do igp=1, neps
              read(10) ! For epsRDyn only
            end do
            ! subspace matrices
            do igp=1, neig_sub
              do ig=1, neig_sub
                read(10)
              end do
            end do
          end do ! iq
        ELSE
          read(10) itrash, neps
        END IF
        !XXXXXXXX
        call close_file(10)
        call open_file(unit=10,file='eps0mat',form='unformatted',status='old')
        inquire(file="epsmat", exist=file_exists)
        if (file_exists) then
          sig%igamma = 0
        else
          sig%igamma = 1
        endif
        if(sig%igamma/=0) then ! Gamma-only calculation
          sig%nq1 = 0
        else
          call open_file(unit=11,file='epsmat',form='unformatted',status='old')
          read(11)
          ierr = 0
          read(11,IOSTAT=ierr) ngarbage1, nFreq, &
                   sig%subspace, sig%matrix_in_subspace_basis, sig%keep_full_eps_static
          IF(ierr .NE. 0) THEN
            sig%subspace = .FALSE.
            sig%matrix_in_subspace_basis = .FALSE.
            sig%keep_full_eps_static = .FALSE.
          END IF
! make sure that in the case of subspace MPI and scalapack are linked
        IF(sig%subspace) THEN
          call die('Subspace method only works with MPI and SCALAPACK')
        END IF
          if (nFreq/=sig%nFreq) then
            call die('nFreq mismatch between eps0mat and epsmat')
          endif
          if (sig%freq_dep/=2 .and. sig%subspace) then
            call die('Subspace approximation implemented only for freq_dep = 2')
          end if
          read(11)
          if (sig%freq_dep==2.or.sig%freq_dep==3) then
            allocate(dFreqGrid (sig%nFreq))
            allocate(dFreqBrd (sig%nFreq))
            allocate(is_imag (sig%nFreq))
            read(11) dFreqGrid, dFreqBrd
            is_imag = (dabs(dFreqGrid)<TOL_ZERO)
            is_imag(1) = .false.
            if (sig%nfreq_imag/=count(is_imag)) then
              call die('nfreq_imag mismatch between eps0mat and epsmat')
            endif
            if(allocated(dFreqGrid))then;deallocate(dFreqGrid);endif
            if(allocated(dFreqBrd))then;deallocate(dFreqBrd);endif
            if(allocated(is_imag))then;deallocate(is_imag);endif
          else
            read(11)
          endif
          read(11)
          read(11)
          read(11)
          read(11) sig%nq1
          read(11)
          do iq=1,sig%nq1
            read(11) itrash, nmtx
            read(11)
            read(11)
            IF(sig%subspace .AND. sig%matrix_in_subspace_basis) THEN
              read(11) ! vcoul
              ! full epsilon omega zero
              IF(sig%keep_full_eps_static) THEN
                do igp=1, nmtx
                  read(11) ! For epsRDyn only
                end do
              END IF
              ! eigenvectors
              read(11) neig_sub
              do igp=1, nmtx
                read(11) ! For epsRDyn only
              end do
              ! subspace matrices
              do igp=1, neig_sub
                do ig=1, neig_sub
                  read(11)
                end do
              end do
              ! write(*,*) neig_sub, nmtx
              if(sig%neig_sub_max < neig_sub) sig%neig_sub_max=neig_sub
            ELSE
              if (sig%freq_dep==0.or.sig%freq_dep==1) then
                do ig=1,nmtx
                  read(11)
                enddo
              endif
              if (sig%freq_dep==2.or.sig%freq_dep==3) then
                do igp=1,nmtx
                  do ig=1,nmtx
                    read(11) ! For epsRDyn
                  enddo
                enddo
              endif
            END IF
            if(neps<nmtx) neps = nmtx
          enddo
          call close_file(11)
          call open_file(unit=11,file='epsmat',form='unformatted',status='old')
        endif
      endif ! peinf%inode==0
    if (sig%freq_dep==2.or.sig%freq_dep==3) then
      allocate(sig%dFreqGrid (sig%nFreq))
      allocate(sig%dFreqBrd (sig%nFreq))
    endif
  else ! sig%freq_dep>=0.and.sig_freq_dep<=2
    neps = 0
    sig%nFreq = 0
  endif
  if (sig%nq0/=1.and..not.sig%subsample) call die('epscopy_init: nq0/=1')
  if (peinf%inode==0) then
    write(6,'(/1x,a,i0)') 'Total number of frequencies in the dielectric matrix: ', sig%nfreq
    write(6,'(1x,a,i0/)') 'Number of imag. frequencies in the dielectric matrix: ', sig%nfreq_imag
    !XXXXX
    if(sig%subspace_q0) then
      write(6,'(1x,a)') 'Eps0mat calculated with subspace truncation method'
      IF(sig%matrix_in_subspace_basis_q0) THEN
        write(6,'(1x,a)') 'Eps0mat matrices written within the subspace of omega 0'
      END IF
      IF(sig%keep_full_eps_static_q0) THEN
        write(6,'(1x,a)') 'Eps0mat retains the full EpsInv for omega=0'
      END IF
      write(6,*)
    end if
    if(sig%subspace .AND. sig%igamma==0) then
      write(6,'(1x,a)') 'Epsmat calculated with subspace truncation method'
      IF(sig%matrix_in_subspace_basis) THEN
        write(6,'(1x,a)') 'Epsmat matrices written within the subspace of omega 0'
      END IF
      IF(sig%keep_full_eps_static) THEN
        write(6,'(1x,a)') 'Epsmat retains the full EpsInv for omega=0'
      END IF
      write(6,*)
    end if
    !XXXXX
  endif
  ! check if we can perform the calc in the subspace basis
  if(sig%do_sigma_subspace) then
    if(sig%igamma==0) then
      ! k-point, we need both
      if(sig%matrix_in_subspace_basis_q0 .AND. sig%matrix_in_subspace_basis) then
        if (peinf%inode==0) write(6,'(1x,a)') 'DIRECT evalution of sigma in subspace basis'
      else
        if (peinf%inode==0) write(6,'(1x,a)') 'Epsilon matrices not in subspace basis!'
        if (peinf%inode==0) write(6,'(1x,a)') 'Direct evalution of sigma in subspace basis turned off'
        sig%do_sigma_subspace = .false.
      end if
    else
      if(sig%matrix_in_subspace_basis_q0) then
        if (peinf%inode==0) write(6,'(1x,a)') 'DIRECT evalution of sigma in subspace basis'
      else
        if (peinf%inode==0) write(6,'(1x,a)') 'Epsilon matrix not in subspace basis!'
        if (peinf%inode==0) write(6,'(1x,a)') 'Direct evalution of sigma in subspace basis turned off'
        sig%do_sigma_subspace = .false.
      end if
    end if
  end if
 
end subroutine epscopy_init
subroutine epscopy(crys,gvec,sig,neps,epsmpi,epshead,iunit_eps,fne)
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (siginfo), intent(inout) :: sig
  integer, intent(in) :: neps !< number of G-vectors up to epsilon cutoff defined in sigma.inp
  type (epsmpiinfo), intent(inout) :: epsmpi
  real(DP), intent(out) :: epshead !< for full frequency, this is the retarded static head.
  integer, intent(in) :: iunit_eps
  character*20, intent(in) :: fne
  real(DP), allocatable :: eps(:)
  complex(DPC), allocatable :: epsRDyn(:,:)
  complex(DPC), allocatable :: epsADyn(:,:)
  logical :: is_static
 
  epshead = 0.0d0
  is_static = sig%freq_dep/=2 .AND. sig%freq_dep/=3
  ! FHJ: Initialize buffers
  if (is_static) then
    allocate(epsmpi%eps (neps,epsmpi%ngpown,sig%nq))
  else
    if(sig%do_sigma_subspace) then
      ! if doing calcultion within subspace allocate stuff
      allocate(sig%epssub%eps_sub (sig%neig_sub_max, sig%epssub%Nbas_own, sig%nFreq, sig%nq))
      !XXXX allocate(sig%epssub%eigenvec_sub (sig%epssub%ngpown_sub, sig%neig_sub_max, sig%nq))
      allocate(sig%epssub%eigenvec_sub (neps, sig%epssub%Nbas_own, sig%nq))
      allocate(sig%epssub%eps_wings_rows (neps, sig%nFreq, sig%nq))
      allocate(sig%epssub%eps_wings_cols (neps, sig%nFreq, sig%nq))
      allocate(sig%epssub%eps_wings_correction_rows (neps, sig%nFreq))
      allocate(sig%epssub%eps_wings_correction_cols (neps, sig%nFreq))
      sig%epssub%eps_sub = (0.0d0,0.0d0)
      sig%epssub%eigenvec_sub = (0.0d0,0.0d0)
      sig%epssub%eps_wings_rows = (0.0d0,0.0d0)
      sig%epssub%eps_wings_cols = (0.0d0,0.0d0)
      sig%epssub%eps_wings_correction_rows = (0.0d0,0.0d0)
      sig%epssub%eps_wings_correction_cols = (0.0d0,0.0d0)
      allocate(sig%epssub%eps_sub_info (3, 2, sig%nq))
      allocate(sig%epssub%wing_pos (sig%nq))
      sig%epssub%wing_pos = 0
      allocate(sig%epssub%vcoul_sub (neps, sig%nq))
      sig%epssub%vcoul_sub = 0.0d0
      !XXX deallocate and reallocate epsmpi%epsR with a single frequency
      !XXX if(associated(epsmpi%epsR))then;deallocate(epsmpi%epsR);nullify(epsmpi%epsR);endif
      allocate(epsmpi%epsR (neps,epsmpi%ngpown,1,sig%nq))
    else
      ! standard case
      allocate(epsmpi%epsR (neps,epsmpi%ngpown,sig%nFreq,sig%nq))
      if (sig%need_advanced) then
        allocate(epsmpi%epsA (neps,epsmpi%ngpown,sig%nFreq,sig%nq))
      endif
    end if
  endif
  ! FHJ: Temp buffers.
  if (is_static) then
    allocate(eps (neps))
  else
    allocate(epsRDyn (neps,sig%nFreq))
    if (sig%need_advanced) then
      allocate(epsADyn (neps,sig%nFreq))
    endif
  endif
  !------------------------------------------------
  ! FHJ: Read dielectric matrices from eps0mat(.h5)
  !------------------------------------------------
  call read_epsmat(.true.)
  if (dble(epshead)<1d-3 .and. sig%iscreen==SCREEN_SEMICOND .and. peinf%inode==0) then
    write(0,'(a)') 'WARNING: You are using semiconductor screening, but the'
    write(0,'(a)') 'head of epsilon inverse is very small and seems metallic.'
  endif
  !------------------------------------------------
  ! FHJ: Read dielectric matrices from epsmat(.h5)
  !------------------------------------------------
  if (sig%igamma==0) call read_epsmat(.false.)
  epsmpi%qk(:,:) = sig%qpt(:,:)
  ! FHJ: Free temp. buffers
  if (is_static) then
    if(allocated(eps))then;deallocate(eps);endif
  else
    if(allocated(epsRDyn))then;deallocate(epsRDyn);endif
    if (sig%need_advanced) then
      if(allocated(epsADyn))then;deallocate(epsADyn);endif
    endif
  endif
 
  return
contains
  subroutine read_epsmat(is_q0)
    logical, intent(in) :: is_q0
    character(len=16) :: fname
    integer :: qoffset, nq_file, iunit
    integer :: freq_dep_flag, nFreq, nfreq_imag, ng_old, ngq, nmtx
    integer :: ii, ig, igp, nq_tmp, iq, qgrid(3), gg(3), iout, iw, dest, igp_loc
    real(DP) :: ecuts, qk(3), ekin
    real(DP), allocatable :: ekold(:)
    integer, allocatable :: isrtq(:), isrtqi(:), isrtold(:), gvecs_old(:,:)
    type(progress_info) :: prog_info !< a user-friendly progress report
    character(len=6) :: sname
    character(len=11) :: sdate
    logical :: matrix_in_subspace_basis, keep_full_eps_static
    real(DP), allocatable :: vcoul(:)
   
    if (is_q0) then
      fname = 'eps0mat'
      qoffset = 0
      nq_file = sig%nq0
      iunit = 10
    else
      fname = 'epsmat'
      qoffset = sig%nq0
      nq_file = sig%nq1
      iunit = 11
    endif
    if (sig%use_hdf5) fname = TRUNC(fname) // '.h5'
    if (peinf%inode==0) then
      allocate(ekold (gvec%ng))
      allocate(isrtold (gvec%ng))
      allocate(isrtq (gvec%ng))
      allocate(isrtqi (gvec%ng))
    endif
    matrix_in_subspace_basis = .FALSE.
    keep_full_eps_static = .FALSE.
    if (is_q0) then
      if(sig%subspace_q0) then
        matrix_in_subspace_basis = sig%matrix_in_subspace_basis_q0
        keep_full_eps_static = sig%keep_full_eps_static_q0
      end if
    else
      if(sig%subspace) then
        matrix_in_subspace_basis = sig%matrix_in_subspace_basis
        keep_full_eps_static = sig%keep_full_eps_static
      end if
    end if
!------------------------------------------------------------------------------
! Read q-grid, q-points, G-vectors and freq. grid
!------------------------------------------------------------------------------
    if (peinf%inode==0) then
        read(iunit) sname, sdate
        read(iunit) freq_dep_flag, nFreq
        read(iunit) qgrid(1:3)
        if (.not.is_static.and.is_q0) then
          read(iunit) sig%dFreqGrid, sig%dFreqBrd
        else
          read(iunit)
        endif
        read(iunit)
        read(iunit)
        read(iunit) ecuts
        read(iunit) nq_tmp, ((sig%qpt(ii,qoffset+iq), ii=1,3), iq=1,nq_tmp)
        read(iunit) ng_old
        backspace(iunit)
        allocate(gvecs_old (3, ng_old))
        read(iunit) ng_old, ((gvecs_old(ii,iq), ii=1,3), iq=1,ng_old)
      if (nq_tmp/=nq_file) call die('nq mismatch for '//TRUNC(fname))
      ! FHJ: only substitute sig%qgrid (used in voronoi averages) if sig%qgrid
      ! is empty and this is epsmat(.h5), if this file is available
      if (all(sig%qgrid(1:3)==0).and.((.not.is_q0).or.(sig%igamma==1))) then
        sig%qgrid(:) = qgrid(:)
      endif
    endif ! peinf%inode==0
!------------------------------------------------------------------------------
! Bcast q-points and allocate/bcast full-frequency and epsilon buffers
!------------------------------------------------------------------------------
    if (is_q0) then
      sig%q0vec(:) = sig%qpt(:,1)
      ! We need to manually set the q0 point in sig%qpt to zero
      if (.not.sig%subsample) sig%qpt(:,:sig%nq0) = 0d0
      if (.not.is_static) then
      endif
    endif !is_q0
!------------------------------------------------------------------------------
! Loop over all q-points, map G-vectors and read/store dielectric matrices
!------------------------------------------------------------------------------
    call progress_init(prog_info, 'reading '//TRUNC(fname), 'q-point', nq_file)
    do iq=1,nq_file
      call progress_step(prog_info, iq)
      call logit('Storing eps to memory')
      ! FHJ: Map old isrt (from epsilon gspace) to new isrt (to gvec gspace)
      !--------------------------------------------------------=------------
      if (peinf%inode==0) then
          read(iunit) ngq, nmtx, (isrtold(ig),ii,ig=1,ngq)
          read(iunit) (ekold(ig),ig=1,ngq)
          read(iunit) (qk(ii),ii=1,3)
          ! in case of subspace read the coulomb potential
          if(matrix_in_subspace_basis) then
            allocate(vcoul (nmtx))
            read(iunit) (vcoul(ig),ig=1,nmtx)
          end if
        isrtq = 0
        isrtqi = 0
        qk = sig%qpt(:,iq+qoffset)
        if (is_q0) qk(:) = 0d0
        ! FHJ: the following check seems to be important when the WFN file used in the
        ! epsilon calculation has a cutoff larger than that of the WFN file used here.
        ngq = min(ngq, gvec%ng)
        do ig=1,ngq
          ! FHJ: note that size(ekold) is gvec%ng
          if (isrtold(ig) <= gvec%ng) then
            if (ekold(isrtold(ig))<=sig%ecutb) then
              gg(:) = gvecs_old(:, isrtold(ig))
              call findvector(iout, gg, gvec)
              isrtq(ig) = iout
              isrtqi(iout) = ig
              if (ig==1) then
                ! just check the first so we do not waste time
                ekin = DOT_PRODUCT(gvec%components(:,iout)+qk(:), MATMUL(crys%bdot, gvec%components(:,iout)+qk(:)))
                if (abs(ekin-ekold(isrtold(ig))) > TOL_Zero) then
                  write(0,*) 'iq = ', iq, ' ig = ',ig, ' qk = ', qk
                  write(0,*) 'epsmat: ekold(isrtold(i)) = ', ekold(isrtold(ig)), ' ekin = ', ekin
                  call die("Incorrect kinetic energies in epsmat.")
                endif
              endif
            endif ! if (ekold(isrtold(i))<=sig%ecutb)
          endif
        enddo ! do ig=1,ngq
      endif ! if (peinf%inode==0)
      if (peinf%inode==0) then
        epsmpi%isrtq(:,iq+qoffset) = isrtq(:)
        epsmpi%isrtqi(:,iq+qoffset) = isrtqi(:)
      endif
      epsmpi%nmtx(iq+qoffset) = nmtx
      IF(matrix_in_subspace_basis) THEN
        ! broadcast coulomb potential
        IF(peinf%inode/=0) THEN
          allocate(vcoul (nmtx))
          vcoul = 0.0_DP
        END IF
        IF(sig%need_advanced) THEN
          write(6,'(1x,a)') &
          'Subspace truncation method only implemented for sig%freq_dep==2 and sig%freq_dep_method==2'
          call die('need_advanced is an invalid option for subspace + cd_integration_method.')
        END IF
        if (is_q0.and.iq==1.and.peinf%inode==0) then
          epshead = epsmpi%epsR(1,1,1,1)
        endif
        if(allocated(vcoul))then;deallocate(vcoul);endif
      ELSE ! subspace read
          ! FHJ: Read & store epsinv === Fortran binary ===
          !-------------------------------------------------
          do igp=1,nmtx
            if (is_static) then
              if (peinf%inode==0) then
                call timing%start(timing%epscopy_io)
                read(iunit) (eps(ig),ig=1,nmtx)
                call timing%stop(timing%epscopy_io)
                if (is_q0.and.iq==1.and.igp==1) epshead = eps(1)
                call timing%start(timing%epscopy_comm)
                call timing%stop(timing%epscopy_comm)
              endif
              call timing%start(timing%epscopy_comm)
              call timing%stop(timing%epscopy_comm)
            else ! is_static
              if (peinf%inode==0) then
                do ig=1,nmtx
                  call timing%start(timing%epscopy_io)
                  read(iunit) (epsRDyn(ig,iw),iw=1,sig%nFreq) ! Retarded part
                  call timing%stop(timing%epscopy_io)
                  if (is_q0.and.iq==1.and.igp==1) epshead = epsRDyn(1,1)
                  call timing%start(timing%epscopy_comm)
                  call timing%stop(timing%epscopy_comm)
                enddo
              endif
              call timing%start(timing%epscopy_comm)
              call timing%stop(timing%epscopy_comm)
            endif ! is_static
            dest = INDXG2P(igp, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
            igp_loc = INDXG2L(igp, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
            call timing%start(timing%epscopy_comm)
            if (peinf%pool_rank==dest) then
              if (is_static) then
                epsmpi%eps(:,igp_loc,iq+qoffset) = eps(:)
              else
                epsmpi%epsR(:,igp_loc,:,iq+qoffset) = epsRDyn(:,:)
                if (sig%need_advanced) then
                  epsmpi%epsA(:,igp_loc,:,iq+qoffset) = epsADyn(:,:)
                endif
              endif ! is_static
            endif ! if (peinf%pool_rank==dest)
            call timing%stop(timing%epscopy_comm)
          enddo ! do ig=1,nmtx
      END IF ! subspace read
    enddo ! iq
    call progress_free(prog_info)
!------------------------------------------------------------------------------
! Done!
!------------------------------------------------------------------------------
    if (peinf%inode==0) then
        call close_file(iunit)
      if(allocated(ekold))then;deallocate(ekold);endif
      if(allocated(isrtold))then;deallocate(isrtold);endif
      if(allocated(isrtq))then;deallocate(isrtq);endif
      if(allocated(isrtqi))then;deallocate(isrtqi);endif
    endif
   
  end subroutine read_epsmat
end subroutine epscopy
end module epscopy_m
