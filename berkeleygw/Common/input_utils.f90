!===============================================================================
!
! Modules:
!
! input_utils_m Originally By DAS
!
! Several routines useful for analysis after input of wavefunctions.
!
!===============================================================================
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
module input_utils_m
  use global_m
  use blas_m
  implicit none
  private
  public :: &
    assess_degeneracies, &
    find_efermi, &
    calc_qtot, &
    gvec_index, &
    kinetic_energies, &
    check_trunc_kpts, &
    eps_setup_sizes
contains
  !> calculates the degeneracy of each state in the kp object
  !> if sig is provided, dies if Sigma diagonals cannot be correctly calculated
  subroutine assess_degeneracies(kp, kp_el_extra, nband, efermi, tolerance, sig, ncore_excl)
    type(kpoints), intent(inout) :: kp
    real(DP), intent(in) :: kp_el_extra(:, :) !< (kp%nrk, kp%nspin) !< one more band
    integer, intent(in) :: nband !< number of states
    real(DP), intent(in) :: efermi
    real(DP), intent(in) :: tolerance !< energy tolerance
    type(siginfo), optional, intent(in) :: sig
    integer, optional, intent(in) :: ncore_excl !< number of core states excluded
    integer :: ii, jj, ik, is, band_offset
    real(DP) :: energy_compare
    character :: tmpstr*100,tmpstr1*16,tmpstr2*16
    logical :: written_header
   
    band_offset=0
    if (present(ncore_excl)) then
      band_offset = ncore_excl
    endif
    allocate(kp%degeneracy (nband, kp%nrk, kp%nspin))
    kp%degeneracy(:,:,:) = 1
    written_header = .false.
    !FHJ: do we really need two nested loop over nband??
    do ii = 1, nband
      do jj = 1, nband + 1
        if(jj == ii) cycle ! self-degeneracy
        do ik = 1, kp%nrk
          do is = 1, kp%nspin
            if(jj > nband) then
              energy_compare = kp_el_extra(ik, is)
            else
              energy_compare = kp%el(jj+band_offset, ik, is)
            endif
             ! DVF : we add band_offset to ii here because the energies in kp%el
             ! are referenced to the case including the core states since they are
             ! setup during wfn i/o.
            if(abs(kp%el(ii+band_offset, ik, is) - energy_compare) .lt. tolerance) then
              kp%degeneracy(ii, ik, is) = kp%degeneracy(ii, ik, is) + 1
              if(abs(kp%el(ii+band_offset, ik, is) - efermi / ryd) < tolerance) then
                if(.not. written_header) then
                  write(0,'(a)') "WARNING: Degeneracies at Fermi level:"
                  written_header = .true.
                endif
                write(0,'(a,i6,a,i6,a,i6)') 'k ', ik, ', spin ', is, ', band ', ii
              endif
              if(present(sig)) then
              ! There is a problem only if: we are trying to use q-symmetry,
              ! we are calculating Sigma for this k-point and band, and the
              ! other degenerate states are not being included.
              ! DVF : we add band_offset to ii here because the matrix elements
              ! are referenced to the case including the core states.
                if(sig%qgridsym .and. any(sig%indkn(1:sig%nkn) == ik) .and. &
                  any(sig%diag(1:sig%ndiag) .eq. ii+band_offset) .and. all(sig%diag(1:sig%ndiag) .ne. jj+band_offset)) then
                  write(tmpstr1,660) ii
                  write(tmpstr2,660) jj
660 format(i16)
                  write(tmpstr,'(a,a,a,a)') "Cannot correctly calculate Sigma for band ", TRUNC(tmpstr1), &
                    " without its degenerate partner ", TRUNC(tmpstr2)
                  write(0,'(a)') "Run degeneracy_check.x for allowable numbers or set no_symmetries_q_grid."
                  call die(tmpstr)
                endif
              endif
            endif
          enddo ! is
        enddo ! ik
      enddo ! ii
    enddo ! jj
    if(written_header) write(0,*)
   
    return
  end subroutine assess_degeneracies
  !---------------------------------------------------------------------------------------------------
  !> only fine unshifted grids use should_search = .true., should_update = .true.
  !! coarse grids have should_search, should update true or false depending on eqp corrections vs fine grid
  !! shifted grids use should_search = .false. (since may have only valence bands), should_update = .false.
  subroutine find_efermi(rfermi, efermi, efermi_input, kp, nband, minband, label, should_search, should_update, &
    write7, dont_die_consistency)
    logical, intent(in) :: rfermi !< relative or absolute Fermi level
    real(DP), intent(inout) :: efermi !< The computed Fermi energy
    real(DP), intent(in) :: efermi_input !< Will be used to set the Fermi energy (rfermi=.false.)
                                          !! or shift the Fermi energy (rfermi=.true.)
    type(kpoints), intent(inout) :: kp !< Electronic structure from which the Fermi level is obtained
    integer, intent(in) :: nband !< total number of bands to consider
    integer, intent(in) :: minband !< lowest band to consider
    character(len=*), intent(in) :: label !< name to be used in output
    logical, intent(in) :: should_search !< if true, efermi is calculated
    logical, intent(in) :: should_update !< if true, efermi is reset
    logical, intent(in) :: write7 !< write to unit 7 as well as unit 6
    !> If .true., don`t die if the consistency of ifmin/ifmax are wrong.
    !! This is only ok in inteqp with unrestricted_transformation. Default=.false.
    logical, intent(in), optional :: dont_die_consistency
    integer :: ik, is, ib, consistency_err
    real(DP) :: emiddle, efermi_temp
    logical :: should_reset_fermi, should_warn_occ, first_occ_warning, die_consistency
    real(DP), allocatable :: vbm(:,:), cbm(:,:)
   
    die_consistency = .true.
    if (present(dont_die_consistency)) then
      die_consistency = .not. dont_die_consistency
    endif
    should_reset_fermi = .false.
    if(should_update .and. .not. should_search) then
      call die("BUG: cannot call find_efermi with should_update but not should_search")
    endif
    if(nband < 1 .or. nband > kp%mnband) then
      call die("find_efermi: nband out of bounds")
    endif
    if(peinf%inode == 0) then
      write(6,900) TRUNC(label), maxval(kp%ifmax(1:kp%nrk, 1:kp%nspin))
      if(write7) write(7,900) TRUNC(label), maxval(kp%ifmax(1:kp%nrk, 1:kp%nspin))
900 format(1x,'Highest occupied band (', a, ') = ',i0)
    endif
    if(all(kp%ifmax(1:kp%nrk, 1:kp%nspin) <= 0)) then
      call die("All k-points have no occupied bands.")
    endif
    if(should_search) then
      if(any(kp%ifmax(1:kp%nrk, 1:kp%nspin) < nband)) then
        allocate(vbm (1:kp%nrk, 1:kp%nspin))
        allocate(cbm (1:kp%nrk, 1:kp%nspin))
        do is=1,kp%nspin
          do ik=1,kp%nrk
            vbm(ik, is) = maxval(kp%el(minband:kp%ifmax(ik, is), ik, is))
            cbm(ik, is) = minval(kp%el(kp%ifmax(ik, is)+1:nband, ik, is))
          enddo
        enddo
        emiddle = (maxval(vbm) + minval(cbm)) / 2.0d0
        if(peinf%inode == 0) then
          write(6,901) 'Valence max   ', TRUNC(label), maxval(vbm)*ryd
          write(6,901) 'Conduction min', TRUNC(label), minval(cbm)*ryd
          write(6,901) 'Middle energy ', TRUNC(label), emiddle*ryd
          if(write7) then
            write(7,901) 'Valence max   ', TRUNC(label), maxval(vbm)*ryd
            write(7,901) 'Conduction min', TRUNC(label), minval(cbm)*ryd
            write(7,901) 'Middle energy ', TRUNC(label), emiddle*ryd
          endif
901 format(1x,a,' (', a, ') = ',f0.6,' eV')
        endif
        ! check for consistency: i.e. no valence band or conduction band is on the wrong side of the Fermi energy
        consistency_err = 0 !0 = no error
        do is = 1, kp%nspin
          do ik = 1, kp%nrk
            if(kp%ifmax(ik, is) > 0) then
              if(any(kp%el(minband:kp%ifmax(ik, is), ik, is) > emiddle + TOL_Zero)) then
                consistency_err = ior(consistency_err, 1)
              endif
            endif
            if(kp%ifmax(ik, is) + 1 <= nband) then
              if(any(kp%el(kp%ifmax(ik, is) + 1:nband, ik, is) < emiddle - TOL_Zero)) then
                consistency_err = ior(consistency_err, 2)
              endif
            endif
          enddo
        enddo
        if (consistency_err>0 .and. die_consistency) then
          if (peinf%inode==0) then
            write(0,*)
            write(0,*) 'WFN ifmin/ifmax fields are inconsistent:'
            if (iand(consistency_err, 1)==1) then
              write(0,*) ' - there is a valence state above the middle energy'
            endif
            if (iand(consistency_err, 2)==2) then
              write(0,*) ' - there is a conduction state below the middle energy'
            endif
            write(0,*) 'Possible causes are:'
            write(0,*) '(1) Your k-point sampling is too coarse and cannot resolve the Fermi energy.'
            write(0,*) ' Try to carefully inspect your mean-field energies, and consider using a finer k-grid.'
            write(0,*) '(2) You are using eqp.dat and the QP corrections change the character of some states'
            write(0,*) ' from valence<->conduction. In this case, you should use another mean-field theory'
            write(0,*) ' that gives the same ground state as your GW calculation.'
            write(0,*) '(3) You are running inteqp, but you are either shifting the Fermi energy or using '
            write(0,*) ' restricted transformation.'
          endif
          call die("WFN ifmin/ifmax fields are inconsistent", only_root_writes = .true.)
        else if (consistency_err>0 .and. peinf%inode==0) then
          write(0,*)
          write(0,*) 'WARNING: ifmin/ifmax fields are inconsistent. Beware that the reported Fermi '
          write(0,*) ' and middle energy are probably wrong, as well as the order of the states.'
          write(0,*)
        endif
        ! adjust by level set in input file if appropriate
        if(rfermi) then
          efermi_temp = emiddle * ryd + efermi_input
        else
          efermi_temp = efermi_input
        endif
        if(peinf%inode == 0) then
          write(6,902) TRUNC(label), efermi_temp
          if(write7) write(7,902) TRUNC(label), efermi_temp
902 format(1x,'Fermi  energy  (', a, ') = ',f0.6,' eV',/)
        endif
        if(should_update) then
          efermi = efermi_temp
        endif
        should_reset_fermi = (abs(efermi / ryd - emiddle) .gt. TOL_Small)
      else
        if(peinf%inode == 0) write(0,'(a)') &
          'WARNING: There are only valence bands present; cannot determine Fermi energy.'
        if ((rfermi.and.efermi_input/=0).or..not.rfermi) then
          call die('No conduction bands; cannot determine nor shift the Fermi energy.', &
            only_root_writes=.true.)
        endif
        ! pick a large arbitrary value to initialize and make sure we consider all the bands fully occupied
        if(should_update) efermi = INF
      endif
    endif ! should_search
    ! reset ifmax if Fermi level was moved by input file, or using Fermi level from other wfns
    ! if neither of these is true, we should have died above if any occs are inconsistent with Fermi level!
    ! FHJ: reset occupations if we are dealing with...
    ! 1) a regular WFN_fi file (should_search/update==.true.) and the FE changed; or
    ! 2) a WFNq_fi file (should_search/update==.false.), in any condition; or
    ! 3) BSE/input_co.f90 (should_search==.true.) with eqp_co (should_update=.true.)
    ! and the FE changed => equivalent to (1)
    if ( (should_update.and.should_reset_fermi).or.(.not.should_search) ) then
!FHJ: never warn about resetting occupations if we manually moved the FE
      should_warn_occ = (peinf%inode==0) .and. ((dabs(efermi_input)<TOL_SMALL).and.rfermi)
      first_occ_warning = .true.
      do is = 1, kp%nspin
        do ik = 1, kp%nrk
          if (kp%ifmin(ik, is) == 0 .and. kp%el(1, ik, is) .lt. efermi / ryd) then
            kp%ifmin(ik, is) = 1
          endif
          if(kp%ifmin(ik, is) == 0) cycle
          do ib = max(minband, kp%ifmin(ik, is)), nband - 1
            if (kp%el(ib, ik, is) .lt. efermi / ryd .and. kp%el(ib + 1, ik, is) .gt. efermi / ryd) then
              if(kp%ifmax(ik, is) /= ib) then
                if (should_warn_occ) then
                  if (first_occ_warning) then
                    write(0,'(3a)') ' WARNING: resetting occupations of ', TRUNC(label), ' according to Fermi level.'
                    write(0,'(a)') ' Make sure your shifted and unshifted WFNs are consistent. The following kpoints'
                    write(0,'(a)') '  had their occupations reset (spin, kpoint, original band, energy, reset band, energy):'
                    first_occ_warning = .false.
                  endif
                  write(0,'(3i9,f12.6,i9,f12.6)') is, ik, kp%ifmax(ik, is), &
                    kp%el(kp%ifmax(ik,is), ik, is)*ryd, ib, kp%el(ib, ik, is)*ryd
                endif
                kp%ifmax(ik, is) = ib
                kp%occ(ib+1:, ik, is) = 0.0
                kp%occ(:ib, ik, is) = 1.0
                cycle !FHJ: we can`t reset a kpt twice!
              endif
            endif
          enddo
        enddo
      enddo
      if (.not.first_occ_warning) write(0, *)
    endif !should_reset_fermi .or.
!FHJ: TODO - create input flag to output occupations.
! call write_occupations(kp)
   
    return
  end subroutine find_efermi
  !> Write LDA energies, shifted energies and occupations to file OCCUPATIONS.
  !!
  !! Useful to debug whether scissors operators/spline adjustment to band
  !! structure is correct.
  subroutine write_occupations(kp)
    type (kpoints), intent(in) :: kp
    integer :: is, ik, nb
   
    nb=8 !TODO create an input flag to set number of energies to output
    if (peinf%inode==0) then
      call open_file(unit=60, file="OCCUPATIONS", status='replace', form='formatted')
      write(60, '(a)') '# kx ky kz  spin  ifmin ifmax'
      write(60, '(a)') '#   LDA energies before shift (in eV)'
      write(60, '(a)') '#   LDA energies after scissors/spline shift (in eV)'
      write(60, '(a)') '#   occupations'
      do is=1, kp%nspin
        do ik = 1, kp%nrk
          write(60, '(3(F9.5,1x),4x,I1,2x,I3,1x,I3)') kp%rk(:, ik), is, kp%ifmin(ik, is), kp%ifmax(ik, is)
          write(60, '(2x, 8(F12.5,1x))') kp%elda(1:8, ik, is)
          write(60, '(2x, 8(F12.5,1x))') kp%el(1:8, ik, is)*ryd
          write(60, '(2x, 8(F12.5,1x))') kp%occ(1:8, ik, is)
        enddo
      enddo
      call close_file(60)
    endif
   
    return
  end subroutine write_occupations
  !---------------------------------------------------------------------------------------------------
  ! Calculate charge in the cell, from ifmax and from occupations. Do some checks,
  ! and write results to output.
  subroutine calc_qtot(kp, celvol, efermi, qtot, omega_plasma, write7)
    type (kpoints), intent(in) :: kp
    real(DP), intent(in) :: celvol
    real(DP), intent(in) :: efermi
    real(DP), intent(out) :: qtot
    real(DP), intent(out) :: omega_plasma
    logical, intent(in) :: write7
    integer :: is, ik
    real(DP) :: qkpt, qkpt_occ, qtot_occ
   
    qtot = 0.0d0
    qtot_occ = 0d0
    do is = 1, kp%nspin
      do ik = 1, kp%nrk
        qkpt_occ = sum(kp%occ(:, ik, is))
        if(kp%ifmax(ik, is) == 0) then
          if(kp%el(1, ik, is) < efermi / ryd + TOL_Degeneracy) then
            qkpt = 0.5d0
          else
            qkpt = 0.0d0
          endif
        else
          qkpt = kp%ifmax(ik, is) - kp%ifmin(ik, is)
          if(kp%el(kp%ifmax(ik, is), ik, is) > efermi / ryd - TOL_Degeneracy) then
            qkpt = qkpt + 0.5d0
          else
            qkpt = qkpt + 1.0d0
          endif
        endif
        qtot = qtot + qkpt * kp%w(ik)
        qtot_occ = qtot_occ + qkpt_occ * kp%w(ik)
      enddo
    enddo
    qtot = qtot * 2.0d0 / dble(kp%nspin*kp%nspinor)
    qtot_occ = qtot_occ * 2.0d0 / dble(kp%nspin*kp%nspinor)
    omega_plasma = 4.0d0 * sqrt(PI_D * qtot / celvol)
    if(peinf%inode == 0) then
      write(6,904) qtot
      if(write7) write(7,904) qtot
904 format(1x,'Number of electrons per unit cell (from ifmax) = ',f0.6)
      if(abs(qtot - nint(qtot)) > TOL_Small * 100) then
        write(0,'(a)') " WARNING: Fractional number of electrons per unit cell (from ifmax)."
      endif
      if(qtot < TOL_Small) then
        write(0,'(a)') " WARNING: No electrons per unit cell (from ifmax)!"
      endif
      write(6,906) qtot_occ
      if(write7) write(7,906) qtot_occ
906 format(1x,'Number of electrons per unit cell (from occupations) = ',f0.6)
      if(abs(qtot_occ - nint(qtot_occ)) > TOL_Small * 100) then
        write(0,'(a)') " WARNING: Fractional number of electrons per unit cell (from occupations)."
      endif
      if(qtot_occ < TOL_Small) then
        write(0,'(a)') " WARNING: No electrons per unit cell (from occupations)!"
      endif
      if(abs(qtot - qtot_occ) > TOL_Small * 100) then
        write(0,'(a)') " WARNING: Discrepancy between number of electron per unit cell from ifmax and from occupations."
      endif
      write(6,905) omega_plasma
      if(write7) write(7,905) omega_plasma
905 format(1x,'Plasma Frequency = ',f0.6,' Ry',/)
    endif
   
    return
  end subroutine calc_qtot
!> FOR HISTORICAL INTEREST: SIB describes original scheme, GSM describes new scheme
!! SIB: gvec%nFFTgridpts is product of 1+2*FFTgrid(i) over all directions i
!! gvec%index_vec(gvec%nFFTgridpts) is allocated.
!!
!! gvec%index_vec(:) is a table of g-vector addresses: put in address,
!! and get out g-vector index in gvec%components(1:3,index) so that
!! index = gvec%index_vec(address) points to gvec%components(:,index).
!!
!! What does this imply about the organization of gvec%components(:,:) ?
!! The address of gvec%components(:,i) is given by
!! address=((gx_i+gxmax)*(2*gymax+1)+gy_i+gymax)*(2*gzmax+1)+gz_i+gzmax+1
!! so that this means z-fastest, then y, then x. Also the g-vectors
!! live on a grid that is 2*gmax+1 in each axial direction.
!!
!! gsm: saving on memory: FFTgrid is FFT grid size, gvec%nFFTgridpts=FFTgrid(1)*FFTgrid(2)*FFTgrid(3),
!! address=((g(1)+FFTgrid(1)/2)*FFTgrid(2)+g(2)+FFTgrid(2)/2)*FFTgrid(3)+g(3)+FFTgrid(3)/2+1
  !---------------------------------------------------------------------------------------------------
  !> Compute index_vec indices relating G-vectors in reduced coordinates to positions in the FFT grid
  subroutine gvec_index(gvec)
    type(gspace), intent(inout) :: gvec
    integer :: ig, iadd
   
    gvec%nFFTgridpts = product(gvec%FFTgrid(1:3))
    allocate(gvec%index_vec (gvec%nFFTgridpts))
    gvec%index_vec(:)=0
    do ig = 1, gvec%ng
      ! if a mean-field code does not use the appropriate convention, this could happen.
      if(any(2 * gvec%components(1:3, ig) >= gvec%FFTgrid(1:3) .or. 2 * gvec%components(1:3, ig) < -gvec%FFTgrid(1:3))) &
        call die("gvectors must be in the interval [-FFTgrid/2, FFTgrid/2)")
      iadd = ((gvec%components(1,ig)+gvec%FFTgrid(1)/2)*gvec%FFTgrid(2)+gvec%components(2,ig)+ &
        gvec%FFTgrid(2)/2)*gvec%FFTgrid(3)+gvec%components(3,ig)+gvec%FFTgrid(3)/2+1
      gvec%index_vec(iadd) = ig
    enddo
   
    return
  end subroutine gvec_index
  !-----------------------------------------------------------------
  !> This routine calculates the kinetic energies |G+q|^2 or |G|^2 for all
  !! the G-vectors in gvec%components using the reciprocal metric bdot:
  !! ekin(ig) = \sum_{m,n} G(m,ig) B(m,n) G(n,ig)
  !! We perform the sum by first performing the Cholesky decomposition of B,
  !! B := U^T U. Then, we write V := U G and write ekin = V^T V
  !! Using Cholesky decomposition has the same flop count as using dgemms, but
  !! it`s easier for the compiler to vectorize.
  !!
  !! \param gvec gspace structure that contains all the gvectors gvec%components
  !! \param bdot reciprocal metric
  !! \param ekin array holding the output kinetic energies
  !! \param qvec use this to compute |q+G|^2 instead of |G|^2
  !! \param qing_q2m optional parameter which specifies the order to pick the
  !! master G-vectors at q-point qvec (this is like the array "isort")
  subroutine kinetic_energies(gvec, bdot, ekin, qvec, gind_q2m)
    type(gspace), intent(in) :: gvec
    real(DP), intent(in) :: bdot(3, 3)
    real(DP), intent(out) :: ekin(:) !< (gvec%ng)
    real(DP), optional, intent(in) :: qvec(3)
    integer, optional, intent(in) :: gind_q2m(:)
    integer :: ig, info
    real(DP) :: q0(3), qkv(3), vmid(3), U(3,3) ! FHJ: stack allocation is faster!
   
    q0 = 0d0
    if (present(qvec)) then
      q0 = qvec
    endif
    U(1:3, 1:3) = bdot(1:3, 1:3)
    ! FHJ: Cholesky decomposition of the metric: bdot = U^T U
    call dpotrf('U', 3, U, 3, info)
    if (present(gind_q2m)) then
      ! FHJ: OpenMP decreases the performance of this loop by 10x.
      ! Maybe that`s just the OpenMP overhead?
      ! FHJ2: OpenMP will break threading, unless we use taskloop.
      !!disabled PARALLEL DO PRIVATE(qkv,vmid,ig)
      do ig = 1, size(gind_q2m)
        qkv = q0 + gvec%components(:,gind_q2m(ig))
        vmid(1) = U(1,1)*qkv(1) + U(1,2)*qkv(2) + U(1,3)*qkv(3)
        vmid(2) = U(2,2)*qkv(2) + U(2,3)*qkv(3)
        vmid(3) = U(3,3)*qkv(3)
        ekin(ig) = vmid(1)**2 + vmid(2)**2 + vmid(3)**2
      enddo
      !!disabled END PARALLEL DO
    else
      ! FHJ: OpenMP decreases the performance of this loop by 10x.
      ! Maybe that`s just the OpenMP overhead?
      ! FHJ2: OpenMP will break threading, unless we use taskloop.
      !!disabled PARALLEL DO PRIVATE(qkv,vmid,ig)
      do ig = 1, gvec%ng
        qkv = q0 + gvec%components(:,ig)
        vmid(1) = U(1,1)*qkv(1) + U(1,2)*qkv(2) + U(1,3)*qkv(3)
        vmid(2) = U(2,2)*qkv(2) + U(2,3)*qkv(3)
        vmid(3) = U(3,3)*qkv(3)
        ekin(ig) = vmid(1)**2 + vmid(2)**2 + vmid(3)**2
      enddo
      !!disabled END PARALLEL DO
    endif
   
  end subroutine kinetic_energies
  !-----------------------------------------------------------------
  !> Write a warning if any k-point is nonzero in a truncated direction.
  subroutine check_trunc_kpts(itruncflag, kp)
    integer, intent(in) :: itruncflag
    type(kpoints), intent(in) :: kp
    if(peinf%inode /= 0) return
   
    select case(itruncflag)
    case(0) ! none
    case(2) ! spherical
      if(any(abs(kp%rk(1:3,:)) > TOL_Zero)) &
        write(0,'(a)') 'WARNING: spherical truncation should not be done with k-sampling in any direction.'
      ! there is one exception: Hartree-Fock with the Spencer-Alavi scheme (Phys. Rev. B 77, 193110 (2008))
    case(4) ! cell_wire
      if(any(abs(kp%rk(1:2,:)) > TOL_Zero)) &
        write(0,'(a)') 'WARNING: cell_wire truncation should not be done with k-sampling in the x- or y-directions.'
    case(5) ! cell_box
      if(any(abs(kp%rk(1:3,:)) > TOL_Zero)) &
        write(0,'(a)') 'WARNING: cell_box truncation should not be done with k-sampling in any direction.'
    case(6) ! cell_slab
      if(any(abs(kp%rk(3,:)) > TOL_Zero)) &
        write(0,'(a)') 'WARNING: cell_slab truncation should not be done with k-sampling in the z-direction.'
    case default
      write(0,*) 'itruncflag = ', itruncflag
      call die("Unknown truncation type.")
    end select
   
    return
  end subroutine check_trunc_kpts
  !-----------------------------------------------------------------
  !> Set the number of matrices in the epsmat file depending on the matrix,
  !! type, frequencyd dependency, and flavor. Assumes that pol%freq_dep and
  !! pol%matrix_type were set.
  subroutine eps_setup_sizes(pol, flavor, nspin)
    type(polarizability), intent(inout) :: pol
    integer, intent(in) :: flavor
    integer, intent(in) :: nspin
   
    ! FHJ: old bevavior: FF calculations always compute complex dielectric
    ! matrices, and need two matrices (retarded+advanced).
    pol%has_advanced = .false.
    if (pol%freq_dep/=0 .and. flavor==2) pol%has_advanced = .true.
    ! FHJ: New behavior: store just the retarded matrix + v(q). Get advanced
    ! matrix from symmetry relations
    if (pol%use_hdf5) pol%has_advanced = .false.
    pol%matrix_flavor = flavor
    if (pol%freq_dep/=0) pol%matrix_flavor = 2
    ! FHJ: Set number of matrices depending on nspin (chimat is spin resolved)
    pol%nmatrix = 1
    if (pol%has_advanced) pol%nmatrix = 2
    if (pol%matrix_type==2) pol%nmatrix = pol%nmatrix*nspin
   
  end subroutine eps_setup_sizes
end module input_utils_m
