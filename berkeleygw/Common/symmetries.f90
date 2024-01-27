!===============================================================================
!
! Module:
!
! symmetries_m Originally By DAS 12/20/2011
!
! Find symmetry operations from lattice vectors and atomic positions,
! using spglib 1.0.9. For use in mean-field codes such as EPM and SIESTA wrapper.
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
module symmetries_m
  use global_m
  use misc_m
  use sort_m
  implicit none
  private
  public :: get_symmetries
  real(DP), parameter :: symprec = 1d-5
  ! this is the same tolerance used by the Quantum ESPRESSO routines
  interface
    ! these functions are defined in spglib_f.c
    subroutine spg_get_multiplicity_f(size, lattice, position, types, num_atom, symprec)
      use global_m
      implicit none
      integer, intent(out) :: size
      real(DP), intent(in) :: lattice
      real(DP), intent(in) :: position
      integer, intent(in) :: types
      integer, intent(in) :: num_atom
      real(DP), intent(in) :: symprec
    end subroutine spg_get_multiplicity_f
    subroutine spg_get_symmetry_f(nsym, rotation, translation, max_size, lattice, position, types, num_atom, symprec)
      use global_m
      implicit none
      integer, intent(out) :: nsym
      integer, intent(out) :: rotation
      real(DP), intent(out) :: translation
      integer, intent(in) :: max_size
      real(DP), intent(in) :: lattice
      real(DP), intent(in) :: position
      integer, intent(in) :: types
      integer, intent(in) :: num_atom
      real(DP), intent(in) :: symprec
    end subroutine spg_get_symmetry_f
    subroutine spg_get_international_f(spacegroup, symbol, lattice, position, types, num_atom, symprec)
      use global_m
      implicit none
      integer, intent(out) :: spacegroup
      character*11, intent(out) :: symbol
      real(DP), intent(in) :: lattice
      real(DP), intent(in) :: position
      integer, intent(in) :: types
      integer, intent(in) :: num_atom
      real(DP), intent(in) :: symprec
    end subroutine spg_get_international_f
  end interface
contains
  subroutine get_symmetries(nat, atyp, apos, avec, nfft, cell_symmetry, ntran, mtrx, tnp, spacegroup, symbol)
    integer, intent(in) :: nat !< number of atoms
    integer, intent(in) :: atyp(:) !< atomic species
    real(DP), intent(in) :: apos(:,:) !< (1:3,1:nat) atomic positions in crystal coordinates
    real(DP), intent(in) :: avec(3,3) !< lattice vectors in real space
    integer, intent(in) :: nfft(3) !< FFT grid (if = 0, not used)
    integer, intent(out) :: cell_symmetry !< 0 = cubic, 1 = hexagonal
    integer, intent(out) :: ntran !< number of symmetry operations
    integer, intent(out) :: mtrx(:,:,:) !< rotation matrices (3, 3, 48)
    real(DP), intent(out) :: tnp(:,:) !< frational translations (3, 48)
    integer, intent(out) :: spacegroup !< spacegroup international index
    character, intent(out) :: symbol*21 !< spacegroup symbol in Schoenflies notation
    integer :: isym, ii, jj, ntran_temp, ntran_temp2, identity(3,3)
    real(DP) :: C_avec(3,3), frac_fft(3)
    integer, allocatable :: mtrx_inv(:,:,:) ! allocatable since more than 48 ops may be returned
    real(DP), allocatable :: tnp_temp(:,:)
    logical :: use_this_sym, disable_frac, found_identity
   
    ! transpose input for C call
    do ii = 1, 3
      C_avec(ii, 1:3) = avec(1:3, ii)
    enddo
    call spg_get_international_f(spacegroup, symbol, C_avec(1, 1), apos(1, 1), atyp(1), nat, symprec)
    ! http://en.wikipedia.org/wiki/Trigonal_crystal_system: 143-167
    ! http://en.wikipedia.org/wiki/Hexagonal_crystal_system: 168-194
    ! All others are cubic.
    if(spacegroup >= 143 .and. spacegroup <= 194) then
      cell_symmetry = 1 ! hexagonal
    else
      cell_symmetry = 0 ! cubic
    endif
    disable_frac = .false.
    call spg_get_multiplicity_f(ntran_temp, C_avec(1,1), apos(1,1), atyp(1), nat, symprec)
    allocate(mtrx_inv (3, 3, ntran_temp))
    allocate(tnp_temp (3, ntran_temp))
    ! we need to check that it is not a supercell, as in the QE routine (sgam_at)
    ! they disable fractional translations if the identity has one, because the sym ops might not form a group.
    ! spglib may return duplicate operations in this case!
    call spg_get_symmetry_f(ntran_temp2, mtrx_inv(1,1,1), tnp_temp(1,1), ntran_temp, C_avec(1,1), apos(1,1), atyp(1), nat, symprec)
    if(ntran_temp2 /= ntran_temp) call die("Inconsistent number of symmetries from spglib. Internal error.")
! for debugging: write all ops
! do isym = 1, ntran_temp2
! write(6,'(10i3, 3f6.2)') isym, mtrx_inv(:,:,isym), tnp_temp(:,isym)
! enddo
    if(ntran_temp > 48) then
      disable_frac = .true.
      write(0,'(a,i6,a)') "Number of symmetry operations = ", ntran_temp, " > 48"
    endif
    found_identity = .false.
    identity = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), shape(identity))
    do isym = 1, ntran_temp
      if(all(mtrx_inv(1:3,1:3,isym) == identity(1:3, 1:3))) then
        found_identity = .true.
        if(any(abs(tnp_temp(1:3, isym)) > TOL_Zero)) then
          disable_frac = .true.
          write(0,'(a,3f12.6)') 'Identity has a fractional translation ', tnp_temp(1:3, isym)
        endif
      endif
    enddo
    if(.not. found_identity) call die("Symmetries internal error: Identity is missing from symmetry operations.")
    if(disable_frac) then
      write(0,'(a)') "WARNING: Disabling fractional translations. System appears to be a supercell."
    endif
    ! spglib does not consider an FFT grid. below is based on QE routine sgam_at.
    ntran = 0
    do isym = 1, ntran_temp
      if(all(nfft(1:3) /= 0)) then
        ! check that rotation matrix is compatible with FFT grid
        use_this_sym = &
          mod(mtrx_inv(2, 1, isym) * nfft(1), nfft(2)) == 0 .and. &
          mod(mtrx_inv(3, 1, isym) * nfft(1), nfft(3)) == 0 .and. &
          mod(mtrx_inv(1, 2, isym) * nfft(2), nfft(1)) == 0 .and. &
          mod(mtrx_inv(3, 2, isym) * nfft(2), nfft(3)) == 0 .and. &
          mod(mtrx_inv(1, 3, isym) * nfft(3), nfft(1)) == 0 .and. &
          mod(mtrx_inv(2, 3, isym) * nfft(3), nfft(2)) == 0
        ! check that fractional translation is compatible with FFT grid
        frac_fft(1:3) = tnp_temp(1:3, isym) * nfft(1:3)
        use_this_sym = use_this_sym .and. all(abs(frac_fft(1:3) - nint(frac_fft(1:3))) / nfft(1:3) < symprec)
      else
        ! if FFT grid is supplied as zero, we accept all operations
        use_this_sym = .true.
      endif
      ! make sure fractional translations are in the right range, just a convention
      ! this makes the results agree with the ESPRESSO routines, but really makes no practical difference
      do jj = 1, 3
        if (tnp_temp(jj, isym).ge.TOL_Zero+0.5d0) &
          tnp_temp(jj, isym)=tnp_temp(jj, isym)-dble(int(tnp_temp(jj, isym)+0.5d0))
        if (tnp_temp(jj, isym).lt.TOL_Zero-0.5d0) &
          tnp_temp(jj, isym)=tnp_temp(jj, isym)-dble(int(tnp_temp(jj, isym)-0.5d0))
      enddo
      if(disable_frac) then
        if(any(abs(tnp_temp(1:3,isym)) > TOL_Zero)) use_this_sym = .false.
      endif
      if(use_this_sym) then
        ntran = ntran + 1
        if(ntran > 48) call die("Internal error: There are more than 48 accepted symmetry operations.")
        ! this could only happen for a supercell, and we are supposed to have handled that situation already
        ! if accepted, add this operation to the list
        call invert_matrix_int(mtrx_inv(1:3, 1:3, isym), mtrx(1:3, 1:3, ntran))
        tnp(1:3, ntran) = 2 * PI_D * tnp_temp(1:3, isym)
      endif
    enddo
    if(allocated(mtrx_inv))then;deallocate(mtrx_inv);endif
    if(allocated(tnp_temp))then;deallocate(tnp_temp);endif
    call make_identity_symmetry_first(ntran, mtrx, tnp)
   
    return
  end subroutine get_symmetries
end module symmetries_m
