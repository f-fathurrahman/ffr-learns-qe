!=================================================================================
!
! Module read_matrix
!
! (1) read_matrix_d() Originally by JRD Last Modified 5/1/2008 (JRD)
!
! This program reads a distributed matrix like chimat or epsmat to file.
!
! (2) read_matrix_f() Originally by JRD Last Modified 9/10/2010 (gsm)
!
! Modification of read_matrix_d for full-frequency.
!
!=================================================================================

module read_matrix_m
  use global_m
  use scalapack_m
  implicit none
  private
  public :: &
    read_matrix_d, &
    read_matrix_d_hdf5, &
    read_matrix_f, &
    read_matrix_f_hdf5
contains
subroutine read_matrix_d(scal,matrix,nmtx,iunit)
  type (scalapack), intent(in) :: scal
  real(DP), intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: iunit
 
  call read_matrix_d_(scal,matrix,nmtx,iunit=iunit)
 
end subroutine read_matrix_d
subroutine read_matrix_d_hdf5(scal,matrix,nmtx,fname,iq,is)
  type (scalapack), intent(in) :: scal
  real(DP), intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  character(len=*), intent(in) :: fname
  integer, intent(in) :: iq
  integer, intent(in) :: is
 
  call read_matrix_d_(scal,matrix,nmtx,fname=fname,iq=iq,is=is)
 
end subroutine read_matrix_d_hdf5
subroutine read_matrix_d_(scal,matrix,nmtx,iunit,fname,iq,is)
  type (scalapack), intent(in) :: scal
  real(DP), intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in), optional :: iunit
  character(len=*), intent(in), optional :: fname
  integer, intent(in), optional :: iq
  integer, intent(in), optional :: is
  integer :: ii, jj
  logical :: use_hdf5
 
  if (.not.present(iunit).and..not.(present(fname).and.present(iq))) then
    call die("Not enough arguments to read_matrix_d_", only_root_writes=.true.)
  endif
  if (present(iunit).and.(present(fname).or.present(iq))) then
    call die("Too many arguments to read_matrix_d_", only_root_writes=.true.)
  endif
  if ((present(fname).or.present(iq)).and..not.(present(fname).and.present(iq))) then
    call die("Inconsistent arguments to read_matrix_d_", only_root_writes=.true.)
  endif
  use_hdf5 = present(fname).and.present(iq)
  if (use_hdf5) then
    call die("read_matrix_d_ was not compiled with HDF5 support.", only_root_writes=.true.)
  endif
  if (peinf%verb_debug .and. peinf%inode==0) then
    if (use_hdf5) then
      write(6,*) 'Reading matrix: ', nmtx, fname
    else
      write(6,*) 'Reading matrix: ', nmtx, iunit
    endif
    write(6,*)
  endif
  if(peinf%inode .eq. 0) then
    do jj = 1, nmtx
      if (use_hdf5) then
      else
        read(iunit) (matrix(ii, jj), ii = 1, nmtx)
      endif
    enddo
  endif
 
  return
end subroutine read_matrix_d_
!=================================================================================
!> FHJ: Front end for read_matrix_f_ for Fortran binary files. See that routine for more info.
subroutine read_matrix_f(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, iunit, advanced)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq
  integer, intent(in) :: nfreq_in_group
  complex(DPC), intent(out) :: retarded(:,:,:) !< (nfreq_in_group,scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: nfreq_group
  integer, intent(in) :: iunit
  complex(DPC), optional, intent(out) :: advanced(:,:,:) !< (nfreq_in_group,scal%npr,scal%npc)
 
  call read_matrix_f_(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, iunit=iunit, advanced=advanced)
 
end subroutine read_matrix_f
!> FHJ: Front end for read_matrix_f_ for HDF5 files. See that routine for more info.
subroutine read_matrix_f_hdf5(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, fname, iq, is, advanced)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq
  integer, intent(in) :: nfreq_in_group
  complex(DPC), intent(out) :: retarded(:,:,:) !< (nfreq_in_group,scal%npr,scal%npc)
  integer, intent(in) :: nmtx
  integer, intent(in) :: nfreq_group
  character(len=*), intent(in) :: fname
  integer, intent(in) :: iq
  integer, intent(in) :: is
  complex(DPC), optional, intent(out) :: advanced(:,:,:) !< (nfreq_in_group,scal%npr,scal%npc)
 
  call read_matrix_f_(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, &
    fname=fname, iq=iq, is=is, advanced=advanced)
 
end subroutine read_matrix_f_hdf5
!> FHJ: This routines the full-frequency chiR/epsR matrix from a file, and
!! optionally chiA/epsA (note: you shouldn`t really need chiA, ever...)
!! If using HDF5, we only read the retarded part. If legacy
!! Fortran binary, we read the retarded and skip the advanced. The final
!! matrix will be distributed in a ScaLAPACK layout given by scal. Note that
!! this routine is pretty innefficient, but this is not a core component
!! of BGW as it`s only used if you read_chi or use the eps*omega utility.
subroutine read_matrix_f_(scal, nfreq, nfreq_in_group, retarded, nmtx, &
  nfreq_group, iunit, fname, iq, is, advanced)
  type(scalapack), intent(in) :: scal
  integer, intent(in) :: nfreq
  integer, intent(in) :: nfreq_in_group
  complex(DPC), intent(out) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  integer, intent(in) :: nmtx
  integer, intent(in) :: nfreq_group
  integer, intent(in), optional :: iunit
  character(len=*), intent(in), optional :: fname
  integer, intent(in), optional :: iq
  integer, intent(in), optional :: is
  complex(DPC), intent(out), optional :: advanced(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  integer :: ig_glob, igp_glob, ifreq, ifreq_para, freq_grp_ind
  logical :: use_hdf5, want_advanced
 
  want_advanced = .false.
  if (.not.present(iunit).and..not.(present(fname).and.present(iq))) then
    call die("Not enough arguments to read_matrix_f_", only_root_writes=.true.)
  endif
  if (present(iunit).and.(present(fname).or.present(iq))) then
    call die("Too many arguments to read_matrix_f_", only_root_writes=.true.)
  endif
  if ((present(fname).or.present(iq)).and..not.(present(fname).and.present(iq))) then
    call die("Inconsistent arguments to read_matrix_f_", only_root_writes=.true.)
  endif
  use_hdf5 = present(fname).and.present(iq)
  if (use_hdf5) then
    call die("read_matrix_f_ was not compiled with HDF5 support.", only_root_writes=.true.)
  endif
  if (peinf%verb_debug .and. peinf%inode==0) then
    if (use_hdf5) then
      write(6,*) ' Reading matrix: ', nmtx, fname
    else
      write(6,*) ' Reading matrix: ', nmtx, iunit
    endif
    write(6,*)
  endif
  if (peinf%npes>1) then
  else !USESCALAPACK & peinf%npes>1
    if (use_hdf5) then
    else !use_hdf5
      do igp_glob = 1, nmtx
        do ig_glob = 1, nmtx
          read(iunit) (retarded(ig_glob, igp_glob, ifreq), ifreq = 1, nfreq)
        enddo
      enddo
    endif !use_hdf5
  endif !peinf%npes==1
 
contains
  !> Saves the global `tempcolR`, such as column of the dielectric matrix,
  !! to the local array corresponding to the distributed matrix `retarted`.
  subroutine save_local_matrix(colR, colA)
    complex(DPC), intent(in) :: colR(:)
    complex(DPC), intent(in), optional :: colA(:)
    integer :: igp_loc, ig_glob, ig_loc
   
    if ((freq_grp_ind==peinf%igroup_f) .and. (scal%mypcol==&
         INDXG2P(igp_glob, scal%nbl, scal%mypcol, 0, scal%npcol))) then
      igp_loc = INDXG2L(igp_glob, scal%nbl, scal%mypcol, 0, scal%npcol)
      do ig_loc = 1, scal%npr
        ig_glob = INDXL2G(ig_loc, scal%nbl, scal%myprow, 0, scal%nprow)
        retarded(ig_loc, igp_loc, ifreq_para) = colR(ig_glob)
        if (want_advanced) then
          advanced(ig_loc, igp_loc, ifreq_para) = colA(ig_glob)
        endif
      enddo
    endif
   
  end subroutine save_local_matrix
end subroutine read_matrix_f_
end module read_matrix_m
