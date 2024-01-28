!=============================================================================
!
! Routines:
!
! (1) inread_kernel()
!
! Read input parameters from file kernel.inp
!
! input: none
!
! output: xct%nvb_co
! xct%ncb_co
! xct%ecute energy cutoff used to calculate the
! interaction matrices (see Rohlfing and Louie PRB 2000 eq. 43-45)
! xct%read_kpoints tells whether the list of k-points is taken
! from "WFN_fi" or from the input "kpoints"
! (useful if you are working with partial samplings)
! flagbz = 0 use all symmetries to unfold the Brillouin zone
! = 1 do not unfold the Brillouin zone (default)
!
!=============================================================================
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
module inread_kernel_m
  use algos_kernel_m
  use global_m
  use inread_common_m
  use references_m
  implicit none
  private
  public :: &
    inread_kernel
contains
subroutine inread_kernel(xct,flagbz,qg)
  type (xctinfo), intent(out) :: xct
  integer, intent(out) :: flagbz
  type (grid), intent(out) :: qg
  character*256 :: keyword,line,errmsg,blockword
  integer :: iostat
  integer :: ii,jj,itestq
  real(DP) :: qpt_read(3,MAX_KPTS),div
  logical :: found
 
!-------------------------
! Set default values
  xct%ecute=0.d0
  xct%ecutg=0.d0
  xct%read_kpoints = .false.
  xct%icutv=TRUNC_NONE
  xct%iscreen=SCREEN_SEMICOND
  xct%iwritecoul=0
  xct%truncval(:)=0.0d0
  xct%bLowComm=.false.
  xct%shift(:)=0.d0
  xct%qflag=1
  xct%finiteq(:)=0.d0
  xct%ilowmem=0
  flagbz=1
  xct%efermi_input=0.0d0
  xct%rfermi=.true.
  xct%freplacebz=.false.
  xct%fwritebz=.false.
  xct%degeneracy_check_override=.true.
  xct%die_outside_sphere=.false.
  xct%extended_kernel=.false.
  xct%theory=0
  xct%coulomb_mod%short_range_frac_fock=1.0d0
  xct%coulomb_mod%long_range_frac_fock=1.0d0
  xct%coulomb_mod%screening_length=0.0d0
  xct%coul_mod_flag=.false.
  xct%patched_sampling_co=.false.
  xct%energy_loss=.false.
  xct%screen_exchange=.false.
  xct%use_hdf5 = .false.
  xct%use_wfn_hdf5 = .false.
  ! default is to read at least 128 Mb of data, using wfn_hdf5_min_band_block a different
  ! minimum value of band can be set
  xct%wfn_hdf5_min_band_block = -1
  xct%nspin = -1 ! Determined in input_kernel.f90
  xct%nspinor = -1 ! Determined in input_kernel.f90
! Set default values for algo acceleration
  call set_algos_to_cpu()
!-----------------------------------
! Never ending loop...
  do while(0.eq.0)
! Actually the loop ends when the end of the file is reached
    read(8,'(a256)',iostat=iostat) line
    if(iostat < 0) exit
! Skip comment lines
    if(len_trim(line).eq.0) cycle
    if(line(1:1).eq.'#') cycle
! Determine keyword:
    keyword=line(1:scan(line," ")-1)
    line=adjustl(line(scan(line," ")+1:256))
    if(trim(keyword).eq.'verbosity') then
      read(line,*,err=110) peinf%verbosity
    elseif(trim(keyword).eq.'number_val_bands') then
      read(line,*,err=110) xct%nvb_co
    elseif(trim(keyword).eq.'number_cond_bands') then
      read(line,*,err=110) xct%ncb_co
    elseif(trim(keyword).eq.'dont_use_hdf5') then
      xct%use_hdf5 = .false.
    elseif(trim(keyword).eq.'extended_kernel') then
      xct%extended_kernel=.true.
    elseif(trim(keyword).eq.'patched_sampling_co') then
      xct%patched_sampling_co=.true.
    elseif(trim(keyword).eq.'exciton_Q_shift') then
      read(line,*,err=110) xct%qflag, xct%finiteq(1),xct%finiteq(2),xct%finiteq(3)
      xct%energy_loss=.true.
    elseif(trim(keyword).eq.'energy_loss') then
      xct%energy_loss=.true.
    elseif(trim(keyword).eq.'screened_coulomb_cutoff') then
      read(line,*,err=110) xct%ecute
    elseif(trim(keyword).eq.'bare_coulomb_cutoff') then
      read(line,*,err=110) xct%ecutg
    elseif(trim(keyword).eq.'read_kpoints') then
      xct%read_kpoints = .true.
    elseif(trim(keyword).eq.'write_vcoul') then
      xct%iwritecoul=1
    elseif(trim(keyword).eq.'low_comm') then
      xct%bLowComm=.true.
    elseif(trim(keyword).eq.'low_memory') then
      xct%ilowmem=1
    elseif(trim(keyword).eq.'high_memory') then
      xct%ilowmem=-1
    elseif(trim(keyword).eq.'fermi_level') then
      read(line,*,err=110) xct%efermi_input
    elseif(trim(keyword).eq.'fermi_level_absolute') then
      xct%rfermi=.false.
    elseif(trim(keyword).eq.'fermi_level_relative') then
      xct%rfermi=.true.
    elseif(trim(keyword).eq.'no_symmetries_coarse_grid') then
      flagbz = 1
    elseif(trim(keyword).eq.'use_symmetries_coarse_grid') then
      flagbz = 0
    elseif(trim(keyword).eq.'fullbz_replace') then
      xct%freplacebz=.true.
    elseif(trim(keyword).eq.'fullbz_write') then
      xct%fwritebz=.true.
    elseif(trim(keyword).eq.'use_wfn_hdf5') then
    elseif(trim(keyword).eq.'wfn_hdf5_min_band_block') then
      read(line,*,err=110) xct%wfn_hdf5_min_band_block
! There is no problem with degeneracy in kernel, actually.
! elseif(trim(keyword).eq.'degeneracy_check_override') then
! xct%degeneracy_check_override=.true.
    elseif(trim(keyword).eq.'die_outside_sphere') then
      xct%die_outside_sphere=.true.
    elseif(trim(keyword).eq.'ignore_outside_sphere') then
      xct%die_outside_sphere=.false.
    elseif(try_inread_truncation(trim(keyword), trim(line), xct%icutv, xct%truncval(1))) then
      ! subroutine already does the job
    elseif(try_inread_screening(trim(keyword), trim(line), xct%iscreen)) then
      ! subroutine already does the job
    else
      call algos_inread(keyword, line, found)
      if(.not.found) then
        write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found in kernel.inp.'
        call die(errmsg, only_root_writes = .true.)
      end if
    end if
  enddo
! FHJ: How many kernel blocks to compute?
  if (xct%extended_kernel) then
    if (peinf%inode==0) write(6,'(1x,a)') &
      "Calculating the kernel for all possible (n1,n2)->(n1',n2') transitions."
    xct%n1b_co = xct%nvb_co + xct%ncb_co
    xct%n2b_co = xct%n1b_co
    if (xct%ilowmem>=0.and.peinf%inode==0) then
      write (6,'(1x,a)') 'Note: this type of kernel calculation requires the high_memory flag.'
    endif
    xct%ilowmem=-1
  else
    if (peinf%inode==0) write(6,'(1x,a)') &
      "Calculating the kernel only for (v,c)->(v',c') transitions."
    xct%n1b_co = xct%nvb_co
    xct%n2b_co = xct%ncb_co
  endif
  if (peinf%inode==0) write(6,*)
  if (xct%extended_kernel) then
    if (peinf%inode==0) then
      write(6,'(1x,a)') 'WARNING: GPU support for extended-kernel still under development, CPU version will be used'
      write(6,'(1x,a)')
    end if
    call set_algos_to_cpu()
  end if
  call peinfo_set_verbosity()
! JRD: Make a note if we have finite Q
  if (peinf%inode==0) then
    if (xct%qflag==0 .or. xct%qflag==2) then
      write(6,'(1x,a)') 'We are doing a finite-Q calculation'
      write(6,'(1x,a,3(1x,f12.9))') 'Q_shift =', xct%finiteq(:)
      write(6,'(1x,a/)') 'where Q_shift is the NEGATIVE of the exciton center-of-mass (COM) momentum'
      write(6,'(1x,a)')
    else
      write(6,'(1x,a/)') 'We are doing a calculation with zero exciton momentum'
      write(6,'(1x,a)')
    endif
  endif
  ! JRD: What screening is present?
  if (peinf%inode==0) then
    select case (xct%iscreen)
      case (SCREEN_SEMICOND)
        write(6,'(1x,a/)') 'Running with semiconductor screening'
      case (SCREEN_GRAPHENE)
        write(6,'(1x,a/)') 'Running with graphene screening'
      case (SCREEN_METAL)
        write(6,'(1x,a/)') 'Running with metal screening'
      case default
        call die('Unknown screening type', only_root_writes=.true.)
    endselect
  endif
  if(peinf%inode == 0) then
    if(peinf%npes > 1) then
      write(6,803)
    else
      write(6,805)
    endif
  endif
803 format(1x,'We are communicating via MPI',/)
805 format(1x,'We are not communicating',/)
  if (peinf%inode == 0) then
    if (xct%ilowmem == 1) then
      write(6,'(1x,a)') 'We are using the low-memory option.'
    elseif (xct%ilowmem == -1) then
      write(6,'(1x,a)') 'We are using the high-memory option.'
    endif
    if (xct%ilowmem==-1) then
      write(6,'(1x,a)') 'We will reuse wavefunction FFTs.'
    else
      write(6,'(1x,a)') 'We will not reuse wavefunction FFTs.'
    endif
    if(xct%bLowComm) then
      write(6,'(1x,a)') 'We are using the low-communication option.'
    endif
  endif
  call print_truncation_summary(xct%icutv, xct%truncval(1))
  if (peinf%inode==0 .and. (xct%ecute>TOL_ZERO .or. xct%ecutg>TOL_ZERO)) then
    write(6,'(/1x,a/)') 'NOTE: `screened_coulomb_cutoff` and `bare_coulomb_cutoff` are now optional flags.'
  endif
  call verify_gpu_settings()
  call require_reference(REF_Deslippe2012)
  call require_reference(REF_Rohlfing2000)
  ! Truncation of the Coulomb potential: slab and write
  if (xct%icutv==TRUNC_SLAB .or. xct%icutv==TRUNC_WIRE) call require_reference(REF_IsmailBeigi2016)
  ! Finite-Q excitons
  if (xct%qflag/=1) call require_reference(REF_Qiu2015)
  ! Non-TDA
  if (xct%extended_kernel) call require_reference(REF_Shao2016)
 
  return
110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', &
      trim(keyword), '. '
  call die(errmsg, only_root_writes = .true.)
end subroutine inread_kernel
end module inread_kernel_m
