!===================================================================
!
! Modules:
!
! 1. message_m Originally By DAS
!
! die routine "gracefully" kills the computation.
! alloc_check writes warnings and errors for memory allocation problems.
! write_memory_usage provides a memory report at the end of a run.
! open_file opens a file unit, and writes an error if it does not work.
!
!===================================================================
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
module message_m
  use, intrinsic :: iso_c_binding
  use nrtype_m
  use peinfo_m
  use os_m
  implicit none
  private
  interface
    subroutine abort() bind(c)
    end subroutine abort
  end interface
  public :: &
    die, &
    alloc_check, &
    write_memory_usage, &
    open_file, &
    close_file, &
    get_optional_arg, &
    TRUNC, &
    operator(+)
  !> these names are as short as practical to avoid lines being too long
  integer, public :: alc !< allocation status for safe(de)allocate
  integer*8, public :: sz !< size returned from sizeof for safe(de)allocate
  !> set to .true. to print array size in bytes, .false. in KB/MB/GB
  logical, parameter :: reportsizeexact = .false.
  interface operator (+)
    module procedure cat
  end interface operator (+)
  interface get_optional_arg
    module procedure get_optional_arg_integer, get_optional_arg_logical, &
      get_optional_arg_real, get_optional_arg_complex
  end interface
contains
  !> remove trailing and leading whitespace from a string
  function TRUNC(s)
    character(len=*), intent(in) :: s
    character(len=len_trim(adjustl(s))) :: TRUNC
    TRUNC = trim(adjustl(s))
  end function TRUNC
  !-----------------------------------------------------------
  !> concatenate two strings
  function cat(str1, str2)
    character(len=*), intent(in) :: str1
    character(len=*), intent(in) :: str2
    character(len=len(str1) + len(str2)) :: cat
    cat = str1//str2
  end function cat
  !-----------------------------------------------------------
  subroutine die(str, only_root_writes)
    character (len=*), intent(in) :: str
    logical, optional, intent(in) :: only_root_writes
    logical :: should_write, should_write_prefix, is_open
    should_write = peinf%inode == 0
    should_write_prefix = .false.
    if (present(only_root_writes)) then
      if (.not.only_root_writes) then
        should_write = .true.
        should_write_prefix = peinf%npes > 1
      endif
    endif
    ! FHJ: FLUSH is not really reliable because the OS might cache the stdout.
    ! Sleeping for 1s is the best solution I found to make the output clean,
    ! otherwise the error message would show up before regular output.
    call sleep(1d0)
    ! FHJ: if we are not writing, wait 5s for the root node to get here and
    ! write the error message. If the root doesn`t get here, we all print the
    ! error messsage anyways and die.
    if (.not.should_write) then
      call sleep(5d0)
    endif
    write(0,*)
    if (should_write_prefix) write(0, '(a, i6, a)', advance='no') "From proc ", peinf%inode, ": "
    write(0, '(2a)') "ERROR: ", TRUNC(str)
    write(0,*)
    FLUSH(0)
    ! FHJ: Use libc`s abort funciton.
    ! Fortran`s `stop` or `error stop` don`t give us the traceback!
    call abort()
  end subroutine die
  !---------------------------------------------------------------------------------------------------
  subroutine alloc_check(status, size, name, file, line, flag)
    integer, intent(in) :: status
    !> on some platforms there is a different return value for sizeof if build is 64-bit
    integer*8, intent(in) :: size
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical, intent(in) :: flag
    real(DP) :: sizekb,sizemb,sizegb
    character(len=16) :: prefix
    character(len=32) :: sizestr
    sizekb = dble(size) / dble(1024)
    sizemb = sizekb / dble(1024)
    sizegb = sizemb / dble(1024)
    if (sizekb.le.1.0d1.or.reportsizeexact) then
      write(sizestr,'(i20,1x,"bytes")')size
    elseif (sizemb.le.1.0d1) then
      write(sizestr,'(f20.3,1x,"KB")')sizekb
    elseif (sizegb.le.1.0d1) then
      write(sizestr,'(f20.3,1x,"MB")')sizemb
    else
      write(sizestr,'(f20.3,1x,"GB")')sizegb
    endif
    if (flag) then
      prefix = "Allocation"
      peinf%mymem = peinf%mymem + size
      peinf%mymaxmem = max(peinf%mymaxmem, peinf%mymem)
    else
      prefix = "Deallocation"
      peinf%mymem = peinf%mymem - size
    endif
    if (peinf%verb_debug .and. sizemb>100 .and. peinf%inode==0) then
      write(0,347) trim(prefix), trim(name), TRUNC(sizestr), &
        trim(file), line
    endif
    if(size .lt. 0 .and. peinf%inode .eq. 0) then
      write(0,345) trim(prefix), trim(name), TRUNC(sizestr), &
        trim(file), line
    endif
    if(status .eq. 0) return
    write(0,346) trim(prefix), trim(name), peinf%inode, &
      trim(file), line
    write(0,348) status, TRUNC(sizestr)
    FLUSH(0)
    call die('Allocation failure.')
345 format(1x,"WARNING:",1x,a,1x,"of array",1x,a,1x, &
      "of size",1x,a,/,3x,"in file",1x,a,1x,"at line",i5, &
      1x,"may fail.",/)
346 format(1x,"ERROR:",1x,a,1x,"of array",1x,a,1x, &
      "on processor",i5,/,3x,"in file",1x,a,1x, &
      "at line",i5,1x,"failed.")
347 format(1x,"NOTICE:",1x,a,1x,"of array",1x,a,1x, &
      "of size",1x,a,/,3x,"in file",1x,a,1x,"at line",i5, &
      1x,"occurring.",/)
348 format(3x,"Allocation status =",i4,",",1x,"Array size =",1x,a)
  end subroutine alloc_check
  !---------------------------------------------------------------------------------------------------
  subroutine write_memory_usage()
    ! the memory is not tracked if not in debug mode, so everything would just be zero
    return
  end subroutine write_memory_usage
  !---------------------------------------------------------------------------------------------------
  !> This is a wrapper to the Fortran 'open' statement, to provide clear error-handling.
  !> arguments 'status', 'form', 'position' have the same meaning as for 'open'.
  !> 'iostat', if provided, will make 'iostat' be passed to 'open', and return its value, rather
  !> than writing a message, if there is an error (e.g. status='old' but file does not exist, or
  !> status='new' but file does exist).
  subroutine open_file(unit, file, status, form, position, access, iostat)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: status
    character(len=*), optional, intent(in) :: form
    character(len=*), optional, intent(in) :: position
    character(len=*), optional, intent(in) :: access
    integer, optional, intent(out) :: iostat
    integer :: ierr, unit_other
    character*80 :: form_, position_, name, unit_str, access_, unit_other_str
    character*200 :: string
    logical :: is_open, does_exist
    if(unit == 0) call die("You may not open unit 0, it is reserved for standard error.")
    if(unit == 5) call die("You may not open unit 5, it is reserved for standard input.")
    if(unit == 6) call die("You may not open unit 6, it is reserved for standard output.")
    ! Cray Fortran has its own reserved units: http://docs.cray.com/books/S-3695-35/html-S-3695-35/pdollsmg.html
    if(unit == 100) call die("You may not open unit 100, it is reserved for standard input (crayftn).")
    if(unit == 101) call die("You may not open unit 101, it is reserved for standard output (crayftn).")
    if(unit == 102) call die("You may not open unit 102, it is reserved for standard error (crayftn).")
    ! these issues would be caught below too, but we can give more helpful messages than just an error code
    inquire(unit = unit, opened = is_open, name = name)
    if(is_open) then
      write(string,'(3a,i6,3a)') "Cannot open file '", TRUNC(file), "' on unit ", unit, &
        ": unit already open for file '", TRUNC(name), "'."
      call die(string)
    endif
    if((trim(status) == 'old' .or. trim(status) == 'OLD') .and. .not. present(iostat)) then
      inquire(file = TRUNC(file), exist = does_exist, opened = is_open, number = unit_other)
      if(.not. does_exist) call die("Cannot open file '" // TRUNC(file) // "' for reading: does not exist.")
      if(is_open) then
        write(unit_str,*) unit
        write(unit_other_str,*) unit_other
        call die("Cannot open file '" // TRUNC(file) // "' for reading on unit " // TRUNC(unit_str) &
          // ": already opened on unit " // TRUNC(unit_other_str) // ".")
      endif
! From the Fortran 95 standard, Section 9.3.4:
!
! If a file is already connected to a unit, execution of an OPEN
! statement on that file and a different unit is not permitted.
!
! From the Fortran 77 Standard, Section 12.3.2:
!
! A unit must not be connected to more than one file at the same time,
! and a file must not be connected to more than one unit at the same time.
    endif
    if((trim(status) == 'new' .or. trim(status) == 'NEW') .and. .not. present(iostat)) then
      inquire(file = TRUNC(file), exist = does_exist)
      if(does_exist) call die("Cannot open file '" // TRUNC(file) // "' for writing as 'new': already exists.")
    endif
    form_ = 'formatted'
    if(present(form )) form_ = form
    position_ = 'asis'
    if(present(position)) position_ = position
    access_ = 'sequential'
    if(present(access)) access_ = access
    ! passing the optionals to 'open' if not given to this routine does not work!
    open(unit=unit, file = TRUNC(file), form=trim(form_), access=trim(access_), &
      position=trim(position_), status=trim(status), iostat=ierr)
    if(present(iostat)) then
      iostat = ierr
    else if(ierr /= 0) then
      write(string,'(5a,i4)') "Failed to open file '", TRUNC(file), "' as status ", trim(status), " with error ", ierr
      call die(string)
    endif
    return
  end subroutine open_file
  !---------------------------------------------------------------------------------------------------
  subroutine close_file(unit, delete)
    integer, intent(in) :: unit
    logical, optional, intent(in) :: delete
    character*80 :: string, status
    logical :: is_open
    integer :: ierr
    if(unit == 0) call die("You may not close unit 0, it is reserved for standard error.")
    if(unit == 5) call die("You may not close unit 5, it is reserved for standard input.")
    if(unit == 6) call die("You may not close unit 6, it is reserved for standard output.")
    ! Cray Fortran has its own reserved units: http://docs.cray.com/books/S-3695-35/html-S-3695-35/pdollsmg.html
    if(unit == 100) call die("You may not close unit 100, it is reserved for standard input (crayftn).")
    if(unit == 101) call die("You may not close unit 101, it is reserved for standard output (crayftn).")
    if(unit == 102) call die("You may not close unit 102, it is reserved for standard error (crayftn).")
    ! these issues would be caught below too, but we can give more helpful messages than just an error code
    inquire(unit = unit, opened = is_open, iostat = ierr)
    if(ierr /= 0) then
      write(string,'(a,i6,a,i4)') "inquire in close_file failed for unit ", unit, " with error ", ierr
      call die(string)
    endif
    if(.not. is_open) then
      write(string,'(a,i6,a)') "Cannot close unit ", unit, ": not open."
      call die(string)
    endif
    status = 'keep'
    if(present(delete)) then
      if(delete) status = 'delete'
    endif
    close(unit=unit, status=trim(status), iostat=ierr)
    if(ierr /= 0) then
      write(string,'(a,i6,a,i4)') "Failed to close unit ", unit, " with error ", ierr
      call die(string)
    endif
    return
  end subroutine close_file
  !> Return `optional_argument` if present, otherwise `default_value`
  function get_optional_arg_integer(default_value, optional_argument) result(output)
    integer, intent(in) :: default_value
    integer, optional, intent(in) :: optional_argument
    integer :: output
    output = default_value
    if (present(optional_argument)) output = optional_argument
  end function get_optional_arg_integer
  !> Return `optional_argument` if present, otherwise `default_value`
  function get_optional_arg_logical(default_value, optional_argument) result(output)
    logical, intent(in) :: default_value
    logical, optional, intent(in) :: optional_argument
    logical :: output
    output = default_value
    if (present(optional_argument)) output = optional_argument
  end function get_optional_arg_logical
  !> Return `optional_argument` if present, otherwise `default_value`
  function get_optional_arg_real(default_value, optional_argument) result(output)
    real(DP), intent(in) :: default_value
    real(DP), optional, intent(in) :: optional_argument
    real(DP) :: output
    output = default_value
    if (present(optional_argument)) output = optional_argument
  end function get_optional_arg_real
  !> Return `optional_argument` if present, otherwise `default_value`
  function get_optional_arg_complex(default_value, optional_argument) result(output)
    complex(DPC), intent(in) :: default_value
    complex(DPC), optional, intent(in) :: optional_argument
    complex(DPC) :: output
    output = default_value
    if (present(optional_argument)) output = optional_argument
  end function get_optional_arg_complex
end module message_m
