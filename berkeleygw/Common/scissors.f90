!===============================================================================
!
! Modules:
!
! scissors_m Originally By DAS
!
! Routines for scissors corrections to mean-field eigenvalues.
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
module scissors_m
  use global_m
  use splines_m
  implicit none
  private
  public :: &
    scissors_shift, &
    scissors_write, &
    scissors_zero, &
    scissors_function, &
    scissors_inread
contains
  !---------------------------------------------------------------------------------------------------
  !> Scissors operator to get quasi-particle spectrum from LDA eigenvalues
  !! evs, ecs, ev0, ec0 are supposed to be in eV
  !! kp%el is taken to be in Ry however
  !! spl_tck (spline coefficients) in Ry
  subroutine scissors_shift(kp, scis, spl_tck)
    type(kpoints), intent(inout) :: kp
    type(scissors_t), intent(in) :: scis
    type(spline_tck), optional, intent(in) :: spl_tck
    integer :: is, irk, ib
    logical :: has_spline, use_spline
    real(DP) :: spline_emin, spline_emax
   
    has_spline=.false.
    if (present(spl_tck)) then
      if (spl_tck%n>0) then
        has_spline = .true.
        spline_emin = spl_tck%t(1)/ryd
        spline_emax = spl_tck%t(spl_tck%n)/ryd
      endif
    endif
    do is=1,kp%nspin
      do irk=1,kp%nrk
        do ib=1,kp%mnband
          use_spline=.false.
          if (has_spline) then
            if ( (kp%el(ib,irk,is) > spline_emin) .and. &
                 (kp%el(ib,irk,is) < spline_emax) ) then
                 use_spline=.true.
            endif
          endif
          if (use_spline) then
            ! use spline interpolation for scissors ops
            call splev_shift( spl_tck, kp%el(ib,irk,is) )
          else
            ! use regular (linear) scissors ops
            if(ib <= kp%ifmax(irk,is)) then
              kp%el(ib,irk,is) = scissors_function(kp%el(ib, irk, is), scis%val)
            else
              kp%el(ib,irk,is) = scissors_function(kp%el(ib, irk, is), scis%cond)
            endif
          endif
        enddo ! ib (band)
      enddo ! irk (kpoint)
    enddo ! is (spin)
   
    return
  end subroutine scissors_shift
  !---------------------------------------------------------------------------------------------------
  subroutine scissors_write(iunit, scis, suffix)
    integer, intent(in) :: iunit !< unit to which to write
    type(scissors_t), intent(in) :: scis
    character(len=*), optional, intent(in) :: suffix !< identifier to be printed
   
    write(iunit,'(/1x,a)',advance='no') "Scissors parameters"
    if(present(suffix)) then
      write(iunit,'(3a)') ' (', TRUNC(suffix), '):'
    else
      write(iunit,'(a)') ":"
    endif
    write(iunit,'(1x,3(a,f8.4))') &
      '- Valence:    es = ', scis%val%es, ' eV, e0 = ', scis%val%e0, ' eV, edel = ', scis%val%edel
    write(iunit,'(1x,3(a,f8.4))') &
      '- Conduction: es = ', scis%cond%es, ' eV, e0 = ', scis%cond%e0, ' eV, edel = ', scis%cond%edel
   
    return
  end subroutine scissors_write
  !---------------------------------------------------------------------------------------------------
  subroutine scissors_zero(scis)
    type(scissors_t), intent(out) :: scis
   
    scis%val%es = 0d0
    scis%val%e0 = 0d0
    scis%val%edel = 0d0
    scis%cond%es = 0d0
    scis%cond%e0 = 0d0
    scis%cond%edel = 0d0
   
    return
  end subroutine scissors_zero
  !---------------------------------------------------------------------------------------------------
  real(DP) function scissors_function(energy, scis)
    real(DP), intent(in) :: energy !< should be in Ry
    type(sub_scissors_t), intent(in) :: scis
   
    scissors_function = energy + scis%es / ryd + scis%edel * (energy - scis%e0 / ryd)
   
    return
  end function scissors_function
  !---------------------------------------------------------------------------------------------------
  subroutine scissors_inread(keyword, line, scis, found, suffix)
    character(len=*), intent(in) :: keyword, line
    type(scissors_t), intent(inout) :: scis
    logical, intent(out) :: found
    character(len=*), optional, intent(in) :: suffix
    integer :: iostat, trunc_length
    character*120 :: trunc_keyword
   
    if(present(suffix)) then
      if(len(suffix) > 0) then
        trunc_length = len(trim(keyword))
        if(keyword(trunc_length-len(suffix)+1:trunc_length) == suffix) then
          trunc_keyword = keyword(1:trunc_length-len(suffix))
        else
          found = .false.
         
          return
        endif
      endif
    else
      trunc_keyword = keyword
    endif
    found = .true.
    if(trim(trunc_keyword).eq.'cvfit') then
      read(line,*,iostat = iostat) scis%val%es,scis%val%e0,scis%val%edel, &
        scis%cond%es,scis%cond%e0,scis%cond%edel
    elseif(trim(trunc_keyword).eq.'evs') then
      read(line,*,iostat = iostat) scis%val%es
    elseif(trim(trunc_keyword).eq.'ev0') then
      read(line,*,iostat = iostat) scis%val%e0
    elseif(trim(trunc_keyword).eq.'evdel') then
      read(line,*,iostat = iostat) scis%val%edel
    elseif(trim(trunc_keyword).eq.'ecs') then
      read(line,*,iostat = iostat) scis%cond%es
    elseif(trim(trunc_keyword).eq.'ec0') then
      read(line,*,iostat = iostat) scis%cond%e0
    elseif(trim(trunc_keyword).eq.'ecdel') then
      read(line,*,iostat = iostat) scis%cond%edel
    else
      found = .false.
    endif
    if(found .and. iostat /= 0) call die( &
      'Unexpected characters were found while reading the value for the keyword ' &
      // trim(keyword) // '. ', only_root_writes = .true.)
   
    return
  end subroutine scissors_inread
end module scissors_m
