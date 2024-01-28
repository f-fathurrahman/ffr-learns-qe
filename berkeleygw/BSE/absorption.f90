!===============================================================================
!
! Routines:
!
! (1) absorption Originally by JRD Last Edited: 9/12/2011 (JRD)
!
! This routine calls either diag or haydock.
!
! For more details the README_absorption file:
!
! Please report bugs to: jdeslip@civet.berkeley.edu
!
!================================================================================
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
program absorption
  use bse_init_m
  use diag_m
  use global_m
  use haydock_m
  use inread_m
  use references_m
  use timing_m, only: timing => bse_timing
  use write_program_header_m
  implicit none
  type (eqpinfo) :: eqp
  type (xctinfo) :: xct
  type (flags) :: flag
  integer :: nmax,neig,error
  call peinfo_init()
!---------------------------
! Write header
  call write_program_header('BSE/Absorption', .false.)
!---------------------------
! Read absorption.inp
  call logit('Calling inread')
  call open_file(8,file='absorption.inp',form='formatted',status='old')
  call inread(eqp, xct, flag, nmax, neig)
  call close_file(8)
!----------------------
! Initialize HDF5
! FHJ: Initialize xct%nkpt_co and dimensionality of the problem
  call bse_init(xct,flag)
!----------------------
! Initialize random numbers
  peinf%jobtypeeval = 1
!----------------------
! Initialize wcoul0
  xct%wcoul0 = 0d0
!----------------------
! Initialize timer
  call timing%init()
  call timing%start(timing%total)
!---------------------------
! Initialize files
  if (xct%iwritecoul .eq. 1 .and. peinf%inode .eq. 0) then
    call open_file(19,file='vcoul',form='formatted',status='replace')
  endif
  if (xct%algo == BSE_ALGO_HAYDOCK) then
    call haydock(eqp,xct,flag,nmax)
  else
    call diag(eqp,xct,flag,neig,nmax)
  endif
  if (xct%nspinor == 2) then
    call require_reference(REF_Wu2020)
  endif
  call write_memory_usage()
  call show_references()
  call timing%stop(timing%total)
  call timing%print()
end program absorption
