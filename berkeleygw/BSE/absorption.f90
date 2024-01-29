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
