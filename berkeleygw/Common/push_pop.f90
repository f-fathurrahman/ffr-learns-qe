!==========================================================================
!
! Module push_pop_m Originally By DAS
!
! Create a stack trace of routines entered and exited, for debugging.
! This only takes effect if the code is compiled with -DDEBUG.
! Enable by setting 'debug_level' below, and recompile:
! 0: no debugging trace
! 1: only node 0 writes trace
! 2: all nodes write trace. Very slow.
! Inspired by Octopus messages.F90 (originally revision 6920)
!
!==========================================================================

module push_pop_m
  use message_m
  use nrtype_m
  use peinfo_m
  implicit none
end module push_pop_m
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
