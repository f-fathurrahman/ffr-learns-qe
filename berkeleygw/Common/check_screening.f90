!============================================================================
!
! Routines:
!
! (1) check_screening_trunc Originally by JRD Last Modified: 2/09/2009 (JRD)
!
! Die if screening, truncation, and q0vec are not set consistently.
!
!============================================================================
module check_screening_m
  use global_m
  implicit none
  private
  public :: check_screening_trunc
contains
subroutine check_screening_trunc(itruncflag,iscreen,q0vec,bdot)
  integer, intent(in) :: itruncflag
  integer, intent(in) :: iscreen
  real(DP), intent(in) :: q0vec(3)
  real(DP), intent(in) :: bdot(3,3)
  real(DP) :: q0len
 
  q0len = sqrt(DOT_PRODUCT(q0vec,MATMUL(bdot,q0vec)))
  if (iscreen==SCREEN_METAL .and. q0len<TOL_SMALL) then
    if(peinf%inode == 0) then
      write(0,*) ' '
      write(0,*) 'You want metallic screening but didn''t specify q0vec!!'
    endif
    call die('Inconsistent Screening', only_root_writes = .true.)
  endif
  if (iscreen==SCREEN_GRAPHENE .and. q0len<TOL_SMALL .and. itruncflag==TRUNC_NONE) then
    if(peinf%inode == 0) then
      write(0,*) ' '
      write(0,*) 'You want graphene screening with no truncation'
      write(0,*) 'but didn''t specify q0vec!!'
    endif
    call die('Inconsistent Screening', only_root_writes = .true.)
  endif
  if ((itruncflag==TRUNC_NONE .or. itruncflag==TRUNC_WIRE .or. &
    itruncflag==TRUNC_SLAB) .and. q0len<TOL_SMALL) then
    if(peinf%inode == 0) then
      write(0,*) ' '
      write(0,*) 'You have a divergent Coulomb interaction but didn''t specify q0vec!!'
    endif
    call die('Inconsistent Screening', only_root_writes = .true.)
  endif
  if ((itruncflag/=TRUNC_NONE .and. itruncflag/=TRUNC_WIRE .and. &
    itruncflag/=TRUNC_SLAB) .and. iscreen/=SCREEN_METAL .and. q0len>=TOL_SMALL) then
    if(peinf%inode == 0) then
      write(0,*) ''
      write(0,*) 'You want semiconductor or graphene screening with truncation'
      write(0,*) 'but specified nonzero q0vec!!'
    endif
    call die('Inconsistent Screening', only_root_writes = .true.)
  endif
 
  return
end subroutine check_screening_trunc
end module check_screening_m
