module exposed_ld1x

! to "expose" fields of objects of derived types 
! "Exposed variables" are the global variables defined in this module
! Call expose_* subroutines before accessing them. These subroutines
! simply copy the relevant field to global variables.


use kinds, only: dp
!
implicit none
!
! scalars
integer :: grid_mesh  ! actual number of grid points
real(dp) :: grid_rmax  ! the maximum radial point
real(dp) :: grid_zmesh  ! the ionic charge used for the mesh
real(dp) :: grid_dx      ! the deltax of the linear mesh
!
! arrays
real(dp), allocatable :: grid_r(:) ! radial mesh
real(dp), allocatable :: grid_rab(:) ! dr(x) / dx where x is the linear grid


!-----------------------------------------
  contains
!-----------------------------------------
! copy to global variable


!----------------------------
subroutine expose_grid_mesh()
!----------------------------
  use ld1inc, only: grid
  implicit none
  grid_mesh = grid%mesh
  return
end subroutine

!----------------------------
subroutine expose_grid_rmax()
!----------------------------
  use ld1inc, only: grid
  implicit none
  grid_rmax = grid%rmax
  return
end subroutine

!----------------------------
subroutine expose_grid_zmesh()
!----------------------------
  use ld1inc, only: grid
  implicit none
  grid_zmesh = grid%zmesh
  return
end subroutine


!--------------------------
subroutine expose_grid_dx()
!--------------------------
  use ld1inc, only: grid
  implicit none
  grid_dx = grid%dx
  return
end subroutine




!---------------------------
subroutine expose_grid_r()
!---------------------------
  use ld1inc, only: grid
  implicit none
  !
  if(allocated(grid_r)) deallocate(grid_r)
  allocate( grid_r(grid%mesh) )
  !
  grid_r(:) = grid%r(:) ! copy
  !
  return
end subroutine



!---------------------------
subroutine expose_grid_rab()
!---------------------------
  use ld1inc, only: grid
  implicit none
  !
  if(allocated(grid_rab)) deallocate(grid_rab)
  allocate( grid_rab(grid%mesh) )
  !
  grid_rab(:) = grid%rab(:) ! copy
  !
  return
end subroutine




end module