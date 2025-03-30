!
! Solution of the Poisson's equation on a radial (logarithmic) grid
!
! debug version of subroutine hartree in radial_grids.f90
!---------------------------------------------------------------
subroutine poisson_solve_radial(k, nst, mesh, grid, f, vh)
!---------------------------------------------------------------
  !
  use upf_kinds, only : DP
  use radial_grids, only: radial_grid_type, series
  implicit none
  integer,intent(in)::       & 
       k,   & ! input: the k of the equation
       nst, & ! input: at low r, f goes as r**nst
       mesh   ! input: the dimension of the mesh

  type(radial_grid_type), intent(in) :: &
       grid   ! input: the radial grid
  real(DP), intent(in)::        &
       f(mesh)  ! input: the 4\pi r2 \rho function
  real(DP), intent(out)::       &
       vh(mesh) ! output: the required solution
  !
  ! local variables
  !
  integer ::        &
       k21,  &   ! 2k+1
       nk1,  &   ! nst-k-1
       ierr, &   ! integer variable for allocation control
       i         ! counter

  real(DP)::        &
       c0,c2,c3, & ! coefficients of the polynomial expansion close to r=0
       ch,       & ! dx squared / 12.0
       xkh2,     & ! ch * f
       ei, di,   & ! auxiliary variables for the diagonal and 
                   ! off diagonal elements of the matrix
       f1, fn,   & ! variables used for the boundary condition
       vhim1, vhi  ! variables for the right hand side

  real(DP), allocatable:: &
       d(:), &       ! the diagonal elements of 
                     ! the tridiagonal sys.
       e(:)          ! the off diagonal elements 
                     ! of the trid. sys.
  
  write(*,*) 'Enter poisson_solve_radial'

  !
  ! Allocate space for the diagonal and off diagonal elements
  !
  if(mesh /= grid%mesh) then
    call upf_error('hartree',' grid dimension mismatch',1) 
  endif
  allocate(d(mesh),stat=ierr)
  allocate(e(mesh),stat=ierr)

  if(ierr /= 0) then
    call upf_error('hartree',' error allocating d or e',1)
  endif
  !
  ! Find the series expansion of the solution close to r=0
  !
  k21 = 2*k + 1
  nk1 = nst - k - 1
  write(*,*) 'poisson_solve_radial: mesh = ', mesh
  write(*,*) 'poisson_solve_radial: k21 = ', k21
  write(*,*) 'poisson_solve_radial: nk1 = ', nk1
  if(nk1 <= 0) then
    write(6,100) k,nst
    100 format(5x,'stop in "hartree": k=',i3,'  nst=',i3)
    stop
    !else if(nk1.ge.4) then
    ! not sure whether the following is really correct, but the above wasn't
  elseif( nk1 >= 3) then
    c2 = 0.d0
    c3 = 0.d0
  else
    e(1) = 0.d0
    do i = 1,4
      d(i) = -k21*f(i) / grid%r(i)**nst
    enddo
    call series(d, grid%r, grid%r2, e(nk1))
    c2 = e(1) / (4.d0*k + 6.d0)
    c3 = e(2) / (6.d0*k + 12.d0)
  endif
  write(*,*) 'After series: e(1:4) = ', e(1:4)
  write(*,*) 'c2 = ', c2
  write(*,*) 'c3 = ', c3
  !
  ! Set the main auxiliary parameters
  !
  ch = grid%dx * grid%dx / 12.d0
  xkh2 = ch*(DBLE(k) + 0.5d0)**2
  ei = 1.d0  - xkh2
  di = -(2.d0 + 10.d0*xkh2)
  write(*,*) 'poisson_solve_radial: ch = ', ch
  write(*,*) 'poisson_solve_radial: xkh2 = ', xkh2
  write(*,*) 'poisson_solve_radial: ei = ', ei
  write(*,*) 'poisson_solve_radial: di = ', di
  !
  ! Set the diagonal and the off diagonal elements of the 
  ! linear system, compute a part of the right hand side 
  !
  do i = 2,mesh
    d(i) = -di
    e(i) = -ei
    vh(i) = k21*ch*grid%sqr(i) * f(i)
  enddo
  !
  ! Use the boundary condition to eliminate the value of the 
  ! solution in the first point from the first equation. This 
  ! part for the diagonal element
  !
  f1 = ( grid%sqr(1)/grid%sqr(2) )**k21
  write(*,*) 'poisson_solve_radial: f1 = ', f1
  d(2) = d(2) - ei*f1
  !
  ! Use the boundary condition to eliminate the value of the 
  ! solution in the last point from the last equation
  !
  fn = ( grid%sqr(mesh-1)/grid%sqr(mesh) )**k21
  write(*,*) 'poisson_solve_radial: fn = ', fn
  d(mesh-1) = d(mesh-1) - ei*fn
  !
  ! In the first point vh(1) has the same definition as in the other points
  !
  vhim1 = k21 * ch * grid%sqr(1) * f(1)
  write(*,*) 'poisson_solve_radial: vhim1 = ', vhim1
  !
  ! Compute the right hand side using the auxiliary quantity vh(i).
  !
  do i = 2,mesh-1
     vhi = vh(i)
     vh(i) = vhim1 + 10.d0*vhi + vh(i+1)
     vhim1 = vhi
  enddo
  !
  ! Use the boundary condition to eliminate the value of the solution in the 
  ! first point from the first equation. This part for the right hand side.
  !
  vh(2) = vh(2) - ei*grid%sqr(1)**k21 * ( &
            c2*( grid%r2(2) - grid%r2(1) ) + &
            c3*( grid%r(2)**3 - grid%r(1)**3 ) &
          )
  write(*,*) 'sum d(2:mesh-1) = ', sum(d(2:mesh-1))
  write(*,*) 'sum e(2:mesh-1) = ', sum(e(2:mesh-1))
  write(1111,*) d
  write(1112,*) e
  write(1113,*) vh
  !
  ! solve the linear system with lapack routine dptsv
  !
  call dptsv(mesh-2, 1, d(2), e(2), vh(2), mesh-2, ierr)
  if( ierr /= 0 ) then
    call upf_error('poisson_solve_radial', 'error in lapack', ierr)
  endif
  !
  ! Set the value of the solution at the first and last point
  ! First, find c0 from the solution in the second point
  !
  c0 = vh(2)/grid%sqr(2)**k21 - c2*grid%r2(2) - c3*grid%r(2)*grid%r2(2)
  !
  ! and then use the series expansion at the first point
  !
  vh(1) = grid%sqr(1)**k21 * ( c0 + c2*grid%r2(1) + c3*grid%r(1)**3 )
  !
  ! the solution at the last point is given  by the boundary 
  ! condition
  !
  vh(mesh) = vh(mesh-1)*fn
  !
  ! The solution must be divided by r (from the equation) 
  ! and multiplied by the square root of r (from the log 
  ! mesh transformation)
  !
  do i = 1,mesh
     vh(i) = vh(i) / grid%sqr(i)
  enddo
  ! final vh
  write(1114,*) vh

  deallocate(e)
  deallocate(d)

  return
end subroutine
