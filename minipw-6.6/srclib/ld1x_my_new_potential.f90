!---------------------------------------------------------------
SUBROUTINE my_new_potential &
     (ndm,mesh,grid,zed,vxt,lsd,nlcc,latt,enne,rhoc,rho,vh,vnew,iflag)
!---------------------------------------------------------------
  ! set up the selfconsistent atomic potential 
  ! from the density and the KS wavefunctions
  !
  use constants, only: fpi, e2
  use radial_grids, only: radial_grid_type, hartree
  use kinds, only : DP
  use funct, only : get_iexch, dft_is_meta, dft_is_gradient
  use ld1inc, only : nwf, vx, vxc, exc, excgga, tau, vtau
  use kli, only : compute_kli_potential
  implicit none
  type(radial_grid_type),intent(in):: grid
  integer, intent(in) :: iflag
  LOGICAL :: nlcc, gga, oep, meta, kli_
  INTEGER :: ndm,mesh,lsd,latt,i,is,nu, nspin, ierr
  real(DP):: rho(ndm,2),vxcp(2),vnew(ndm,2),vxt(ndm),vh(ndm), rhoc(ndm)
  real(DP):: zed,enne,rh(2),rhc, excp
  real(DP),allocatable:: vgc(:,:), egc(:), rhotot(:)
  ! real(DP),allocatable:: vx(:,:)
  real(DP),allocatable:: dchi0(:,:) ! 

  if( mesh /= grid%mesh ) then
    call errore('new_potential', 'mesh dimension is not as expected', 1)
  endif
  gga = dft_is_gradient()
  meta = dft_is_meta()

  oep = get_iexch() == 4
  kli_= get_iexch() == 10

  nspin = 1
  ! this is the case for spinpol calculation
  if( lsd == 1) then
    nspin = 2
  endif
  !
  ! compute hartree potential with the total charge
  !
  allocate(rhotot(ndm), stat=ierr)
  !
  ! get total density
  do i = 1,ndm
     rhotot(i) = rho(i,1)
  enddo
  !
  if(lsd == 1) then
    do i=1,ndm
      rhotot(i) = rhotot(i) + rho(i,2)
    enddo
  endif

  !call hartree(0, 2, mesh, grid, rhotot, vh)
  write(*,*) 'Input for poisson_solve_radial: sum(rhotot) = ', sum(rhotot)
  call poisson_solve_radial(0, 2, mesh, grid, rhotot, vh)
  ! this subroutine is defined in radial_grids:
  ! hartree(k, nst, mesh, grid, f, vh)
  deallocate(rhotot)

  write(*,*) 'sum(vh) (in Ha) = ', sum(vh)*0.5d0

  !
  ! add exchange and correlation potential: LDA or LSDA only
  !
  rhc = 0.0_DP
  do i= 1, mesh
    vh(i) = e2 * vh(i) ! convert to Ry
    do is = 1, nspin
      rh(is) = rho(i,is)/grid%r2(i)/fpi
    enddo
    ! XXX: This is probably relevant
    IF( nlcc ) then
      rhc = rhoc(i)/grid%r2(i)/fpi
    endif
    !
    IF( meta ) THEN
      !
      ! Workaround: the meta-GGA XC functional already contains the LDA part
      !
      vxcp(:) = 0.0_dp
      exc(i) = 0.0_dp
      !print *, "meta gga"
    ELSE
      vxcp = 0
      call vxc_t(lsd, rh, rhc, excp, vxcp)
      exc(i) = excp
    ENDIF
    !
    do is =1, nspin
      vxc(i, is) = vxcp(is)
      vnew(i,is)= -zed*e2/grid%r(i) + vxt(i) + vh(i) + vxcp(is)
    enddo
  endDO
  !
  ! add exchange and correlation potential: GGA only
  !
  if( gga ) then
    allocate(vgc(ndm,2), stat=ierr)
    allocate(egc(ndm), stat=ierr)
    !
    call errore('new_potential', 'allocating vgc and egc', ierr)

    call vxcgc(ndm, mesh, nspin, grid%r, grid%r2, rho, rhoc, &
               vgc, egc, tau, vtau, iflag)
    do is=1,nspin
      do i = 1,mesh
        vxc(i,is) = vxc(i,is) + vgc(i,is)
        vnew(i,is) = vnew(i,is) + vgc(i,is)
        excgga(i) = egc(i)*fpi*grid%r2(i)
      enddo
    enddo
    deallocate(egc)
    deallocate(vgc)
  else
    excgga = 0.0_DP
  endIF
  !
  ! add OEP exchange 
  !
  IF( oep ) THEN
    ! write (*,*) ndm, nwf
    ALLOCATE(dchi0(ndm,nwf))
    DO nu=1,nwf ! num wave functions
      CALL dvex(nu,dchi0(1,nu))
    ENDDO 
    CALL dfx_new(dchi0, vx)
    ! vx contains the oep term
    ! ADD OEP VX
    vnew(:,1:nspin) = vnew(:,1:nspin)  + vx(:,1:nspin)
    DEALLOCATE( dchi0 )
  ENDIF 

  IF( kli_ ) THEN
    CALL compute_kli_potential(grid%mesh,vx)
    vnew(:, 1:nspin) = vnew(:, 1:nspin ) + vx(:,1:nspin)
  ENDIF

  !
  ! latter correction
  !
  IF( latt /= 0) THEN
    DO is = 1,nspin
      DO i = 1,mesh
        vnew(i,is) = min(vnew(i,is),-e2*(zed-enne+1.0_DP)/grid%r(i))
      ENDDO 
    ENDDO 
  ENDIF 

  RETURN 
END SUBROUTINE

