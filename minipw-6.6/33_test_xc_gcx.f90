INCLUDE 'prepare_all.f90'


PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL test_xc_gcx()
END PROGRAM

! XC functional is read from PWINPUT

!-----------------------
subroutine test_xc_gcx()
!-----------------------
  use xc_gga, only: xc_gcx
  IMPLICIT NONE
  integer :: Nrmesh, Nspin
  real(8), allocatable :: arho(:,:), gradx(:,:,:)
  real(8), allocatable :: sx(:), sc(:)
  real(8), allocatable :: v1x(:,:), v2x(:,:), v1c(:,:), v2c(:,:)


  Nrmesh = 1
  Nspin = 1
  allocate( arho(Nrmesh,Nspin), gradx(3,Nrmesh,Nspin) )
  allocate( sx(Nrmesh), sc(Nrmesh) )
  allocate( v1x(Nrmesh,Nspin), v2x(Nrmesh,Nspin) )
  allocate( v1c(Nrmesh,Nspin), v2c(Nrmesh,Nspin) )
  
  arho(1,1) = 5.d0
  gradx(1,1,1) = 2.d0
  gradx(2,1,1) = 3.d0
  gradx(3,1,1) = 4.d0

  ! we will consider PBE function here as an example

  ! In case of using Libxc, the evaluated sx will also contain
  ! the LDA part times arho + GGA (gradient correction)

  ! In case of using internal xc (not using Libxc), this will only evaluate
  ! the contribution from gradient correction only.

  CALL xc_gcx( Nrmesh, Nspin, arho, gradx, sx, sc, v1x, v2x, v1c, v2c )

  write(*,*) 'sx = ', sx
  write(*,*) 'sc = ', sc
  write(*,*) 'v1x = ', v1x
  write(*,*) 'v2x = ', v2x
  write(*,*) 'v1c = ', v1c
  write(*,*) 'v2c = ', v2c

  deallocate( arho, gradx )
  deallocate( sx, sc )
  deallocate( v1x, v2x, v1c, v2c )

  return
end subroutine