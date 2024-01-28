!==============================================================================
!
! Routines:
!
! (1) Xinvert_with_scalapack_d() Originally by JRD Last Modified 02/2015 (FHJ)
!
! This routine inverts a matrix which is already distributed in block
! cyclic form with ScaLAPACK.
!
! (2) Xinvert_serial() Originally by JRD Last Modified 02/2015 (FHJ)
!
! Inverts a matrix using LAPACK.
!
!==============================================================================
module inversion_m
  use global_m
  use lapack_m
  use scalapack_m
  implicit none
  private
  public :: &
    dinvert_serial, &
    zinvert_serial
contains
!overrules flavor.mk
!===========================================================================
!
! Included from inversion.F90
!
!============================================================================
!---------------------- Use scaLAPACK For Inversion -----------------------------------
!------------------------------------------------------------
subroutine dinvert_serial(nmtx, matrix)
  integer, intent(in) :: nmtx
  real(DP), intent(inout) :: matrix(nmtx,nmtx)
  integer :: ii, info, lwork, ipiv(nmtx)
  real(DP), allocatable :: work(:)
 
  ! FHJ: LU factorization of the matrix
  call dgetrf(nmtx, nmtx, matrix, nmtx, ipiv, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getrf'
    call die('?getrf failed')
  endif
  ! FHJ: tringular inversion of LU decomposition
  allocate(work (10))
  call dgetri(nmtx, matrix, nmtx, ipiv, work, -1, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
    call die('?getri failed for query mode')
  endif
  lwork = max(1,int(work(1)))
  if(allocated(work))then;deallocate(work);endif
  allocate(work (lwork))
  call dgetri(nmtx, matrix, nmtx, ipiv, work, lwork, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
    call die('?getri failed')
  endif
  if(allocated(work))then;deallocate(work);endif
 
  return
end subroutine dinvert_serial


!===========================================================================
!
! Included from inversion.F90
!
!============================================================================
!---------------------- Use scaLAPACK For Inversion -----------------------------------
!------------------------------------------------------------
subroutine zinvert_serial(nmtx, matrix)
  integer, intent(in) :: nmtx
  complex(DPC), intent(inout) :: matrix(nmtx,nmtx)
  integer :: ii, info, lwork, ipiv(nmtx)
  complex(DPC), allocatable :: work(:)
 
  ! FHJ: LU factorization of the matrix
  call zgetrf(nmtx, nmtx, matrix, nmtx, ipiv, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getrf'
    call die('?getrf failed')
  endif
  ! FHJ: tringular inversion of LU decomposition
  allocate(work (10))
  call zgetri(nmtx, matrix, nmtx, ipiv, work, -1, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
    call die('?getri failed for query mode')
  endif
  lwork = max(1,int(work(1)))
  if(allocated(work))then;deallocate(work);endif
  allocate(work (lwork))
  call zgetri(nmtx, matrix, nmtx, ipiv, work, lwork, info)
  if (info/=0) then
    if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
    call die('?getri failed')
  endif
  if(allocated(work))then;deallocate(work);endif
 
  return
end subroutine zinvert_serial
end module inversion_m
