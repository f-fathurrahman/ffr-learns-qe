!===========================================================================
!
! Modules:
!
! scalapack_m Originally By DAS
!
! Functions, types, and interfaces for ScaLAPACK/BLACS.
! Interfaces are from http://www.netlib.org/scalapack/tools, double, complex16
! and from http://www.netlib.org/blacs/BLACS/QRef.html (entered manually...)
! Every ScaLAPACK/BLACS function used in the code should be listed here, and this
! module should be used in every routine containing ScaLAPACK/BLACS calls to ensure
! the argument types are correct.
!
!============================================================================

module scalapack_m
  use global_m
  implicit none
  private
  public :: &
    scalapack, &
    blacs_setup, &
    layout_scalapack, &
    blacs_pinfo, &
    iceil, &
    descinit, &
    descset, &
    pdgesv, &
    pzgesv, &
    pdsyevx, &
    pdgeqrf, &
    pzgeqrf, &
    pzheevx, &
    pdgemr2d, &
    pzgemr2d, &
    blacs_get, &
    blacs_gridinit, &
    blacs_gridmap, &
    blacs_gridexit, &
    blacs_exit
!-----------------------------
  type scalapack
    integer :: nprow !< the number of processors in a row of your processor grid
    integer :: npcol !< the number of processors in a column of your processor grid
    integer :: nbl !< the linear dimension of a block of a distributed matrix
    integer :: myprow !< the processor`s row coordinate in your processor grid
    integer :: mypcol !< the processor`s column coordinate in your processor grid
    integer :: npr !< the number of rows of the matrix the processor owns
    integer :: npc !< the number of columns of the matrix the processor owns
    integer :: icntxt !< BLACS context; see BLACS documentation
    integer, pointer :: npcd(:) !< global list of the number of cols of the matrix owned by all processors
    integer, pointer :: nprd(:) !< global list of the number of rows of the matrix owned by all processors
    integer, pointer :: isrtxrow(:) !< isrtxrow/isrtxcol give the sorted index of the gvector in a given block
    integer, pointer :: isrtxcol(:) !! owned by a processor in terms of the whole list of gvectors
    integer, pointer :: imycol(:) !< imyrow/imycol give the row/column index of a g-vector owned by a given
    integer, pointer :: imyrow(:) !! processor in the whole matrix
    integer, pointer :: imycolinv(:) !< inverse of imycol
    integer, pointer :: imyrowinv(:) !! inverse of imyrow
    integer, pointer :: imycold(:,:) !< imycold/imyrowd are global lists of the row/column index of g-vectors
    integer, pointer :: imyrowd(:,:) !! owned by all the processors in the whole matrix
  end type scalapack
!> SCALAPACK
  interface
    subroutine BLACS_PINFO( MYPNUM, NPROCS )
      implicit none
      INTEGER MYPNUM, NPROCS
    end subroutine
  end interface
  interface
    INTEGER FUNCTION ICEIL( INUM, IDENOM )
      implicit none
      INTEGER IDENOM, INUM
    end FUNCTION ICEIL
  end interface
  interface
    SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO )
      implicit none
      INTEGER ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
      INTEGER DESC( * )
    end SUBROUTINE DESCINIT
  end interface
  interface
    SUBROUTINE DESCSET( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD )
      implicit none
      INTEGER ICSRC, ICTXT, IRSRC, LLD, M, MB, N, NB
      INTEGER DESC( * )
    end SUBROUTINE DESCSET
  end interface
  interface
    SUBROUTINE PDGESV( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
      implicit none
      INTEGER IA, IB, INFO, JA, JB, N, NRHS
      INTEGER DESCA( * ), DESCB( * ), IPIV( * )
      DOUBLE PRECISION A( * ), B( * )
    end SUBROUTINE PDGESV
  end interface
  interface
    SUBROUTINE PZGESV( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB, DESCB, INFO )
      implicit none
      INTEGER IA, IB, INFO, JA, JB, N, NRHS
      INTEGER DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX*16 A( * ), B( * )
    end SUBROUTINE PZGESV
  end interface
  interface
    SUBROUTINE PDSYEVX( JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL, &
      VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, Z, IZ, JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, IFAIL, &
      ICLUSTR, GAP, INFO )
      implicit none
      CHARACTER JOBZ, RANGE, UPLO
      INTEGER IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LWORK, M, N, NZ
      DOUBLE PRECISION ABSTOL, ORFAC, VL, VU
      INTEGER DESCA( * ), DESCZ( * ), ICLUSTR( * ), IFAIL( * ), IWORK( * )
      DOUBLE PRECISION A( * ), GAP( * ), W( * ), WORK( * ), Z( * )
    end SUBROUTINE PDSYEVX
  end interface
  interface
    SUBROUTINE PDGEQRF( M, N, A, IA, JA, DESCA, TAU, WORK, LWORK, INFO )
      implicit none
      INTEGER IA, INFO, JA, LWORK, M, N
      INTEGER DESCA( * )
      DOUBLE PRECISION A( * ), TAU( * ), WORK( * )
    end SUBROUTINE PDGEQRF
  end interface
  interface
    SUBROUTINE PZGEQRF( M, N, A, IA, JA, DESCA, TAU, WORK, LWORK, INFO )
      implicit none
      INTEGER IA, INFO, JA, LWORK, M, N
      INTEGER DESCA( * )
      COMPLEX*16 A( * ), TAU( * ), WORK( * )
    end SUBROUTINE PZGEQRF
  end interface
  interface
    subroutine PZHEEVX( JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL, &
      VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, Z, IZ, &
      JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, &
      LIWORK, IFAIL, ICLUSTR, GAP, INFO )
      implicit none
      character JOBZ, RANGE, UPLO
      integer IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LRWORK, LWORK, M, N, NZ
      double precision ABSTOL, ORFAC, VL, VU
      integer DESCA( * ), DESCZ( * ), ICLUSTR( * ), IFAIL( * ), IWORK( * )
      double precision GAP( * ), RWORK( * ), W( * )
      complex*16 A( * ), WORK( * ), Z( * )
    end subroutine PZHEEVX
  end interface
  interface
    subroutine PDGEMR2D( M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICNTXT)
      implicit none
      integer M, N, IA, JA, IB, JB, DESCA( * ), DESCB( * ), ICNTXT
      double precision A( * ), B( * )
    end subroutine PDGEMR2D
  end interface
  interface
    subroutine PZGEMR2D( M, N, A, IA, JA, DESCA, B, IB, JB, DESCB, ICNTXT)
      implicit none
      integer M, N, IA, JA, IB, JB, DESCA( * ), DESCB( * ), ICNTXT
      double complex A( * ), B( * )
    end subroutine PZGEMR2D
  end interface
!> BLACS
  interface
    subroutine blacs_get(icontxt, what, val)
      implicit none
      integer, intent(in) :: icontxt
      integer, intent(in) :: what
      integer, intent(out) :: val
    end subroutine blacs_get
  end interface
  interface
    subroutine blacs_gridinit(icontxt, order, nprow, npcol)
      implicit none
      integer, intent(inout) :: icontxt
      character, intent(in) :: order
      integer, intent(in) :: nprow
      integer, intent(in) :: npcol
    end subroutine blacs_gridinit
  end interface
  !> note: args are out of order so that ldumap,npcol are declared
  !! prior to their usage as dimensions of usermap.
  interface
    subroutine blacs_gridmap(icontxt, usermap, ldumap, nprow, npcol)
      implicit none
      integer, intent(inout) :: icontxt
      integer, intent(in) :: ldumap
      integer, intent(in) :: nprow
      integer, intent(in) :: npcol
      integer, intent(in) :: usermap(ldumap,npcol)
    end subroutine blacs_gridmap
  end interface
  interface
    subroutine blacs_gridexit(icontxt)
      implicit none
      integer, intent(in) :: icontxt
    end subroutine blacs_gridexit
  end interface
  interface
    subroutine blacs_exit(icontxt)
      implicit none
      integer, intent(in) :: icontxt
    end subroutine blacs_exit
  end interface
  interface
    subroutine blacs_gridinfo(icontxt, nprow, npcol, myprow, mypcol)
      implicit none
      integer, intent(in) :: icontxt
      integer, intent(out) :: nprow
      integer, intent(out) :: npcol
      integer, intent(out) :: myprow
      integer, intent(out) :: mypcol
    end subroutine blacs_gridinfo
  end interface
contains
!>--------------------------------------------------------------------------
!! Originally by AC, last modified 6/12/2008 (JRD)
!! Figures out a p by q processor grid layout for the scalapack library.
!! This p by q grid is used to partition the matrix with a block size b.
!! The goal is to get a processor grid which is as close to "square" as
!! possible. For more details, see scalapack documentation.
!!
!! Input nproc number of processors
!! matsize size of matrix
!!
!! Output nbl block size
!! nprow processor grid row
!! npcol processor grid column
  subroutine layout_scalapack(matsize, nbl, nproc, nprow, npcol)
    integer, intent(in) :: matsize
    integer, intent(out) :: nbl
    integer, intent(in) :: nproc
    integer, intent(out) :: nprow, npcol
    integer :: i
   
!------------------
! Find processor grid
    nprow = int(sqrt(dble(nproc) + 1.0d-6))
    do i = nprow, 1, -1
      if(mod(nproc, i) .eq. 0) exit
    enddo
    nprow = i
    npcol = nproc/nprow
!-------------------
! Now for the block size
    nbl = min(32, matsize/(max(nprow, npcol)))
!-------------------
! Ensure nonzero
    nbl = max(nbl, 1)
   
    return
  end subroutine layout_scalapack
!--------------------------------------------------------------------------
  subroutine blacs_setup(scal, matsize, is_row_order,nppgroup_f,nfreq_group,np_left)
    type(scalapack), intent(inout) :: scal !< other elements might have been set earlier
    integer, intent(in) :: matsize
    logical, intent(in) :: is_row_order
    integer, intent(in),optional :: nppgroup_f !< # of proc. per freq. group
    integer, intent(in),optional :: nfreq_group !< number of parallel frequencies
    integer, intent(in),optional :: np_left !< # of proc. leftover for parallel frequencies
    character :: order
    integer :: iw,ir,ic,npe
    logical :: custom_grid
    integer,allocatable :: usermap(:,:)
   
    scal%npr=matsize
    scal%npc=matsize
    scal%nbl=matsize
    scal%nprow=1
    scal%npcol=1
    scal%myprow=0
    scal%mypcol=0
   
    return
  end subroutine blacs_setup
end module scalapack_m
