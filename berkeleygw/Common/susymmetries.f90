!============================================================================
!
! Routines:
!
! (1) susymmetries Originally by BAB Last Modified: 4/01/2014 (BAB)
!
! Find the 2x2 matrix to rotate spinor components of wavefunctions,
! given the integer rotation matrix in lattice coordinates.
! This function must be called in all genwf files for spinor calculations.
!
!============================================================================
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
module susymmetries_m
  use global_m
  use random_m
  use misc_m
  use blas_m
  implicit none
  private
  public :: susymmetries, rot_axis_angle, su2_mat, random_rotation
contains
  !> This routine reads in a an SO(3) rotation matrix associated with a space group symmetry
  !! and generates the corresponding rotation matrix in the SU(2) group (in the postive branch).
  !! This matrix in SU(2), umtrx, is then used to rotate the spinors of the wavefunction.
  !! This routine is an essential part of all genwf routines, which use symmetry operations to
  !! reconstruct the wavefunctions in the little group of non-zero k- or q-vectors.
  subroutine susymmetries(bvec, mtrx, umtrx, itqq)
    real(DP), intent(in) :: bvec(3,3) !< crystal%bvec
    integer, intent(in) :: mtrx(3,3) !< rotation matrices
    complex(DPC), intent(out) :: umtrx(2,2) !< spinor rotation matrix
    integer, intent(in) :: itqq !< symmetry index
    real(DP) :: axis(3) !< rotation axis
    real(DP) :: angle, bvecinv(3,3)
    real(DP) :: mtrxtemp(3,3),det
   
    ! NOTE: mtrx usually has a third index, indicating symmetry type
    ! Check if identity
    if (itqq.eq.1) then
      umtrx(1,1) = (1.0d0,0.0d0)
      umtrx(1,2) = (0.0d0,0.0d0)
      umtrx(2,1) = (0.0d0,0.0d0)
      umtrx(2,2) = (1.0d0,0.0d0)
     
      return
    endif
    ! calculate determinant of matrix mtrx
    ! for umtrx, improper rotations should be turned into proper rotations
    ! mtrxtemp = TRANSPOSE(mtrx)
    mtrxtemp = mtrx
    call compute_det(mtrxtemp,det)
    if (det.lt.TOL_Zero) then
      mtrxtemp = -mtrxtemp
    endif
    ! convert rotation matrix for lattice vectors to cartesian coordinates
    ! R_cart = B_Transpose * R * (B_Transpose)_Inverse
    ! Note: bvec is already transposed
    call invert_matrix(bvec, bvecinv)
    mtrxtemp = MATMUL(bvec,(MATMUL(mtrxtemp,bvecinv)))
    ! get rotation axis and rotation angle
    call rot_axis_angle(axis,angle,mtrxtemp)
    call su2_mat(axis,angle,umtrx)
   
    return
  end subroutine susymmetries
!===============================================================================
  !> Given a 3x3 orthogonal matrix, this routine gives the rotation axis and angle
  !!
  !! This algorithm is based on the publication
  !! F. Landis Markley, Journal of Guidance, Control, and Dynamics 31 (2008), p. 440
  subroutine rot_axis_angle(axis, angle, mtrx)
    real(DP), intent(in) :: mtrx(3,3) !< rotation matrices
    real(DP), intent(out) :: axis(3) !< rotation axis
    real(DP), intent(out) :: angle
    integer :: ind, ii, jj, indmax, cyclic(3)
    real(DP) :: tr, vecnorm, maxnorm, quat(4), xquat(4,4), normx(4)
   
    quat(:) = 0.0d0
    ! Ordering for quaternion:
    ! (axis(1) sin(angle/2), axis(2) sin(angle/2), axis(3) sin(angle/2), cos(angle/2))
    tr = 0.0d0
    do ind=1,3
      tr = tr + mtrx(ind,ind)
    enddo
    cyclic(1) = 2
    cyclic(2) = 3
    cyclic(3) = 1
    do ii=1,4
      do jj=1,4
        if (ii.eq.jj .and. jj.lt.4) then
          xquat(ii,jj) = 1 + 2*mtrx(ii,jj) - tr
        else if (ii.ne.jj .and. ii.lt.4 .and. jj.lt.4) then
          xquat(ii,jj) = mtrx(jj,ii) + mtrx(ii,jj)
        else if (ii.eq.jj .and. jj.eq.4) then
          xquat(ii,jj) = 1 + tr
        else if (ii.ne.jj .and. jj.eq.4) then
          xquat(ii,jj) = mtrx(6-cyclic(ii)-ii,cyclic(ii)) - mtrx(cyclic(ii),6-cyclic(ii)-ii)
        endif
      enddo
    enddo
    do ii=1,4
      do jj=ii+1,4
        xquat(jj,ii) = xquat(ii,jj);
      enddo
    enddo
    normx(:) = 0.0d0
    do ii=1,4
      do jj=1,4
        normx(ii) = normx(ii) + xquat(ii,jj)**2;
      enddo
      normx(ii) = sqrt(normx(ii));
    enddo
    ! Note that maxnorm is guaranteed to be at least 1.0d0
    maxnorm = 0.0d0
    indmax = 0
    do ii=1,4
      if (normx(ii) .gt. maxnorm) then
        maxnorm = normx(ii)
        indmax = ii
      endif
    enddo
    quat(:) = xquat(:,indmax)/maxnorm
    vecnorm = 0.0d0
    do ii=1,3
      vecnorm = vecnorm + quat(ii)**2
    enddo
    vecnorm = sqrt(vecnorm)
    angle = 2*atan2(vecnorm,quat(4))
    if (vecnorm.eq.0) then
      ! probably does not happen since unity matrix is treated separately
      axis(1) = 0.0d0
      axis(2) = 0.0d0
      axis(3) = 1.0d0
    else
      axis(:) = quat(1:3)/vecnorm
    endif
   
    return
  end subroutine rot_axis_angle
!===============================================================================
  !> This routine takes in an axis and an angle from an SO(3) matrix and generates
  !! the corresponding SU(2) matrix.
  subroutine su2_mat(axis,angle,umtrx)
    real(DP), intent(in) :: axis(3) !< rotation axis
    real(DP), intent(in) :: angle
    complex(DPC), intent(out) :: umtrx(2,2) !< SU(2) rotation matrices
    real(DP) :: cosa, sina
   
    cosa = COS(angle*0.5d0)
    sina = SIN(angle*0.5d0)
    ! Note: this is actually the transpose of the usual U matrix, for use with genwf
    umtrx(1,1) = cmplx(cosa,-axis(3)*sina,kind=DPC)
    umtrx(2,1) = cmplx(-axis(2)*sina,-axis(1)*sina,kind=DPC)
    umtrx(1,2) = -conjg(umtrx(2,1))
    umtrx(2,2) = conjg(umtrx(1,1))
   
    return
  end subroutine su2_mat
!===============================================================================
  !> This routine is for use in wfnmix_spinor. It creates a random rotation axis and angle
  !! and constructs the associated random rotation matrix.
  !!
  !! As implemented in wfnmix_spinor, we start off with the non-spinor wfn coeffs, c(G).
  !! We then create pairs of bands, c1(G) = c(G)*alpha and c2(G) = c(G)*beta, with
  !! alpha and beta the eigenvectors of Paul matrix sigma_z.
  !! We then rotate c1 and c2 by a random angle about a random axis according to
  !! cos(theta_rand / 2.0d0) - i * sigma-dot-axis_rand sin(theta_rand / 2.0d0).
  !!
  !! One tricky aspect is degenerate subspaces (from the spin-less calculation).
  !! These have to have the same axis and angle, as the coefficients were already
  !! carefully pre-orthogonalized. Having degenerate bands have different spin-space
  !! character would lead to artificial overlaps, generally.
  subroutine random_rotation(umtrx)
    complex(DPC), intent(out) :: umtrx(2,2) !< SU(2) rotation matrices
    real(DP) :: rangle, rnd(4), axis_norm
    real(DP) :: raxis(3)
    integer :: ii
   
    ! generate random numbers between 0 and 1
    call genrand_init(put=5000)
    do ii=1,4
      call genrand_real4(rnd(ii))
    enddo
    ! create random axis and normalize
    raxis=(/rnd(1), rnd(2), rnd(3)/)
    axis_norm=blas_nrm2(3, raxis, 1)
    if (axis_norm .gt. TOL_ZERO) then
      raxis = raxis / axis_norm
    else
      rangle = 0.0d0
      raxis = (/0.0d0,0.0d0,1.0d0/)
    endif
    ! create random angle from rnda, from -pi to pi
    rangle=2.0d0*PI_D*(rnd(4)-0.5)
    call su2_mat(raxis, rangle, umtrx)
   
    return
  end subroutine random_rotation
end module susymmetries_m
