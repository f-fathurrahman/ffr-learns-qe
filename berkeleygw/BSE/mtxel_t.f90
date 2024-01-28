!==============================================================================
!
! Routines:
!
! (1) mtxel_t() Originally By MLT Last Modified 02/13/2013 (FHJ)
!
! Calculated matrix elements needed for interpolation of kernel (dcn,dvn)
!
! input: gvec, wfnc_co, wfnc_fi, wfnv_co, wfnvq_fi types
! igumk index of umklapp, translation, vector
! ikpt index of the current k-point
!
! output: dcc, dvv transformation matrices between a k-point in the
! fine grid and one in the coarse grid
!
!==============================================================================
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
module mtxel_t_m
  use global_m
  implicit none
  public :: mtxel_t
  private
contains
  !> Calculates the transformation matrix between the coarse and fine
  !! wavefunction. Can use full or restricted transformation, where restricted
  !! means cond states are expanded in terms of other cond. states, etc.
  !! TODO: this whole function can be BLAS-ified
  subroutine mtxel_t(gvec,wfnc_co,wfnc_fi,wfnv_co,wfnvq_fi,dcn,dvn,igumk,igumkq,restricted,qflag)
    type (gspace), intent(in) :: gvec
    type (wavefunction), intent(in) :: wfnc_co,wfnc_fi,wfnv_co,wfnvq_fi
    !> Restricted: (wfnc_fi%nband, wfnc_co%nband, wfnc_fi%nspin)
    !! Unrestric.: (wfnc_fi%nband, wfnc_co%nband+wfnv_co%nband, wfnc_fi%nspin)
    real(DP), intent(inout), target :: dcn(:,:,:)
    !> Restricted: (wfnvq_fi%nband, wfnv_co%nband, wfnc_fi%nspin)
    !! Unrestric.: (wfnvq_fi%nband, wfnv_co%nband+wfnc_co%nband, wfnc_fi%nspin)
    real(DP), intent(inout), target :: dvn(:,:,:)
    integer, intent(in) :: igumk,igumkq, qflag
    logical, intent(in) :: restricted
    integer :: ig,ig_co
    integer, allocatable :: isorti(:), isortiq(:)
    real(DP), pointer :: dmm(:,:,:)
   
    dcn(:,:,:)=0.0d0
    dvn(:,:,:)=0.0d0
    ! Compute inverse array to wfnc_co%isort
    allocate(isorti (gvec%ng))
    isorti=0
    do ig=1,gvec%ng
      isorti(wfnc_co%isort(ig))=ig
    enddo
    allocate(isortiq (gvec%ng))
    isortiq=0
    do ig=1,gvec%ng
      isortiq(wfnv_co%isort(ig))=ig
    enddo
    ! Loop over G-vectors
    ! Compute matrix element <c_co,k_co|exp(i(k_co-k_fi).r)|c_fi,k_fi>
    ! fine grid = valence
    ! v_fi <- v_co
    dmm => dvn(:,:wfnv_co%nband,:)
    do ig=1,wfnvq_fi%ng
      call map_ig(gvec, ig, ig_co, igumkq, isortiq, wfnvq_fi, wfnv_co)
      if (ig_co==0) cycle
      call mtxel_dmm(ig, ig_co, wfnvq_fi%nspinor, wfnvq_fi, wfnv_co, dmm)
    enddo
    if (qflag==1 .and. .not. restricted) then
      ! v_fi <- c_co
      dmm => dvn(:,wfnv_co%nband+1:,:)
      do ig=1,wfnvq_fi%ng
        call map_ig(gvec, ig, ig_co, igumkq, isortiq, wfnvq_fi, wfnc_co)
        if (ig_co==0) cycle
        call mtxel_dmm(ig, ig_co, wfnvq_fi%nspinor, wfnvq_fi, wfnc_co, dmm)
      enddo
    endif
    ! fine grid = conduction
    if (restricted) then
      ! c_fi <- c_co
      dmm => dcn(:,:,:)
      do ig=1,wfnc_fi%ng
        call map_ig(gvec, ig, ig_co, igumk, isorti, wfnc_fi, wfnc_co)
        if (ig_co==0) cycle
        call mtxel_dmm(ig, ig_co, wfnc_fi%nspinor, wfnc_fi, wfnc_co, dmm)
      enddo
    else
      ! c_fi <- c_co
      dmm => dcn(:,wfnv_co%nband+1:,:)
      do ig=1,wfnc_fi%ng
        call map_ig(gvec, ig, ig_co, igumk, isorti, wfnc_fi, wfnc_co)
        if (ig_co==0) cycle
        call mtxel_dmm(ig, ig_co, wfnc_fi%nspinor, wfnc_fi, wfnc_co, dmm)
      enddo
      if (qflag==1) then
        ! c_fi <- v_co
        dmm => dcn(:,:wfnv_co%nband,:)
        do ig=1,wfnc_fi%ng
          call map_ig(gvec, ig, ig_co, igumk, isorti, wfnc_fi, wfnv_co)
          if (ig_co==0) cycle
          call mtxel_dmm(ig, ig_co, wfnc_fi%nspinor, wfnc_fi, wfnv_co, dmm)
        enddo
      endif
    endif
    if(allocated(isorti))then;deallocate(isorti);endif
    if(allocated(isortiq))then;deallocate(isortiq);endif
   
    return
  end subroutine mtxel_t
  !> Maps the G vector (ig) from a fine WFN to the coarse one,
  !! including umklapp. Resulting index is saved to ig_co (0 if invalid)
  subroutine map_ig(gvec, ig, ig_co, ig_umklapp, isorti_co, wfn_fi, wfn_co)
    type(gspace), intent(in) :: gvec
    integer, intent(in) :: ig
    integer, intent(out) :: ig_co
    integer, intent(in) :: ig_umklapp
    integer, intent(in) :: isorti_co(:)
    type(wavefunction), intent(in) :: wfn_fi, wfn_co
    integer :: igadd, igaddn, kn(3), kaddn
    ! no push/pop, called too frequently
    ig_co = 0
    if (ig_umklapp==1) then
      ig_co = isorti_co(wfn_fi%isort(ig))
    else
      igadd = wfn_fi%isort(ig)
      ! Compute address of gn=g+gumk:
      ! If the result is beyond planewaves included in wavefunction, skip
      kn(1:3) = gvec%components(1:3,igadd) + &
        gvec%components(1:3,ig_umklapp) + gvec%FFTgrid(1:3)/2 + 1
      if (any(kn(1:3) < 1) .or. any(kn(1:3) > gvec%FFTgrid(1:3))) then
        ig_co = 0
        return
      endif
      kaddn = ((kn(1) - 1)*gvec%FFTgrid(2) + kn(2) - 1)*gvec%FFTgrid(3) + kn(3)
      igaddn = gvec%index_vec(kaddn) ! relate the cube to the sphere
      ig_co = isorti_co(igaddn)
    endif
    if (ig_co > wfn_co%ng) ig_co = 0
  end subroutine map_ig
  !> Calculates the overlap between the coarse and fine wavefunction, store
  !! info the dmm transformation matrix. Assumes you know the mapping between
  !! the the fine g-vector ig and the coarse ig_co
  subroutine mtxel_dmm(ig, ig_co, nspinor, wfn_fi, wfn_co, dmm)
    integer, intent(in) :: ig, ig_co
    integer, intent(in) :: nspinor
    type(wavefunction), intent(in) :: wfn_fi, wfn_co
    real(DP), intent(inout) :: dmm(:,:,:) !< (wfn_fi%nband, wfn_co%nband, nspin)
    integer :: im_fi, im_co, is, ispinor
    ! no push/pop, called too frequently
    do is=1,wfn_co%nspin
      do im_co=1,wfn_co%nband
        do im_fi=1,wfn_fi%nband
          do ispinor=1,nspinor
            dmm(im_fi, im_co, is) = dmm(im_fi, im_co, is) + &
              wfn_fi%cg(ig, im_fi, is*ispinor) * (wfn_co%cg(ig_co, im_co, is*ispinor))
          enddo
        enddo
      enddo
    enddo
  end subroutine mtxel_dmm
end module mtxel_t_m
