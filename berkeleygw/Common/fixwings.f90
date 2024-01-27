!============================================================================
!
! MODULE: fixwings_m
!
!> Rescale epsmat to make it compatible with W averaging
!
! DESCRIPTION:
!> This module contains routines to rescale the wings and \b head of the
!! epsilon matrix so that later we can compute <W> = <epsinv(q) * v(q)>.
!! Note that <W> should take into consideration the analytical form of epsmat(q)
!! and v(q) for small q, for each type of truncation/screening.
!
! REVISION HISTORY:
! 15 Jan 2009 - Initial version (JRD and GSM)
! 13 Oct 2011 - Renamed internal variables (FHJ)
!
! (1) fixwings() Originally by JRD Last Modified: 2/09/2009 (JRD)
!
! (2) fixwings_dyn() Originally by GSM Last Modified: 15/09/2009 (JRD)
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
module fixwings_m
  use global_m
  implicit none
  private
  public :: &
    fixwings, &
    fixwings_dyn
contains
  !> Despite its name, this routine fixes both the wings and \b head of the
  !! GPP epsmat. The goal is to rescale epsmat so that we get the appropriate
  !! W averaging, i.e., \f$ W_0 = \varepsilon^{-1}(q) v(q) \f$.
  !!
  !! \param vcoul Bare potential (v) for the given q point, taking into
  !! consideration the truncation. For 3D SC, this should be 8*PI/q^2
  !! \param wcoul0 Screened potential W_0. This is actually the value that W_0
  !! should have, and \f$\varepsilon^{-1}(q)\f$ is rescaled in this function so
  !! that \f$W_0 = \varepsilon^{-1}(q) v(q)\f$
  !! \param epstemp A specific column (#icol::icol) of the dielectric matrix,
  !! for a specific q point
  !! \param icutv Truncation
  !! \param iscreen Screening type
  !! \param icol Column of the whole dielectric constant that we are dealing with
  !! \param nmtx Number of cols/rows in the epsmat
  !! \param irow_G0 Which row of epstmp represents the G=0 vector?
  !! \param oneoverq WARNING: This is actually 8*PI/q!!!
  !! \param q0flag Is this the q0 point?
  !! \param averagew Should we do W averaging?
  !! \param bdot Reciprocal metric
  !!
  !! \sa fixwings_dyn - the FF version
  subroutine fixwings(vcoul,wcoul0,epstemp,icutv,iscreen, &
    icol,nmtx,irow_G0,q0len,oneoverq,fact,q0flag,averagew,bdot)
    integer, intent(in) :: icol,nmtx,irow_G0,iscreen,icutv
    real(DP), intent(in) :: vcoul,oneoverq,q0len,fact
    real(DP), intent(inout) :: epstemp(nmtx)
    real(DP), intent(in) :: wcoul0
    logical, intent(in) :: q0flag,averagew
    real(DP), intent(in) :: bdot(3,3)
    integer :: i
    real(DP) :: zc
!if irow_G0 < avgcut \=icol
!epstemp()*oneverq*q0len
!if icol < avgcut i\=icol
!epstemp(i)*oneoverq/(vcoul(icol)*qicol00)
!----------------------
! No Truncation
   
    if (icutv==TRUNC_NONE) then
      if (icol.ne.irow_G0) then ! wing` (Gp/=0)
        if (iscreen==SCREEN_SEMICOND) then
          epstemp(irow_G0) = epstemp(irow_G0)*oneoverq*q0len/(8D0*PI_D)
        endif
        ! JRD Zero out q0 wings
        if (q0flag .and. iscreen==SCREEN_SEMICOND) epstemp(irow_G0) = 0D0
      else
        do i=1,nmtx
          if (i .ne. irow_G0) then ! wing (G/=0)
            if (iscreen==SCREEN_SEMICOND) then
              epstemp(i) = epstemp(i)*fact*oneoverq/(vcoul*q0len)
            endif
            if (iscreen==SCREEN_GRAPHENE) then
              epstemp(i) = epstemp(i)*fact*8D0*PI_D/(vcoul*q0len**2)
            endif
            ! JRD Zero out q0 wings
            if (q0flag .and. iscreen==SCREEN_SEMICOND) epstemp(i) = 0D0
          else ! Head
            if (q0flag .and. averagew) then
               epstemp(i) = wcoul0/vcoul
            endif
          endif
        enddo
      endif
    endif
!----------------------
! Cell Wire Truncation
! May not be implemented correctly for graphene screening... I`m not even sure there is a 1D system with linear DOS...
    if (icutv==TRUNC_WIRE) then
      if (icol.ne.irow_G0) then ! wing` (Gp/=0)
        ! JRD We zero q0 wings
        if (q0flag .and. iscreen==SCREEN_SEMICOND) then
          epstemp(irow_G0) = 0d0
        endif
      else
        do i=1,nmtx
          if (i .ne. irow_G0) then ! wing (G/=0)
            ! JRD We zero q0 wings
            if (q0flag .and. iscreen==SCREEN_SEMICOND) then
              epstemp(i) = 0d0
            endif
          else ! Head
            if (q0flag .and. averagew .and. iscreen/=SCREEN_METAL) then
              epstemp(i) = wcoul0/vcoul
            endif
          endif
        enddo
      endif
    endif
!----------------------
! Cell Slab Truncation
    if (icutv==TRUNC_SLAB) then
      zc=2D0*PI_D/(sqrt(bdot(3,3))*2D0)
      if (icol.ne.irow_G0) then ! wing` (Gp/=0)
        ! JRD We zero q0 wings
        if (q0flag .and. iscreen==SCREEN_SEMICOND) epstemp(irow_G0) = 0d0
      else
        do i=1,nmtx
          if (i .ne. irow_G0) then ! wing (G/=0)
            if (iscreen/=SCREEN_METAL) then
              epstemp(i) = epstemp(i) * 8D0 * PI_D * fact * zc &
                / (vcoul * q0len)
              ! JRD We zero q0 wings
              if (q0flag .and. iscreen==SCREEN_SEMICOND) epstemp(i) = 0d0
            endif
          else ! Head
            if (q0flag .and. averagew .and. iscreen/=SCREEN_METAL) then
               epstemp(i) = wcoul0/vcoul
            endif
          endif
        enddo
      endif
    endif
   
    return
  end subroutine fixwings
!============================================================================
  !> Full frequency version of #fixwings
  !!
  !! \see fixwings - the GPP version
  subroutine fixwings_dyn(vcoul,epstemp,icutv,iscreen,icol, &
    nfreq,nmtx,irow_G0,q0len,oneoverq,fact,q0flag,bdot)
    integer, intent(in) :: icol,nfreq,nmtx,irow_G0,iscreen,icutv
    real(DP), intent(in) :: vcoul,oneoverq,q0len,fact
    complex(DPC), intent(inout) :: epstemp(nmtx,nfreq)
    logical, intent(in) :: q0flag
    real(DP), intent(in) :: bdot(3,3)
    real(DP) :: zc
    integer :: i
!----------------------
! No Truncation
   
! JRD This routine has horrible locality. Need to loop over iw on outside
    if (icutv==TRUNC_NONE) then
      if (icol.ne.irow_G0) then
        if (iscreen==SCREEN_SEMICOND) then
          epstemp(irow_G0,:) = epstemp(irow_G0,:)*oneoverq*q0len/(8D0*PI_D)
        endif
        if (q0flag .and. iscreen==SCREEN_SEMICOND) then
          epstemp(irow_G0,:) = 0d0
        endif
      else
        do i=1,nmtx
          if (i .ne. irow_G0) then
            if (iscreen==SCREEN_SEMICOND) then
              epstemp(i,:) = epstemp(i,:)*oneoverq*fact/(vcoul*q0len)
            endif
            if (iscreen==SCREEN_GRAPHENE) then
              epstemp(i,:) = epstemp(i,:)*fact*8D0*PI_D/(vcoul*q0len**2)
            endif
            if (q0flag .and. iscreen==SCREEN_SEMICOND) then
              epstemp(i,:) = 0d0
            endif
          endif
        enddo
      endif
    endif
!----------------------
! Cell Wire Truncation
! May not be implemented correctly for graphene screening... I`m not even sure there is a 1D system with linear DOS...
    if (icutv==TRUNC_WIRE) then
      if (icol.ne.irow_G0) then
        if (q0flag .and. iscreen==SCREEN_SEMICOND) then
          epstemp(irow_G0,:) = 0d0
        endif
      else
        do i=1,nmtx
          if (i .ne. irow_G0) then
            if (q0flag .and. iscreen==SCREEN_SEMICOND) then
              epstemp(i,:) = 0d0
            endif
          endif
        enddo
      endif
    endif
!----------------------
! Cell Slab Truncation
    zc=2D0*PI_D/(sqrt(bdot(3,3))*2D0)
    if (icutv==TRUNC_SLAB) then
      if (icol.ne.irow_G0) then
        if (q0flag .and. iscreen==SCREEN_SEMICOND) then
          epstemp(irow_G0,:) = 0d0
        endif
      else
        do i=1,nmtx
          if (i .ne. irow_G0) then
            if (iscreen/=SCREEN_METAL) then
              epstemp(i,:) = epstemp(i,:) * 8D0 * PI_D * fact * zc &
                / (vcoul * q0len)
              if (q0flag .and. iscreen==SCREEN_SEMICOND) then
                epstemp(i,:) = 0d0
              endif
            endif
          endif
        enddo
      endif
    endif
   
    return
  end subroutine fixwings_dyn
end module fixwings_m
