!===========================================================================
!
! Routines()
!
! (1) shiftenergy() Originally by ? Last Edited: Apr/2016 (FHJ)
!
! Computes and symmetrizes the quasiparticle spectrum, and
! perform linear interpolation/extrapolation to solve Dyson`s equation.
! Note that, since we may use more than two frequency points to evaluate
! Sigma(omega), our "eqp1" in general goes beyond the linear correction
! scheme in Eq. (37) of Hybertsen & Louie PRB.
!
!===========================================================================
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
module shiftenergy_m
  use global_m
  implicit none
  private
  public :: &
    shiftenergy
contains
subroutine shiftenergy(sig, wfnk, alda, asx, ach, achcor, ach_n1, achcor_n1, &
  ax, efsto, asig, enew, zrenorm, nfreqgpp, ncore_excl, akih) ! ZL: akih attached as optional
  type (siginfo), intent(in) :: sig
  type (wfnkstates), intent(in) :: wfnk
  integer, intent(in) :: nfreqgpp
  real(DP), intent(inout) :: &
    alda(sig%ndiag+sig%noffdiag,sig%nspin), & !< vxc
    asx(nfreqgpp,sig%ndiag+sig%noffdiag,sig%nspin), & !< sx
    ach(nfreqgpp,sig%ndiag+sig%noffdiag,sig%nspin) !< ch
  complex(DPC), intent(inout) :: achcor(sig%ndiag+sig%noffdiag,sig%nspin) !< static remainder
  real(DP), intent(inout) :: ach_n1(sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin)
  real(DP), intent(inout) :: achcor_n1(sig%ntband,sig%ndiag+sig%noffdiag,sig%nspin)
  real(DP), intent(inout) :: ax(sig%ndiag+sig%noffdiag,sig%nspin) !< x
  real(DP), intent(out) :: efsto(sig%ndiag,sig%nspin) !< eqp0
  real(DP), intent(inout) :: asig(sig%ndiag+sig%noffdiag,sig%nspin) !< sig
  real(DP), intent(out) :: enew(sig%ndiag,sig%nspin) !< eqp1
  real(DP), intent(out) :: zrenorm(sig%ndiag,sig%nspin) !< Znk
  integer, intent(in) :: ncore_excl !< number of core states excluded
  real(DP), intent(inout), optional :: akih(sig%ndiag+sig%noffdiag,sig%nspin) !< kih
  integer :: ii,jj,istart,istop,nsubs,ispin,iwlda0,iwlda1,iwlda2,iw
  integer :: ndeg(sig%ndiag)
  real(DP) :: fact,dek,dele,diffmin,delta_omega,diff,e_lk,freq0
  real(DP) :: aldai,axi, akihi ! ZL: add akihi
  real(DP), dimension(nfreqgpp) :: asigi, asxi, achi
  real(DP), dimension(sig%ntband) :: ach_n1i, achcor_n1i
  complex(DPC) :: achcori
 
  ! ZL: check optional akih
  if(sig%use_kihdat) then
    if (.not.present(akih)) call die('In shiftenergy, to use KIH, pass akih into this subroutine.')
  endif
  do ispin=1,sig%nspin
    nsubs = 1
    ndeg(nsubs) = 1
    ! This loop is used to find degeneracy
    do ii=2,sig%ndiag
! DVF : ncore_excl has to be substracted here because wfnk%elda is defined in the
! read_wavefunction subroutine in input.f90 to be referenced to the case with
! no core states. Same with wfnk%el/elda below.
      dek = wfnk%elda(sig%diag(ii)-ncore_excl,ispin) - wfnk%elda(sig%diag(ii-1)-ncore_excl,ispin)
      if (abs(dek)<sig%tol .and. sig%symmetrize) then
        ! Band ii is degenerate to ii-1
        ndeg(nsubs) = ndeg(nsubs) + 1
      else
        ! Start a new degenerate subspace
        nsubs = nsubs + 1
        ndeg(nsubs) = 1
      endif
    enddo
    istop = 0
    ! loop degeneracy
    do ii=1,nsubs
      istart = istop + 1
      istop = istart + ndeg(ii) - 1
      aldai = 0.0d0
      akihi = 0.0d0 ! ZL: for KIH
      axi = 0.0d0
      asxi = 0.0d0
      achi = 0.0d0
      achcori = (0.0d0, 0.0d0)
      ach_n1i(:) = 0.0d0
      achcor_n1i(:) = 0.0d0
      ! start real loop over bands
      do jj=istart,istop
        ! average VXC in the same degen group: combine first
        aldai = aldai + alda(jj,ispin)
        ! "i" is intermediate steps in Ry
        ! ZL: do the same for KIH as VXC
        akihi = akihi + akih(jj,ispin)
        ! Note: while using Vxc or KIH, the other array would just be zeros
        axi = axi + ax(jj,ispin)
        asxi(:) = asxi(:) + asx(:,jj,ispin)
        achi(:) = achi(:) + ach(:,jj,ispin)
        achcori = achcori + achcor(jj,ispin)
        ach_n1i(:) = ach_n1i(:) + ach_n1(:,jj,ispin)
        achcor_n1i(:) = achcor_n1i(:) + achcor_n1(:,jj,ispin)
      enddo
      fact = ryd / dble(ndeg(ii))
      do jj=istart,istop
        ! then AVERAGE
        ! common quantities in eV
        alda(jj,ispin) = aldai * fact
        ! ZL: do the same for KIH
        akih(jj,ispin) = akihi * fact
        ax(jj,ispin) = axi * fact
        asx(:,jj,ispin) = asxi(:) * fact
        ach(:,jj,ispin) = achi(:) * fact
        achcor(jj,ispin) = achcori * fact
        ach_n1(:,jj,ispin) = ach_n1i(:) * fact
        achcor_n1(:,jj,ispin) = achcor_n1i(:) * fact
        ! intermediate the quantity, sigma in Ry, unnecessary i
        asigi(:) = ax(jj,ispin) + asx(:,jj,ispin) + ach(:,jj,ispin)
        ! FHJ: Determine the expansion points for the linear interpolation: iwlda0
        ! is where we evaluate eqp0, and iwlda1, iwlda2 are the expansion points
        ! we used in the finite difference linearization (ie, used to compute the
        ! derivatives)
        if (sig%fdf==-3) then
          ! FHJ: Find iw closest to e_lk
          iwlda1 = 1
          iwlda2 = 2
          diffmin = INF
          e_lk = wfnk%ek(sig%diag(jj)-ncore_excl,ispin)
          ! FHJ: Figure out starting frequency for freq. grid
          if (sig%freq_grid_shift<2) then
            freq0 = sig%freqevalmin
          else
            freq0 = e_lk - sig%freqevalstep*(sig%nfreqeval-1)/2
          endif
          do iw=1,sig%nfreqeval
            diff = abs(freq0 + (iw-1)*sig%freqevalstep - e_lk)
            if (diff .lt. diffmin) then
              diffmin=diff
              iwlda2=iwlda1
              iwlda1=iw
            endif
          enddo
          iwlda0 = iwlda1
          delta_omega = (iwlda2 - iwlda1)*sig%freqevalstep
        else
          !iwlda0 is always 2 here
          iwlda0 = 2
          iwlda1 = 2
          iwlda2 = 2
          select case (sig%fdf)
            case (-1)
              iwlda1 = 1
              iwlda2 = 2
            case (0)
              iwlda1 = 1
              iwlda2 = 3
            case (1,2)
              iwlda1 = 2
              iwlda2 = 3
          end select
          delta_omega = (iwlda2 - iwlda1)*sig%dw
        endif
        asig(jj,ispin) = asigi(iwlda0)
!FHJ: Perform linear interpolation/extrapolation for all finite difference
! schemes at once. enew is the interpolated off-shell answer = "eqp1",
! and efsto is the on-shell answer = "eqp0".
        ! ZL: Vxc vs. KIH
        if(.not. sig%use_kihdat) then
          efsto(jj,ispin) = wfnk%elda(sig%diag(jj)-ncore_excl,ispin) - &
            alda(jj,ispin) + asig(jj,ispin)
        else
          efsto(jj,ispin) = akih(jj,ispin) + asig(jj,ispin)
        endif
        dele = efsto(jj,ispin) - wfnk%ek(sig%diag(jj)-ncore_excl,ispin)
        if (iwlda1/=iwlda2) then
          enew(jj,ispin) = efsto(jj,ispin) + &
            (asigi(iwlda2)-asigi(iwlda1)) / (delta_omega - asigi(iwlda2) + asigi(iwlda1)) * dele
          zrenorm(jj, ispin) = 1d0 / (1d0 - (asigi(iwlda2) - asigi(iwlda1))/delta_omega)
        else
          enew(jj,ispin) = efsto(jj,ispin)
          zrenorm(jj, ispin) = 1d0
        endif
      enddo ! jj
    enddo ! ii
  enddo ! ispin
 
  return
end subroutine shiftenergy
end module shiftenergy_m
