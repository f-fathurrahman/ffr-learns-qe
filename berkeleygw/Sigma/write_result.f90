!===========================================================================
!
! Routines:
!
! (1) write_result() Originally By ? Last Modified 7/3/2008 (JRD)
!
! Writes the quasiparticle spectrum to the output.
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
module write_result_m
  use global_m
  implicit none
  private
  public :: write_result
contains
subroutine write_result(kp, wfnk, sig, ax, asx, ach, achcor, asig, alda, efsto, enew, &
  zrenorm, ikn, ncore_excl, akih)
  type (kpoints), intent(in) :: kp
  type (wfnkstates), intent(in) :: wfnk
  type (siginfo), intent(in) :: sig
  real(DP), intent(in) :: ax(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin)
  real(DP), intent(in) :: asx(:,:,:) !< (nfreqgpp,sig%ndiag+sig%noffdiag,sig%nspin)
  real(DP), intent(in) :: ach(:,:,:) !< (nfreqgpp,sig%ndiag+sig%noffdiag,sig%nspin)
  complex(DPC), intent(in) :: achcor(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin)
  real(DP), intent(in) :: asig(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin)
  real(DP), intent(in) :: alda(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin)
  real(DP), intent(in) :: efsto(:,:) !< (sig%ndiag,sig%nspin)
  real(DP), intent(in) :: enew(:,:) !< (sig%ndiag,sig%nspin)
  real(DP), intent(in) :: zrenorm(:,:) !< (sig%ndiag,sig%nspin)
  integer, allocatable :: iwlda(:)
  integer, intent(in) :: ikn
  integer, intent(in) :: ncore_excl !< number of core states excluded
  real(DP), intent(in), optional :: akih(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin) ZL
  integer :: i, j, ispin, iw
  real(DP) :: eval, e_lk, freq0, diffmin, diff
  !> FHJ: CH without static remainder
  real(DP), dimension(sig%ntband) :: chp_v, chp_c, chp_d
  !> FHJ: Static remainder
  real(DP), dimension(sig%ntband) :: sr_v, sr_c
  !> FHJ: CH with static remainder
  real(DP), dimension(sig%ntband) :: ch_v, ch_c, ch_d
  real(DP) :: extrap_v, extrap_c, extrap_d, extrapp_v, extrapp_c, extrapp_d
! Initialization
 
  ! ZL: check optional akih
  if(sig%use_kihdat) then
    if (.not.present(akih)) call die('In write_result, to use KIH, pass akih into this subroutine.')
  endif
  allocate(iwlda (sig%ndiag))
  do ispin=1,sig%nspin
! Sigma Diagonal
! DVF : ncore_excl has to be substracted here because wfnk%ek is defined in the
! read_wavefunction subroutine in input.f90 to be referenced to the case with
! no core states.
! FHJ: writes the spectrum.dat file, if GPP calculation with many freqs.
    if (sig%fdf .eq. -3) then
      do i=1, sig%ndiag
          diffmin = INF
          e_lk = wfnk%ek(sig%diag(i)-ncore_excl,ispin)
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
              iwlda(i)=iw
            endif
          enddo
          write(8000,2001)
          do iw=1,sig%nfreqeval
            eval = freq0 + (iw-1)*sig%freqevalstep
            write(8000,2000) kp%rk(:,ikn), ispin, sig%diag(i), iw, eval, &
                          dble(asx(iw,i,ispin)+ax(i,ispin)), &
                          dble(ach(iw,i,ispin))+dble(achcor(i,ispin)), &
                          dble(asx(iw,i,ispin)+ax(i,ispin)+ach(iw,i,ispin))+dble(achcor(i,ispin))
          enddo
          write(8000,*) ''
      enddo
    else
      iwlda(:) = 2
    endif
2000 format(3F12.5,2x,3i4,2x,4F12.5)
2001 format("#",6x,"kx",10x,"ky",10x,"kz",7x,"spn",1x,"bnd",2x,"iw",8x,"Ew",3x, &
    "      Re(SX)", "      Re(CH)", "   Re(SX+CH)")
    call write_qp(6)
  enddo ! ispin
  if(allocated(iwlda))then;deallocate(iwlda);endif
 
  return
contains
  subroutine write_qp(iunit)
    integer, intent(in) :: iunit
    write(iunit,979) (kp%rk(j,ikn),j=1,3),ikn,sig%spin_index(ispin)
979 format(/,7x,"k =",3f10.6,1x,"ik =",i4,1x,"spin =",i2)
    write(iunit,'()')
    if(.not.sig%use_kihdat) then
      write(iunit,900) 'n', 'Emf', 'Eo', 'Vxc', 'X', 'Cor', 'Eqp0', 'Eqp1', 'Znk'
    else
      write(iunit,900) 'n', 'Emf', 'Eo', 'KIH', 'X', 'Cor', 'Eqp0', 'Eqp1', 'Znk'
    endif
900 format(a6,8(a9))
! JRD: Not ideal for fdf -3, but user should see spectrum.dat
    do i = 1, sig%ndiag
      ! ZL: Vxc vs. KIH
      if(.not.sig%use_kihdat) then
        write(iunit,977) sig%diag(i), wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
          wfnk%ek(sig%diag(i)-ncore_excl,ispin), dble(alda(i,ispin)), dble(ax(i,ispin)), &
          dble(asx(iwlda(i),i,ispin)+ach(iwlda(i),i,ispin)+achcor(i,ispin)), &
          efsto(i,ispin)+dble(achcor(i,ispin)), &
          enew(i,ispin)+dble(achcor(i,ispin))*zrenorm(i,ispin), zrenorm(i,ispin)
      else
        write(iunit,977) sig%diag(i), wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
          wfnk%ek(sig%diag(i)-ncore_excl,ispin), dble(akih(i,ispin)), dble(ax(i,ispin)), &
          dble(asx(iwlda(i),i,ispin)+ach(iwlda(i),i,ispin)+achcor(i,ispin)), &
          efsto(i,ispin)+dble(achcor(i,ispin)), &
          enew(i,ispin)+dble(achcor(i,ispin))*zrenorm(i,ispin), zrenorm(i,ispin)
      endif
    enddo
977 format(i6,10f9.3)
! Sigma Off-Diagonal
! JRD: Not ideal for fdf -3, but user should see spectrum.dat
    if (sig%noffdiag.gt.0) then
      if(.not.sig%use_kihdat) then
        write(iunit,969)
      else
        write(iunit,9691)
      endif
!969 format(/,5x,"n",5x,"m",5x,"l",6x,6x,"Vxc",8x,"X",5x,"SX-X",7x,"CH",6x, &
! "Sig",6x,"Vxc")
969 format(/,5x,"n",5x,"m",5x,"l",6x,6x,"Vxc",8x,"X",5x,"SX-X",7x,"CH",6x, &
           "Sig") ! ZL fixed the bug: redundant Vxc at the end
9691 format(/,5x,"n",5x,"m",5x,"l",6x,6x,"KIH",8x,"X",5x,"SX-X",7x,"CH",6x, &
           "Sig")
      do i=sig%ndiag+1,sig%ndiag+sig%noffdiag
        if(.not.sig%use_kihdat) then
          write(iunit,968) sig%off1(i-sig%ndiag), sig%off2(i-sig%ndiag), &
            sig%off3(i-sig%ndiag), dble(alda(i,ispin)), dble(ax(i,ispin)), &
            dble(asx(2,i,ispin)), dble(ach(2,i,ispin)+achcor(i,ispin)), &
            dble(ax(i,ispin)+asx(2,i,ispin)+ach(2,i,ispin)+achcor(i,ispin))
        else
          write(iunit,968) sig%off1(i-sig%ndiag), sig%off2(i-sig%ndiag), &
            sig%off3(i-sig%ndiag), dble(akih(i,ispin)), dble(ax(i,ispin)), &
            dble(asx(2,i,ispin)), dble(ach(2,i,ispin)+achcor(i,ispin)), &
            dble(ax(i,ispin)+asx(2,i,ispin)+ach(2,i,ispin)+achcor(i,ispin))
        endif
968 format(3i6,1x,"real",1x,5f9.3)
      enddo
    endif
  end subroutine write_qp
end subroutine write_result
end module write_result_m
