!===========================================================================
!
! Routines:
!
! (1) write_result_hp() Originally By ? Last Modified 7/3/2008 (JRD)
!
! Writes the quasiparticle spectrum to the output
!
!===========================================================================

module write_result_hp_m
  use global_m
  implicit none
  private
  public :: write_result_hp
contains
subroutine write_result_hp(kp,wfnk,sig,ax,asx,ach,achcor,asig,alda,efsto,&
                           enew,zrenorm,ikn,ncore_excl, akih)
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
  integer, intent(in) :: ikn
  integer, intent(in) :: ncore_excl !< number of core states excluded
  real(DP), intent(in), optional :: akih(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin), ZL
  integer :: i, j, ispin, iw, nstart, nend
 
  ! ZL: check optional akih
  if(sig%use_kihdat) then
    if (.not.present(akih)) call die('In write_result_hp, to use KIH, pass akih into this subroutine.')
  endif
  do ispin=1,sig%nspin
! Sigma Diagonal
    write(8,979) (kp%rk(j,ikn),j=1,3),ikn,sig%spin_index(ispin)
    if (sig%freq_dep .eq. -1 .or. &
       (sig%freq_dep .eq. 0 .and. sig%exact_ch .eq. 1) .or. &
       (sig%freq_dep .eq. 1 .and. sig%exact_ch .eq. 0)) then
      ! ZL: Vxc vs. KIH
      if(.not.sig%use_kihdat) then
        write(8,978)
        write(8,977) (sig%diag(i),wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
          wfnk%ek(sig%diag(i)-ncore_excl,ispin),dble(ax(i,ispin)),dble(asx(2,i,ispin)), &
          dble(ach(2,i,ispin)+achcor(i,ispin)), &
          dble(asig(i,ispin)+achcor(i,ispin)),dble(alda(i,ispin)), &
          efsto(i,ispin)+dble(achcor(i,ispin)), &
          enew(i,ispin)+dble(achcor(i,ispin))*zrenorm(i,ispin), &
          zrenorm(i,ispin),i=1,sig%ndiag)
      else
        write(8,9781)
        write(8,977) (sig%diag(i),wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
          wfnk%ek(sig%diag(i)-ncore_excl,ispin),dble(ax(i,ispin)),dble(asx(2,i,ispin)), &
          dble(ach(2,i,ispin)+achcor(i,ispin)), &
          dble(asig(i,ispin)+achcor(i,ispin)),dble(akih(i,ispin)), &
          efsto(i,ispin)+dble(achcor(i,ispin)), &
          enew(i,ispin)+dble(achcor(i,ispin))*zrenorm(i,ispin), &
          zrenorm(i,ispin),i=1,sig%ndiag)
      endif
    else
      if(.not.sig%use_kihdat) then
        write(8,976)
        write(8,975) (sig%diag(i),wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
          wfnk%ek(sig%diag(i)-ncore_excl,ispin),dble(ax(i,ispin)),dble(asx(2,i,ispin)), &
          dble(ach(2,i,ispin)+achcor(i,ispin)), &
          dble(asig(i,ispin)+achcor(i,ispin)),dble(alda(i,ispin)), &
          efsto(i,ispin)+dble(achcor(i,ispin)), &
          enew(i,ispin)+dble(achcor(i,ispin))*zrenorm(i,ispin), &
          dble(ach(2,i,ispin)),dble(asig(i,ispin)), &
          efsto(i,ispin),enew(i,ispin),zrenorm(i,ispin),i=1,sig%ndiag)
      else
        write(8,9761)
        write(8,975) (sig%diag(i),wfnk%elda(sig%diag(i)-ncore_excl,ispin), &
          wfnk%ek(sig%diag(i)-ncore_excl,ispin),dble(ax(i,ispin)),dble(asx(2,i,ispin)), &
          dble(ach(2,i,ispin)+achcor(i,ispin)), &
          dble(asig(i,ispin)+achcor(i,ispin)),dble(akih(i,ispin)), &
          efsto(i,ispin)+dble(achcor(i,ispin)), &
          enew(i,ispin)+dble(achcor(i,ispin))*zrenorm(i,ispin), &
          dble(ach(2,i,ispin)),dble(asig(i,ispin)), &
          efsto(i,ispin),enew(i,ispin),zrenorm(i,ispin),i=1,sig%ndiag)
      endif
    endif
979 format(7x,"k =",3f10.6,1x,"ik =",i4,1x,"spin =",i2)
978 format(/,3x,"n",7x,"  Emf",7x,"   Eo",7x,"    X",7x," SX-X",7x,"   CH", &
                    7x,"  Sig",7x,"  Vxc",7x," Eqp0",7x," Eqp1",7x, "Znk")
9781 format(/,3x,"n",7x,"  Emf",7x,"   Eo",7x,"    X",7x," SX-X",7x,"   CH", &
                     7x,"  Sig",7x,"  KIH",7x," Eqp0",7x," Eqp1",7x, "Znk")
977 format(i4,10f12.6)
976 format(/,3x,"n",7x,"  Emf",7x,"   Eo",7x,"    X",7x," SX-X",7x,"   CH", &
                    7x,"  Sig",7x,"  Vxc",7x," Eqp0",7x," Eqp1",7x,"  CH`", &
                    7x," Sig`",7x,"Eqp0`",7x,"Eqp1`",7x,"  Znk")
9761 format(/,3x,"n",7x,"  Emf",7x,"   Eo",7x,"    X",7x," SX-X",7x,"   CH", &
                     7x,"  Sig",7x,"  KIH",7x," Eqp0",7x," Eqp1",7x,"  CH`", &
                     7x," Sig`",7x,"Eqp0`",7x,"Eqp1`",7x,"  Znk")
975 format(i4,14f12.6)
! Sigma Off-Diagonal
    if (sig%noffdiag.gt.0) then
      if (sig%fdf.eq.-1) then
        nstart = 1
        nend = 2
      elseif (sig%fdf.eq.0) then
        nstart = 1
        nend = 3
      elseif (sig%fdf.eq.1) then
        nstart = 2
        nend = 3
      else
        nstart = 2
        nend = 2
      endif
      do iw = nstart, nend
        if (iw.eq.1) write(8,951)
        if (iw.eq.2) write(8,952)
        if (iw.eq.3) write(8,953)
        if (sig%freq_dep .eq. -1 .or. &
           (sig%freq_dep .eq. 0 .and. sig%exact_ch .eq. 1) .or. &
           (sig%freq_dep .eq. 1 .and. sig%exact_ch .eq. 0)) then
          if(.not.sig%use_kihdat) then
            write(8,969) iw, iw, iw
          else
            write(8,9691) iw, iw, iw
          endif
          do i=sig%ndiag+1,sig%ndiag+sig%noffdiag
            if(.not.sig%use_kihdat) then
              write(8,968) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
                sig%off3(i-sig%ndiag),dble(ax(i,ispin)), &
                dble(asx(iw,i,ispin)),dble(ach(iw,i,ispin)+achcor(i,ispin)), &
                dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
                dble(alda(i,ispin))
            else
              write(8,968) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
                sig%off3(i-sig%ndiag),dble(ax(i,ispin)), &
                dble(asx(iw,i,ispin)),dble(ach(iw,i,ispin)+achcor(i,ispin)), &
                dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
                dble(akih(i,ispin))
            endif
          enddo
        else
          if(.not.sig%use_kihdat) then
            write(8,966) iw, iw, iw, iw, iw
          else
            write(8,9661) iw, iw, iw, iw, iw
          endif
          do i=sig%ndiag+1,sig%ndiag+sig%noffdiag
            if(.not.sig%use_kihdat) then
              write(8,965) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
                sig%off3(i-sig%ndiag),dble(ax(i,ispin)), &
                dble(asx(iw,i,ispin)),dble(ach(iw,i,ispin)+achcor(i,ispin)), &
                dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
                dble(alda(i,ispin)),dble(ach(iw,i,ispin)), &
                dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin))
            else
              write(8,965) sig%off1(i-sig%ndiag),sig%off2(i-sig%ndiag), &
                sig%off3(i-sig%ndiag),dble(ax(i,ispin)), &
                dble(asx(iw,i,ispin)),dble(ach(iw,i,ispin)+achcor(i,ispin)), &
                dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin)+achcor(i,ispin)), &
                dble(akih(i,ispin)),dble(ach(iw,i,ispin)), &
                dble(ax(i,ispin)+asx(iw,i,ispin)+ach(iw,i,ispin))
            endif
          enddo
        endif
      enddo
    endif
951 format(/,3x,"E1 = E0 - dE")
952 format(/,3x,"E2 = E0")
953 format(/,3x,"E3 = E0 + dE")
969 format(3x,"n",3x,"m",3x,"l",21x,"X",4x,"SX(E",i1,")-X",6x, &
      "CH(E",i1,")",5x,"Sig(E",i1,")",9x,"Vxc")
9691 format(3x,"n",3x,"m",3x,"l",21x,"X",4x,"SX(E",i1,")-X",6x, &
      "CH(E",i1,")",5x,"Sig(E",i1,")",9x,"KIH")
968 format(3i4,3x,"real",3x,5f12.6)
967 format(3i4,3x,"imag",3x,5f12.6)
966 format(3x,"n",3x,"m",3x,"l",21x,"X",4x,"SX(E",i1,")-X",6x, &
      "CH(E",i1,")",5x,"Sig(E",i1,")",9x,"Vxc",5x,"CH`(E",i1,")",4x, &
      "Sig`(E",i1,")")
9661 format(3x,"n",3x,"m",3x,"l",21x,"X",4x,"SX(E",i1,")-X",6x, &
      "CH(E",i1,")",5x,"Sig(E",i1,")",9x,"KIH",5x,"CH`(E",i1,")",4x, &
      "Sig`(E",i1,")")
965 format(3i4,3x,"real",3x,7f12.6)
964 format(3i4,3x,"imag",3x,7f12.6)
    write(8,*)
  enddo ! ispin
 
  return
end subroutine write_result_hp
end module write_result_hp_m
