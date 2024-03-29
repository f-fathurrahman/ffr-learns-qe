!===============================================================================
!
! Routines:
!
! 1. checkbz() Originally By gsm Last Modified 7/29/2010 (gsm)
!
! Check that the Brillouin Zone generated by subroutine fullbz is
! identical to the original full Brillouin Zone. Subroutine fullbz
! constructs the Brillouin Zone by unfolding the irreducible wedge
! with all the symmetries of the space group of the crystal. If the
! irreducible wedge has missing k-points (for example from using too
! many symmetries in kgrid.x), the full Brillouin Zone will also have
! missing k-points.
!
! For the unshifted grid, fullbz generates the original full grid
! from the irreducible wedge.
!
! For the grid shifted by half a grid step, fullbz doubles the grid size.
! For fcc-Si, (4 4 4 0.5 0.5 0.5) becomes (8 8 8 0.0 0.0 0.0) where half
! the points uniformly distributed across the grid are missing.
!
! For the randomly-shifted grid, if symmetries are allowed, fullbz
! generates a non-uniform grid with the points clustered together.
! You should never allow symmetries for the randomly-shifted grid.
!
! For the grid shifted by half a grid step, checkbz would print "extra points"
! warning message. There is nothing to worry about, so we set
! allow_half_shift to .true. to suppress this warning message in this case.
! However, if you see "missing points" warning message, it may indicate a
! problem with your k-point sampling.
!
!===============================================================================

module checkbz_m
  use global_m
  use misc_m
  implicit none
  private
  public :: checkbz
contains
subroutine checkbz(nfk,fk,kgrid,kshift,bdot, &
  filename,kqchar,wignerseitz,freplacebz,fwritebz)
  integer, intent(inout) :: nfk
  real(DP), intent(inout) :: fk(:,:) !< (3, nfk)
  integer, intent(in) :: kgrid(3)
  real(DP), intent(in) :: kshift(3)
  real(DP), intent(in) :: bdot(3,3)
  character(len=*), intent(in) :: filename
  character, intent(in) :: kqchar
  logical, intent(in) :: wignerseitz,freplacebz,fwritebz
  logical, parameter :: allow_half_shift = .true.
  integer, parameter :: iunit_checkbz=50
  logical :: f1,f2,f3,flag_half_shift
  logical :: any_warning
  integer :: ii,jj,i1,i2,i3,nk,gpt(3)
  real(DP) :: l1,l2,k1(3),k2(3),kpt(3)
  real(DP), allocatable :: kk(:,:)
  real(DP), allocatable :: kref(:,:)
  character(len=128) :: tmpstr
 
  any_warning = .false.
  if(size(fk, 1) /= 3) then
    write(0,*) 'size(fk, 1) = ', size(fk, 1)
    call die("checkbz internal error: fk must have first dimension = 3")
  endif
  if(size(fk, 2) /= nfk) then
    write(0,*) 'nfk = ', nfk, 'size(fk, 2) = ', size(fk, 2)
    call die("checkbz internal error: fk must have second dimension = nfk")
  endif
! Identify the grid type
! Print a warning message if the grid type is unknown
! FHJ: Is this really necessary?!
  kpt(:)=kshift(:)
  call k_range(kpt, gpt, TOL_Small)
  f1=.true.
  f2=.true.
  f3=.true.
  do ii=1,3
! FHJ: ignore dimensions that have only one kpt
    if (kpt(ii)<2) cycle
  ! the unshifted grid
    f1=f1.and.(abs(kpt(ii)).lt.TOL_Small)
  ! the grid shifted by half a grid step
    f2=f2.and.(abs(kpt(ii)-0.5d0).lt.TOL_Small)
  ! the randomly-shifted grid
    f3=f3.and.(abs(kpt(ii)).gt.TOL_Small.and. abs(kpt(1)-0.5d0).gt.TOL_Small)
  enddo
  if (.not.f1.and..not.f2.and..not.f3) then
    call auto_open_file()
    if (peinf%inode.eq.0) then
      write(0,901) kqchar, trim(filename)
      write(iunit_checkbz,901) kqchar, trim(filename)
    endif
901 format(1x,"WARNING: checkbz: unknown",1x,a,"-grid type in",1x,a,/)
  endif
  flag_half_shift = f2
! Find the number of k-points in the full Brillouin Zone
  nk=product(kgrid(1:3))
  if (nk.le.0) then
    call auto_open_file()
    if (peinf%inode==0) then
      write(0,902) kqchar, trim(filename)
      write(iunit_checkbz,902) kqchar, trim(filename)
    endif
902 format(1x,"WARNING: checkbz: zero",1x,a,"-grid in",1x,a,/)
    call auto_close_file()
   
    return
  endif
! Allocate array for k-points in the full Brillouin Zone
  allocate(kk (3,nk))
! Construct k-points in the full Brillouin Zone
  ii=0
  do i1=0,kgrid(1)-1
    do i2=0,kgrid(2)-1
      do i3=0,kgrid(3)-1
        ii=ii+1
        kk(1,ii)=(dble(i1)+kshift(1))/dble(kgrid(1))
        kk(2,ii)=(dble(i2)+kshift(2))/dble(kgrid(2))
        kk(3,ii)=(dble(i3)+kshift(3))/dble(kgrid(3))
        call k_range(kk(:,ii), gpt, TOL_Small)
      enddo
    enddo
  enddo
! Construct a Wigner-Seitz box
  if (wignerseitz) then
    do ii=1,nk
      l2=INF
      do i1=-ncell+1,ncell
        k1(1)=kk(1,ii)-dble(i1)
        do i2=-ncell+1,ncell
          k1(2)=kk(2,ii)-dble(i2)
          do i3=-ncell+1,ncell
            k1(3)=kk(3,ii)-dble(i3)
            l1=DOT_PRODUCT(k1,MATMUL(bdot,k1))
            if (l1.lt.l2) then
              l2=l1
              k2(:)=k1(:)
            endif
          enddo
        enddo
      enddo
      kk(:,ii)=k2(:)
    enddo
  endif
! Write unfolded BZ and full BZ to files
  if (fwritebz) then
    if (peinf%inode.eq.0) then
      write(tmpstr,801) kqchar, trim(filename)
801 format(a,"_",a,"_unfolded.dat")
      call open_file(14, tmpstr, status='replace', form='formatted')
      write(14,803) nfk
      do ii=1,nfk
        write(14,804) ii,fk(:,ii)
      enddo
      call close_file(14)
      write(tmpstr,802) kqchar, trim(filename)
802 format(a,"_",a,"_full.dat")
      call open_file(14, tmpstr, status='replace', form='formatted')
      write(14,803) nk
      do ii=1,nk
        write(14,804) ii,kk(:,ii)
      enddo
      call close_file(14)
    endif
  endif
! Replace unfolded BZ with full BZ
  if (freplacebz) then
    if (nk.le.nfk) then
      nfk=nk
      fk(1:3,1:nk)=kk(1:3,1:nk)
    else
      call die('checkbz: failed replacebz')
    endif
  endif
! Before comparing k-points translate from Wigner-Seitz box
! to [0,1) interval
  allocate(kref (3,nfk))
  kref(1:3,1:nfk)=fk(1:3,1:nfk)
  if (wignerseitz) then
    do ii=1,nfk
      call k_range(kref(:,ii), gpt, TOL_Small)
    enddo
    do ii=1,nk
      call k_range(kk(:,ii), gpt, TOL_Small)
    enddo
  endif
! Check that kref(1:3,1:nfk) is a subset of kk(1:3,1:nk)
! Print a warning message otherwise
  if(.not.(flag_half_shift.and.allow_half_shift)) then
    f1=.true.
    do ii=1,nfk
      f3=.false.
      do jj=1,nk
        if (all(abs(kref(1:3,ii)-kk(1:3,jj)).lt.TOL_Small)) f3=.true.
      enddo
      if (.not.f3) then
        f1=.false.
        call auto_open_file()
        if (peinf%inode==0) then
          write(iunit_checkbz,'(a,3f12.6)') 'Extra point: ', kref(1:3,ii)
        endif
      endif
    enddo
    if (.not.f1) then
      if (peinf%inode==0) then
        write(0,'(/5a)') 'WARNING: checkbz: unfolded BZ from ', &
          trim(filename), ' has extra ', kqchar, '-points'
        write(0,'(a/)') 'See file checkbz.log for more information.'
      endif
    endif
  endif
! Check that kk(1:3,1:nk) is a subset of kref(1:3,1:nfk)
! Print a warning message otherwise
  f2=.true.
  do ii=1,nk
    f3=.false.
    do jj=1,nfk
      if (all(abs(kk(1:3,ii)-kref(1:3,jj)).lt.TOL_Small)) f3=.true.
    enddo
    if (.not.f3) then
      f2=.false.
      call auto_open_file()
      if (peinf%inode==0) then
        write(iunit_checkbz,'(a,3f12.6)') 'Missing point: ', kk(1:3,ii)
      endif
    endif
  enddo
  if (.not.f2) then
    if (peinf%inode==0) then
      if (peinf%inode==0) then
        write(0,'(/5a)') 'WARNING: checkbz: unfolded BZ from ', &
          trim(filename), ' has missing ', kqchar, '-points'
        if (trim(filename)=='epsilon.inp' .and. kqchar=='q') then
          write(0,'(a)') &
            '(disregard this warning if your epsilon calculation is split by q-points)'
        endif
        write(0,'(a/)') 'See file checkbz.log for more information.'
      endif
    endif
  endif
! Deallocate and finish
  if(allocated(kk))then;deallocate(kk);endif
  if(allocated(kref))then;deallocate(kref);endif
  call auto_close_file()
 
  return
803 format(i5)
804 format(i5,3f13.9)
contains
  subroutine auto_open_file
    character :: adate*11, atime*14
   
    if (.not.any_warning) then
      any_warning = .true.
      if (peinf%inode==0) then
        call date_time(adate, atime)
        call open_file(iunit_checkbz, 'checkbz.log', status='UNKNOWN', form='formatted')
        write(iunit_checkbz,'(/1x,5a/)') 'Output from checkbz on ', &
          trim(adate), ' at ', trim(atime), ':'
      endif
    endif
   
  end subroutine auto_open_file
  subroutine auto_close_file
   
    if (any_warning.and.peinf%inode==0) then
      write(iunit_checkbz,'(/a/)') &
        '--------------------------------------------------------------------------------'
      call close_file(iunit_checkbz)
    endif
   
  end subroutine auto_close_file
end subroutine checkbz
end module checkbz_m
