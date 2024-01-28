!===================================================================================
!
! Routines:
!
! (1) sortbyq Originally By MLT Last Modified: 9/27/2010 (JRD)
!
!===================================================================================
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
module sortbyq_m
  use global_m
  implicit none
  private
  public :: &
    sortbyq
contains
subroutine sortbyq(fqa,qqa,g0a,ifqa,irqa,qg,kg,crys)
  real(DP), intent(out) :: fqa(3,peinf%myown), qqa(peinf%myown)
  integer, intent(out) :: g0a(3,peinf%myown)
  integer, intent(out) :: ifqa(peinf%myown)
  integer, intent(out) :: irqa(peinf%myown)
  type (grid), intent(in) :: qg, kg
  type (crystal), intent(in) :: crys
  integer :: sinv,ik,ic,iv,ikp,icp,ivp,ii,ifq,irq,ifqmin
  real(DP) :: fq(3),difmin,diff(3),length
 
  do ii = 1, peinf%myown
    ik=peinf%ik(peinf%inode+1,ii)
    ic=peinf%ic(peinf%inode+1,ii)
    iv=peinf%iv(peinf%inode+1,ii)
    ikp=peinf%ikp(peinf%inode+1,ii)
    icp=peinf%icp(peinf%inode+1,ii)
    ivp=peinf%ivp(peinf%inode+1,ii)
! If ik.gt.ikp, define q + g0 = kp-k
! all factors of (G0-G) and (q) must get a minus sign,
! and epsinv must be replaced by its complex conjugate
! otherwise, work with q + g0 = k-kp
! sinv keeps track of signs (in the comments below, assume ik.le.ikp)
    sinv = 1
    if (ik.gt.ikp) sinv = -1
    fq(:)=dble(sinv)*( kg%f(:,ik)-kg%f(:,ikp) )
! For a given pair of k-points (ik,ikp), calculate the shortest
! q = k - kp + G_0 , with G_0 a reciprocal lattice vector
    difmin=1.d10
    ifqmin=0
    do ifq=1,qg%nf
      diff(1:3) = fq(1:3) - qg%f(1:3, ifq)
      diff(1:3) = diff(1:3) - anint( diff(1:3) )
      length = DOT_PRODUCT(diff,MATMUL(crys%bdot,diff))
! Found a good matching point. Save it
      if (length <= TOL_Zero) then
        ifqmin = ifq
        difmin = 0.d0
        exit
      endif
! Check other points for the best matching point
      if (length.lt.difmin) then
        difmin = length
        ifqmin = ifq
      endif
    enddo ! ifq
    ifq=ifqmin
    if (difmin.gt.TOL_Zero .or. ifq .eq. 0) then
      write(0,'(a,3f12.6,a,i6,a,i6,a)') 'epsmat does not contain q-point ', fq, ' (= k ', ik, ' - k ', ikp, ')'
      write(0,'(a,i6,a,f10.5)') 'Closest q-point is', ifq, ' of distance ', difmin
      call die("Needed q-point not found. This likely means your WFN_co file is shifted, which is not recommended.")
    endif
    fqa(:,ii) = qg%f(:,ifq)
    g0a(:,ii) = anint(fq(:) - fqa(:,ii))
    ifqa(ii)=ifq
    irqa(ii)=qg%indr(ifq)
    qqa(ii) = sqrt( DOT_PRODUCT(fqa(:,ii),MATMUL(crys%bdot,fqa(:,ii))) )
  enddo ! ii
  allocate(peinf%nxqown (qg%nr))
  allocate(peinf%nxqi (peinf%myown))
  peinf%nxqown = 0
  peinf%nxqi = 0
  do irq = 1, qg%nr
    peinf%nxqown(irq)=0
    do ii = 1, peinf%myown
      if (irqa(ii) .eq. irq) then
        peinf%nxqown(irq)=peinf%nxqown(irq)+1
        peinf%nxqi(peinf%nxqown(irq))=ii
      endif
    enddo
  enddo
 
  return
end subroutine sortbyq
end module sortbyq_m
