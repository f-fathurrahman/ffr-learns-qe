!===============================================================================
!
! Routines:
!
! (1) epsmodel() Originally By MLT Last Modified 7/3/2008 (JRD)
!
! Model dielectric function for GaAs
! see G. Cappellini et al, PRB 47, 9892 (1993)
!
! (2) epsmodel_pac() Orignally By MLT Last Modified 7/3/2008 (JRD)
!
! Model dielectric function for polyacetylene
! see B. Wenzien et al, PRB 51, 14701 (1995)
! omega_plasma = 2.0 (actual value)
! q_ThomasFermi = 1.0
! epsp : dielectric constant perpendicular to chain
! epsc : dielectric constant along chain
!
! (3) epsmodel_cnt() Orignally By Li Yang Last Modified 7/3/2008 (JRD)
!
! Model dielectric function for carbon nanotube
! zcut=inf, rho=15.11, grid 1x1x64
!
! These routines are not used in current version of the code.
! Calculated dielectric function from Epsilon is used instead.
!
!===============================================================================
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
module epsmodel_m
  use global_m
  implicit none
  public :: &
    epsmodel, &
    epsmodel_pac, &
    epsmodel_cnt
contains
subroutine epsmodel(bdot,fq,epsinv)
  real(DP), intent(in) :: bdot(3,3),fq(3)
  real(DP), intent(out) :: epsinv
  real(DP) :: eps0,alpha,beta,q2
 
  q2= DOT_PRODUCT(fq,MATMUL(bdot,fq))
! eps0 = 0.0747
! alpha = 1.37
! beta = 0.76
  eps0 = 0.0971
  alpha = 1.261
  beta = 0.67
  epsinv = 1.d0/(1.d0 + 1.d0/(eps0 + alpha*q2 + beta*q2*q2))
 
  return
end subroutine epsmodel
!==============================================================================
subroutine epsmodel_pac(bdot,fq,epsinv)
  real(DP), intent(in) :: bdot(3,3),fq(3)
  real(DP), intent(out) :: epsinv
  ! epsinv was complex originally, but the calc here seems to give always a real result
  integer :: jj
  real(DP) :: epsp,epsc,alphac,alphap,sum,qc,qp,cost2,sint2,qq
 
  epsp=3.7
  epsc=20.0
  alphac=8.0
  alphap=1.0
  qq= sqrt( DOT_PRODUCT(fq,MATMUL(bdot,fq)) )
  if (abs(qq).lt.1.d-5) then
    epsinv=0.05
  else
    sum=0.0
    do jj=1,3
      sum = sum + fq(jj)*bdot(jj,3)
    enddo
    qc=sum/sqrt(bdot(3,3))
    qp=sqrt(abs(qq**2 - qc**2))
    cost2=(qc/qq)**2
    sint2=(qp/qq)**2
    epsinv = cost2/(1.d0+1.d0/(1.d0/(epsc-1)+alphac*qc*qc+ &
      qq**4/8.d0)) + sint2/(1.d0+1.d0/(1.d0/(epsp-1d0)+ &
      alphap*qp*qp+qq**4/8.d0))
  endif
 
  return
end subroutine epsmodel_pac
!==============================================================================
subroutine epsmodel_cnt(diff,epsinv)
  real(DP), intent(in) :: diff(3)
  real(DP), intent(out) :: epsinv
  real(DP) :: a0,a1,a2,x,length2,qz2,qp2,epsz,epsp,v
 
  a0 = -248.187d0
  a1 = 921.851d0
  a2 = 347.231d0
  x=abs(diff(3))
  if(x.lt.1.0d-15) then
    epsinv=1.d0
  else
    epsinv = (1+a2*x*x)/(1+a0*x*x*log(x*0.7889*15.11)+a1*x*x)
  endif
 
  return
end subroutine epsmodel_cnt
end module epsmodel_m
