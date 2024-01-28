!=================================================================================
!
! Routines:
!
! (1) w_sum() Originally By JRD Last Modified 4/1/2012 (JRD)
!
! Multiply Valence-Valence matrix elements by W to create temporary arrays
! for the head, wing and body.
!
! This routine scales as N^3, but is nested within the the igp loop
! in mtxel_kernel. Thus, it represents an N^4 step. Doing the multiplication
! here is cheaper than doing it in the N^5 g_sum subroutine.
!
!=================================================================================
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
module w_sum_m
  use global_m
  use blas_m
  use algos_kernel_m
  implicit none
  public :: w_sum_cpu, w_sum_openacc
  private
contains
  subroutine w_sum_cpu(xct,wptcol,ofs1,ofs1p,n1,n1p,temph,tempw,tempb,m11p,indinvigp,ng_eps)
    type (xctinfo), intent(in) :: xct
    real(DP), intent(in) :: wptcol(:)
    !> offset (i.e., add this number to map a local index to the global band index)
    integer, intent(in) :: ofs1, ofs1p
    !> number of bands for each wfn
    integer, intent(in) :: n1, n1p
    real(DP), intent(inout) :: tempw(:,:,:,:), tempb(:,:,:,:), temph(:,:,:), m11p(:,:,:,:)
    integer, intent(in) :: indinvigp
    integer, intent(in) :: ng_eps
    real(DP), allocatable :: m11p_conj(:,:,:)
    integer :: isv, ig, i1, i1p, gi1, gi1p
    logical :: use_omp
   
    ! NAG gives me a runtime error with threading here. Did not have time
    ! yet to find out if the error is real or a compiler bug.
    use_omp = .true.
    ! JRD: We allocate a new temporary array in order to get better cache performance
    allocate(m11p_conj (n1,n1p,xct%nspin))
    m11p_conj(:,:,:) = (m11p(indinvigp,:,:,:))
    do isv=1,xct%nspin
      if (indinvigp .eq. 1) then
        temph(ofs1+1:ofs1+n1, ofs1p+1:ofs1p+n1p, isv) = wptcol(1)*m11p_conj(1:n1, 1:n1p, isv)
        !disabled PARALLEL PRIVATE(i1p, gi1p, i1, gi1, ig) DEFAULT(SHARED) IF(use_omp)
        do i1p = 1, n1p
          gi1p = ofs1p+i1p
          do i1 = 1, n1
            gi1 = ofs1+i1
            !disabled DO
            do ig=2,ng_eps
              tempw(ig, gi1, gi1p, isv) = wptcol(ig) * m11p_conj(i1, i1p, isv)
            enddo
            !disabled END DO NOWAIT
          enddo
        enddo
        !disabled END PARALLEL
      else
        tempw(1, ofs1+1:ofs1+n1, ofs1p+1:ofs1p+n1p, isv) = tempw(1, ofs1+1:ofs1+n1, ofs1p+1:ofs1p+n1p, isv) + &
          wptcol(1)*m11p_conj(1:n1, 1:n1p, isv)
        !disabled PARALLEL PRIVATE(i1p, gi1p, i1, gi1, ig) DEFAULT(SHARED) IF(use_omp)
        do i1p = 1, n1p
          gi1p = ofs1p+i1p
          do i1 = 1, n1
            gi1 = ofs1+i1
            !disabled DO
            do ig=2,ng_eps
              tempb(ig, gi1, gi1p, isv) = tempb(ig, gi1, gi1p, isv) + &
                wptcol(ig) * m11p_conj(i1, i1p, isv)
            enddo
            !disabled END DO NOWAIT
          enddo
        enddo
        !disabled END PARALLEL
      endif
    enddo
    if(allocated(m11p_conj))then;deallocate(m11p_conj);endif
   
  end subroutine w_sum_cpu
  subroutine w_sum_openacc(xct, wptcol_mat, m11p_all, temph, tempw, tempb, indinv, ng_eps, &
                       phinv, nepsmin, ngpown_ipe, invband, nbands, ipe)
   type (xctinfo), intent(in) :: xct
   real(DP), intent(in) :: wptcol_mat(:,:)
   real(DP), intent(in), target :: m11p_all(:,:,:,:,:)
   real(DP), intent(inout) :: tempw(:,:,:,:), tempb(:,:,:,:), temph(:,:,:)
   integer, intent(in) :: indinv(:), nbands(:)
   real(DP), intent(in) :: phinv(:)
   integer, intent(in) :: ng_eps, ngpown_ipe, nepsmin, invband, ipe
   
    call die("OpenACC version of w_sum requested, but OpenACC not compiled&
            & into this executable", only_root_writes = .true.)
   
  end subroutine w_sum_openacc
end module w_sum_m
