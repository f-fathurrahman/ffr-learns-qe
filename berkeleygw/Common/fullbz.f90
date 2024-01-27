!-----------------------------------------------------------------------
!
! module fullbz_m
!
! Routines: fullbz
!
! If wigner_seitz = .true. (for BSE, PlotXct, NonLinearOptics)
! Uses a Wigner-Seitz construction to define the Brillouin zone.
! If wigner_seitz = .false. (for Epsilon, Sigma)
! Uses the usual "box" BZ.
!
! input: crys,syms type
! ntran (number of symmetry operations)
! gr%nr
! gr%rk
!
! output: gr type (except gr%nr, gr%rk)
!
!-----------------------------------------------------------------------
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
module fullbz_m
  use global_m
  use misc_m
  implicit none
  private
  public :: fullbz, dealloc_grid
contains
  !> Set up the mapping between the irreducible and the full full Brillouin zone
  !! i.e. the k-point correspondence and the symmetry operations.
  !! Everything is stored into the grid object.
  subroutine fullbz(crys,syms,gr,ntran,skip_checkbz,wigner_seitz,paranoid,do_nothing,nfix)
    type (crystal), intent(in) :: crys
    type (symmetry), intent(in) :: syms
    type (grid), intent(inout) :: gr
    integer, intent(in) :: ntran !< number of sym ops, typically syms%ntrans
    logical, intent(out) :: skip_checkbz
    logical, intent(in) :: wigner_seitz !< do a Wigner-Seitz construction
    logical, intent(in) :: paranoid !< perform paranoia check: no k-pts differ by G-vector
    logical, optional, intent(in) :: do_nothing !< initialize without performing any checks.
    integer, optional, intent(in) :: nfix !< don`t unfold the first nfix points
    integer :: ii,jj,kk,it,ir,if,i1,i2,i3,gpt(3),iostat_c,ntran_
    real(DP) :: tmpf(3),tmpfm(3),length,lmin,fq(3)
    real(DP), allocatable :: fk(:,:),kg0(:,:)
    integer, allocatable :: indr(:),itran(:)
    logical :: found
!
! i r,nr k-points in irr bz
! o fk,gr%nf k-points in full bz
! o sz radius of a spherical subzone equivalent to
! one point in the set fk
!
! Loop over ir-points
!
   
    if(present(do_nothing)) then
      if(do_nothing) then
        if(wigner_seitz .or. paranoid) call die("bug: cannot call fullbz with do_nothing and wigner_seitz or paranoid.")
        gr%nf=gr%nr
        gr%sz=2.0d0*PI_D*(3.0d0/(4.0d0*PI_D*gr%nf*crys%celvol))**(1.0d0/3.0d0)
        allocate(gr%kg0 (3,gr%nf))
        allocate(gr%f (3,gr%nf))
        allocate(gr%itran (gr%nf))
        allocate(gr%indr (gr%nf))
        gr%kg0(:,:)=0
        gr%itran(:)=1
        gr%f(:,:)=gr%r(:,:)
        do ii=1,gr%nf
          gr%indr(ii)=ii
        enddo
        skip_checkbz = .true.
       
        return
      endif
    endif
    allocate(fk (3,gr%nr*ntran))
    allocate(indr (gr%nr*ntran))
    allocate(itran (gr%nr*ntran))
    call open_file(unit=21,file='fullbz.dat',status='old',iostat=iostat_c)
    skip_checkbz = (iostat_c==0)
    if (skip_checkbz) then
      write(6,'(3x,a)') 'Reading the full Brillouin Zone from fullbz.dat'
      write(6,'(3x,a)') 'Will not be checking if the points there form a full zone'
      read(21,*) gr%nf
      do ii=1,gr%nf
        read(21,*) fk(1:3,ii),itran(ii),indr(ii)
      enddo
      call close_file(21)
      ! This is just in case we want to do wigner_seitz or paranoid
    else
      gr%nf=0
      do ir=1,gr%nr
!
! Loop over transformations
!
        ntran_ = ntran
        if (present(nfix)) then
          if (ir<=nfix) ntran_ = 1
        endif
        do it=1,ntran_
!
! Rotate gr%r and put into tmpf
!
          tmpf(:) = matmul(dble(syms%mtrx(:,:,it)),gr%r(:,ir))
          call k_range(tmpf, gpt, TOL_Small)
!
! Compare to other points in full zone
!
          found = .false.
          do if=1,gr%nf
            if(all(abs(tmpf(1:3)-fk(1:3,if)).lt.TOL_Small)) then
              found = .true.
              exit
            endif
          enddo
          if(found) cycle ! skip adding it
!
! Store new kpoint in fbz
!
          gr%nf=gr%nf+1
          if (gr%nf > gr%nr * ntran) call die('fullbz internal error')
          fk(1:3,gr%nf)=tmpf(1:3)
!
! Store index of rotation itran and corresponding IBZ point
!
          itran(gr%nf)=it
          indr(gr%nf)=ir
        enddo !end loop over symmetries
      enddo !end loop over the q-points from epsmat and eps0mat
    endif ! whether reading from fullbz.dat
!
! SIB: Now that we have the full BZ with components [0,1), we
! will move each point into the Wigner-Seitz cell by adding G-vectors
! to it until we get the shortest length vector we can. Then we find
! the appropriate Umklapp vector as well (kg0).
!
    if(wigner_seitz) then
      allocate(kg0 (3,gr%nr*ntran))
      do ii=1,gr%nf
        tmpf(:) = fk(:,ii)
        lmin = 1.0d10
        do i1=-ncell+1,ncell
          do i2=-ncell+1,ncell
            do i3=-ncell+1,ncell
              fq(1) = tmpf(1) - i1
              fq(2) = tmpf(2) - i2
              fq(3) = tmpf(3) - i3
              length = DOT_PRODUCT(fq,MATMUL(crys%bdot,fq))
              if (length.lt.lmin) then
                lmin = length
                tmpfm(:) = fq(:)
              endif
            enddo
          enddo
        enddo
        fk(:,ii) = tmpfm
!
! SIB: Find Umklapp (kg0) a.k.a. translation
!
        tmpf = MATMUL(dble(syms%mtrx(:,:,itran(ii))),gr%r(:,indr(ii)))
        tmpf(:) = fk(:,ii) - tmpf
        do jj=1,3
          if (tmpf(jj).ge.0.0) kg0(jj,ii)=tmpf(jj)+TOL_Small
          if (tmpf(jj).lt.0.0) kg0(jj,ii)=tmpf(jj)-TOL_Small
        enddo
      enddo !over full zone (ii)
      allocate(gr%kg0 (3,gr%nf))
      gr%kg0(1:3,1:gr%nf)=kg0(1:3,1:gr%nf)
      if(allocated(kg0))then;deallocate(kg0);endif
    else
      nullify(gr%kg0)
    endif
!------------------------------------------------------------------------
! SIB: Paranoia check. We will compare all pairs of points in the full
! zone to each other in order to see if any of them differ by a G-vector.
! If they do, we are in big trouble and stop!
    if(paranoid) then
      ii_loop: do ii=1,gr%nf
        do jj=1,ii-1
          tmpf = abs(fk(:,ii)-fk(:,jj))
! after this loop, tmpf contains how far each component of fii-fjj is
! from the closest integer to it.
          do kk=1,3
            it = tmpf(kk)
            tmpf(kk) = tmpf(kk)-it
            if (tmpf(kk).ge.0.5) tmpf(kk)=1.0-tmpf(kk)
          enddo
! if fkii-fkjj is very close to an integer vector (G-vector), trouble!
          if (sum(abs(tmpf)).le.TOL_Small) then
            write(0,123) ii,jj,fk(:,ii)-fk(:,jj)
123 format('equiv kpts',i4,' and ',i4,' diff=',3f10.5)
            call die('equivalent points found in the full BZ')
          endif
        enddo
      enddo ii_loop
    endif
!------------------------------------------------------------------------
! SIB: Now store all gathered information where it belongs
    allocate(gr%f (3,gr%nf))
    allocate(gr%indr (gr%nf))
    allocate(gr%itran (gr%nf))
    gr%f(1:3,1:gr%nf)=fk(1:3,1:gr%nf)
    gr%indr(1:gr%nf)=indr(1:gr%nf)
    gr%itran(1:gr%nf)=itran(1:gr%nf)
    if(allocated(fk))then;deallocate(fk);endif
    if(allocated(indr))then;deallocate(indr);endif
    if(allocated(itran))then;deallocate(itran);endif
! Compute radius of spherical subzone
! assuming bz filled with gr%nf spheres
    gr%sz=2.0d0*PI_D*(3.0d0/(4.0d0*PI_D*gr%nf*crys%celvol))**(1.0d0/3.0d0)
    if (ntran.eq.1 .and. gr%nf.ne.gr%nr .and. paranoid) then
      write(0,*) gr%nf, ' and ', gr%nr
      call die('fullbz paranoia check failed')
    endif
   
    return
  end subroutine fullbz
  !-----------------------------------------------------------------------
  subroutine dealloc_grid(gr)
    type(grid), intent(inout) :: gr
   
    if(associated(gr%r))then;deallocate(gr%r);nullify(gr%r);endif
    if(associated(gr%f))then;deallocate(gr%f);nullify(gr%f);endif
    if(associated(gr%indr))then;deallocate(gr%indr);nullify(gr%indr);endif
    if(associated(gr%itran))then;deallocate(gr%itran);nullify(gr%itran);endif
    if(associated(gr%kg0))then;deallocate(gr%kg0);nullify(gr%kg0);endif
   
    return
  end subroutine dealloc_grid
end module fullbz_m
