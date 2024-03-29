!============================================================================
!
! MODULE: fftw_m, originally by DAS 1/14/2011
!
!> Routines used with FFTW, as well as interfaces for library calls.
!
! DESCRIPTION:
!> No FFTW calls should exist outside this routine: the wrapper routines
!! here should be used everywhere.
!!
!! Interfaces for FFTW2 functions are formulated from fftw-2.1.5/fftw/fftwf77.c
!! http://www.fftw.org/fftw2_doc/fftw_5.html. Every FFTW2 function used should
!! have an interface to ensure the argument types are correct.
!! Contains include file fftw_f77.i which has parameters used in FFTW2 calls.
!!
!! FFTW3 provides its own file of constants and interfaces, fftw3.f03.
!! See http://fftw.org/doc/Overview-of-Fortran-interface.html
!
!============================================================================
module fftw_m
  !> use of this is recommended by FFTW3 documentation. It causes no harm for
  !! FFTW2, but some compilers (e.g. Open64) do not have it available, so we
  !! will keep it hidden for FFTW2 so as not to have to solve that problem yet.
  use, intrinsic :: iso_c_binding
  use global_m
  use timing_m, only: timing => common_timing
  implicit none
  private
!> It is better to use this one which has interfaces too, rather than fftw3.f which has only constants
  include 'fftw3.f'
    integer*8, private :: fft_plan = 0
    integer, private :: ifirst = 0
    integer, private :: num_fft_threads = 0
    integer :: iret
    integer, private :: Nfftold(3) = 0
  public :: &
    check_FFT_size, &
    setup_FFT_sizes, &
    gvec_to_fft_index, &
    put_into_fftbox, &
    get_from_fftbox, &
    do_FFT, &
    conjg_fftbox, &
    multiply_fftboxes, &
    destroy_fftw_plans
  interface put_into_fftbox
    module procedure dput_into_fftbox, zput_into_fftbox
  end interface put_into_fftbox
  interface get_from_fftbox
    module procedure dget_from_fftbox, zget_from_fftbox
  end interface get_from_fftbox
contains
!> Originally by gsm Last Modified: 4/10/2010 (gsm)
!! Best FFT grid dimension is given by 2^a*3^b*5^c*7^d*11^e*13^f
!! where a,b,c,d are arbitrary and e,f are 0 or 1
!! Ref: http://www.fftw.org/fftw2_doc/fftw_3.html
!! On entry
!! Nfft = FFT grid dimension to test
!! Nfac = number of factors to test
!! On exit
!! check_FFT_size = .true. if good FFT grid dimension
  logical function check_FFT_size(Nfft, Nfac)
    integer, intent(in) :: Nfft, Nfac
    integer :: remainder, product, ifac, ipow, maxpow
    integer, parameter :: maxfac = 6
    integer :: pow(maxfac)
    integer, parameter :: fac(maxfac) = (/ 2, 3, 5, 7, 11, 13 /)
   
    if(Nfft .lt. 1 .or. Nfac .lt. 1 .or. Nfac .gt. maxfac) then
      call die('check_FFT_size input')
    endif
    remainder = Nfft
    do ifac = 1, maxfac
      pow(ifac) = 0
    enddo
    do ifac = 1, Nfac
      maxpow = int(log(dble(remainder)) / log(dble(fac(ifac)))) + 1
      do ipow = 1, maxpow
        if (mod(remainder, fac(ifac)) .eq. 0) then
          remainder = remainder / fac(ifac)
          pow(ifac) = pow(ifac) + 1
        endif
      enddo
    enddo
    product = remainder
    do ifac = 1, Nfac
      do ipow = 1, pow(ifac)
        product = product * fac(ifac)
      enddo
    enddo
    if (product .ne. Nfft) then
      call die('Internal error in check_FFT_size; factorization failed')
    endif
    check_FFT_size = remainder .eq. 1 .and. pow(5) .le. 1 .and. pow(6) .le. 1
   
    return
  end function check_FFT_size
!> The former "fft_routines.f90"
!! Sohrab Ismail-Beigi Feb 28 2001
!!
!! There are a set of Fast Fourier-related routines that are used
!! to compute the matrix elements of the type <nk|e^(i*G.r)|mk`>.
!! For many G-vectors, FFTs will be the fastest way to compute them.
!!
!! The FFTW (http://www.fftw.org) suite of routines do the actual work.
!! Most of what is below is interfacing code and routines that simplify
!! small and useful tasks.
 !
!!
!! Given gvec%FFTgrid(1:3) values (in FFTgrid), finds appropriate FFT box
!! sizes to use in Nfft(1:3). scale = 1/(Nfftx*Nffty*Nfftz).
!!
  subroutine setup_FFT_sizes(FFTgrid,Nfft,scale)
    integer, intent(in) :: FFTgrid(3)
    integer, intent(out) :: Nfft(3)
    real(DP), intent(out) :: scale
    integer, parameter :: Nfac = 3
    integer :: i
   
    do i=1,3
      Nfft(i) = FFTgrid(i)
      do while (.not. check_FFT_size(Nfft(i), Nfac))
        Nfft(i) = Nfft(i) + 1
      enddo
    enddo
    scale = 1.0d0/product(Nfft(1:3))
   
    return
  end subroutine setup_FFT_sizes
!> Takes the G-vector g(1:3) and FFT box size Nfft(1:3) and finds the
!! point idx(1:3) in the box corresponding to that G-vector.
!!
  subroutine gvec_to_fft_index(g,idx,Nfft)
    integer, intent(in) :: g(3), Nfft(3)
    integer, intent(out) :: idx(3)
! no push/pop since called too frequently.
    idx(1:3) = g(1:3) + 1
    if (g(1) < 0) idx(1) = Nfft(1) + idx(1)
    if (g(2) < 0) idx(2) = Nfft(2) + idx(2)
    if (g(3) < 0) idx(3) = Nfft(3) + idx(3)
    return
  end subroutine gvec_to_fft_index
!> Do an FFT on the fftbox in place: destroys contents of fftbox
!! and replaces them by the Fourier transform.
!!
!! The FFT done is:
!!
!! fftbox(p) <- sum_j { fftbox(j)*e^{sign*i*j.p} }
!!
!! where j and p are integer 3-vectors ranging over Nfft(1:3).
!!
  subroutine do_FFT(fftbox, Nfft, sign)
    complex(DPC), intent(inout) :: fftbox(:,:,:)
    integer, intent(in) :: Nfft(3)
    integer, intent(in) :: sign
    character(len=100) :: str
   
!JRD To be removed
    !complex(DPC), allocatable :: fftbox2(:,:,:)
    if (peinf%verb_max) then
      write(str,'(a,2(i0," x "),i0,a)') 'Creating ', Nfft(1:3), ' FFTW plans.'
      call logit(str)
    endif
    call timing%start(timing%fft_plan)
    if (sign == 1) then
      call dfftw_plan_dft(fft_plan,3,Nfft,fftbox,fftbox,FFTW_BACKWARD,FFTW_ESTIMATE)
      !call dfftw_plan_dft(fft_plan,3,Nfft,fftbox,fftbox2,FFTW_BACKWARD,FFTW_MEASURE)
      !call dfftw_plan_dft_3d(fft_plan,Nfft(1),Nfft(2),Nfft(3),fftbox,fftbox2,FFTW_BACKWARD,FFTW_MEASURE)
    else if (sign == -1) then
      call dfftw_plan_dft(fft_plan,3,Nfft,fftbox,fftbox,FFTW_FORWARD,FFTW_ESTIMATE)
      !call dfftw_plan_dft(fft_plan,3,Nfft,fftbox,fftbox2,FFTW_FORWARD,FFTW_MEASURE)
      !call dfftw_plan_dft_3d(fft_plan,Nfft(1),Nfft(2),Nfft(3),fftbox,fftbox2,FFTW_FORWARD,FFTW_MEASURE)
    else
      call die('sign is not 1 or -1 in do_FFT')
    endif
    call timing%stop(timing%fft_plan)
    call timing%start(timing%fft_exec)
    call dfftw_execute_dft(fft_plan,fftbox,fftbox)
    call timing%stop(timing%fft_exec)
    ! otherwise there is a memory leak
    call dfftw_destroy_plan(fft_plan)
    !call dfftw_cleanup_threads()
    Nfftold(:) = -1
!JRD To be removed
    !fftbox(:,:,:)=fftbox2(:,:,:)
    !if(allocated(fftbox2))then;deallocate(fftbox2);endif
   
    return
  end subroutine do_FFT
  subroutine destroy_fftw_plans()
    character(len=100) :: str
    ! FFTW plan was never created
    if(all(Nfftold(1:3) == 0)) return
   
    if(all(Nfftold(1:3) == -1)) then
     
      return
      ! call die("Cannot destroy FFTW plan for a second time.")
    endif
    write(str,'(a,2(i0," x "),i0,a)') 'Destroying ', Nfftold(1:3), ' FFTW plans.'
    call logit(str)
    Nfftold(1:3) = -1 ! make clear there is no plan anymore so we do not try to destroy twice
    call dfftw_destroy_plan(fft_plan)
    ! should forget wisdom here, but I cannot figure out how... --DAS
   
    return
  end subroutine destroy_fftw_plans
!> Complex conjugate contents of FFT box
!
  subroutine conjg_fftbox(fftbox,Nfft)
    integer, intent(in) :: Nfft(3)
    complex(DPC), intent(inout) :: fftbox(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
    ! for some reason, absoft segfaults if dims specified for fftbox as above right
    integer :: ix,iy,iz
   
    call timing%start(timing%fft_conjg)
!disabled PARALLEL DO PRIVATE (ix,iy,iz)
    do iz = 1, Nfft(3)
    do iy = 1, Nfft(2)
    do ix = 1, Nfft(1)
      fftbox(ix,iy,iz) = conjg(fftbox(ix,iy,iz))
    enddo
    enddo
    enddo
!disabled END PARALLEL DO
    call timing%stop(timing%fft_conjg)
   
    return
  end subroutine conjg_fftbox
!> Multiply contents of two fft boxes, result into fftbox2
!
  subroutine multiply_fftboxes(fftbox1, fftbox2, Nfft)
    integer, intent(in) :: Nfft(3)
    complex(DPC), intent(in) :: fftbox1(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
    complex(DPC), intent(inout) :: fftbox2(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
    integer :: ix,iy,iz
   
    call timing%start(timing%fft_mltply)
    !forall(ix=1:Nfft(1), iy=1:Nfft(2), iz=1:Nfft(3)) &
    ! fftbox2(ix,iy,iz) = fftbox1(ix,iy,iz) * fftbox2(ix,iy,iz)
!disabled PARALLEL DO PRIVATE (ix,iy,iz)
    do iz = 1, Nfft(3)
    do iy = 1, Nfft(2)
    do ix = 1, Nfft(1)
      fftbox2(ix,iy,iz) = fftbox1(ix,iy,iz) * fftbox2(ix,iy,iz)
    enddo
    enddo
    enddo
!disabled END PARALLEL DO
    call timing%stop(timing%fft_mltply)
   
    return
  end subroutine multiply_fftboxes
! use between inclusions of f_defs.h in template modules
! list here everything defined differently by flavor in f_defs.h
! these undefs prevent lots of warnings from cpp
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
!overrules flavor.mk
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
!============================================================================
!
! Included from fftw.F90
!
!============================================================================
!> This routine takes data(1:ndata) and puts it into the FFT box fftbox(:,:,:).
!! The FFT box is zeroed out first, and the data is entered into it.
!!
!! ndata -- number of data items in data(:)
!! data -- the data set, real or complex, depending on ifdef CPLX
!! ng -- number of g vectors in glist
!! glist -- a master list of g vectors
!! gindex(1:ng) -- which g vector (in the master list) the data(1:ndata)
!! actually refer to: so data(j) is for the g-vector
!! glist(1:3,gindex(j))
!! fftbox(:,:,:) -- 3D complex FFT box where the data is put
!! Nfft(1:3) -- sizes of FFT box Nx,Ny,Nz
subroutine dput_into_fftbox(ndata, data, glist, gindex, fftbox, Nfft)
  use timing_m, only: timing => common_timing
  integer, intent(in) :: ndata
  real(DP), intent(in) :: data(:) !< (ndata) this is to avoid creation of array temporary
  integer, intent(in) :: glist(:,:) !< (3, ng)
  integer, intent(in) :: gindex(:) !< (ng)
  integer, intent(in) :: Nfft(:) !< (3)
  complex(DPC), intent(out) :: fftbox(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
  integer :: j, k, bidx(3)
 
  ! Zero out FFT box and put data into it
  call timing%start(timing%fft_zero)
!disabled PARALLEL
!disabled DO COLLAPSE(2)
  do j=1,Nfft(3)
    do k=1,Nfft(2)
      fftbox(:,k,j) = (0.0d0,0.0d0)
    enddo
  enddo
!disabled END DO
!disabled END PARALLEL
  call timing%stop(timing%fft_zero)
  call timing%start(timing%fft_put)
!disabled PARALLEL PRIVATE(bidx,j) SHARED(fftbox, glist, gindex, data)
!disabled DO
  do j=1,ndata
    call gvec_to_fft_index(glist(:,gindex(j)),bidx,Nfft)
    fftbox(bidx(1),bidx(2),bidx(3)) = data(j)
  end do
!disabled END DO
!disabled END PARALLEL
  call timing%stop(timing%fft_put)
 
  return
end subroutine dput_into_fftbox
!> Does the inverse of the above routine: takes the data in the
!! fftbox(:,:,:) and puts it into the data(1:ndata) array. ndata entries
!! are extracted, and the gindex and glist specify which ones to get:
!! data(j) corresponds to the g-vector glist(:,gindex(j)). The data
!! in fftbox is multiplied by scale before storage into data(:).
!!
!! data(:) is zeroed first and then the data is put into it.
!!
subroutine dget_from_fftbox(ndata, data, glist, gindex, fftbox, Nfft, scale)
  use timing_m, only: timing => common_timing
  integer, intent(in) :: ndata
  real(DP), intent(out) :: data(:) !< (ndata)
  integer, intent(in) :: glist(:,:) !< (3, ng)
  integer, intent(in) :: gindex(:) !< (ng)
  integer, intent(in) :: Nfft(:)
  complex(DPC), intent(in) :: fftbox(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
  real(DP), intent(in) :: scale
  integer :: j, bidx(3)
 
  call timing%start(timing%fft_get)
!disabled PARALLEL PRIVATE(bidx)
!disabled DO
  do j=1,ndata
    !data(j) = 0.0
    call gvec_to_fft_index(glist(:,gindex(j)),bidx,Nfft)
    data(j) = fftbox(bidx(1),bidx(2),bidx(3))*scale
  end do
!disabled END DO
!disabled END PARALLEL
  call timing%stop(timing%fft_get)
 
  return
end subroutine dget_from_fftbox
! use between inclusions of f_defs.h in template modules
! list here everything defined differently by flavor in f_defs.h
! these undefs prevent lots of warnings from cpp
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
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
!============================================================================
!
! Included from fftw.F90
!
!============================================================================
!> This routine takes data(1:ndata) and puts it into the FFT box fftbox(:,:,:).
!! The FFT box is zeroed out first, and the data is entered into it.
!!
!! ndata -- number of data items in data(:)
!! data -- the data set, real or complex, depending on ifdef
!! ng -- number of g vectors in glist
!! glist -- a master list of g vectors
!! gindex(1:ng) -- which g vector (in the master list) the data(1:ndata)
!! actually refer to: so data(j) is for the g-vector
!! glist(1:3,gindex(j))
!! fftbox(:,:,:) -- 3D complex FFT box where the data is put
!! Nfft(1:3) -- sizes of FFT box Nx,Ny,Nz
subroutine zput_into_fftbox(ndata, data, glist, gindex, fftbox, Nfft)
  use timing_m, only: timing => common_timing
  integer, intent(in) :: ndata
  complex(DPC), intent(in) :: data(:) !< (ndata) this is to avoid creation of array temporary
  integer, intent(in) :: glist(:,:) !< (3, ng)
  integer, intent(in) :: gindex(:) !< (ng)
  integer, intent(in) :: Nfft(:) !< (3)
  complex(DPC), intent(out) :: fftbox(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
  integer :: j, k, bidx(3)
 
  ! Zero out FFT box and put data into it
  call timing%start(timing%fft_zero)
!disabled PARALLEL
!disabled DO COLLAPSE(2)
  do j=1,Nfft(3)
    do k=1,Nfft(2)
      fftbox(:,k,j) = (0.0d0,0.0d0)
    enddo
  enddo
!disabled END DO
!disabled END PARALLEL
  call timing%stop(timing%fft_zero)
  call timing%start(timing%fft_put)
!disabled PARALLEL PRIVATE(bidx,j) SHARED(fftbox, glist, gindex, data)
!disabled DO
  do j=1,ndata
    call gvec_to_fft_index(glist(:,gindex(j)),bidx,Nfft)
    fftbox(bidx(1),bidx(2),bidx(3)) = data(j)
  end do
!disabled END DO
!disabled END PARALLEL
  call timing%stop(timing%fft_put)
 
  return
end subroutine zput_into_fftbox
!> Does the inverse of the above routine: takes the data in the
!! fftbox(:,:,:) and puts it into the data(1:ndata) array. ndata entries
!! are extracted, and the gindex and glist specify which ones to get:
!! data(j) corresponds to the g-vector glist(:,gindex(j)). The data
!! in fftbox is multiplied by scale before storage into data(:).
!!
!! data(:) is zeroed first and then the data is put into it.
!!
subroutine zget_from_fftbox(ndata, data, glist, gindex, fftbox, Nfft, scale)
  use timing_m, only: timing => common_timing
  integer, intent(in) :: ndata
  complex(DPC), intent(out) :: data(:) !< (ndata)
  integer, intent(in) :: glist(:,:) !< (3, ng)
  integer, intent(in) :: gindex(:) !< (ng)
  integer, intent(in) :: Nfft(:)
  complex(DPC), intent(in) :: fftbox(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
  real(DP), intent(in) :: scale
  integer :: j, bidx(3)
 
  call timing%start(timing%fft_get)
!disabled PARALLEL PRIVATE(bidx)
!disabled DO
  do j=1,ndata
    !data(j) = 0.0
    call gvec_to_fft_index(glist(:,gindex(j)),bidx,Nfft)
    data(j) = fftbox(bidx(1),bidx(2),bidx(3))*scale
  end do
!disabled END DO
!disabled END PARALLEL
  call timing%stop(timing%fft_get)
 
  return
end subroutine zget_from_fftbox
end module fftw_m
