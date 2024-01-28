!============================================================================
!
! Utilities:
!
! (1) bsebinasc() Originally By MLT Last Modified 7/1/2008
!
! File converter. Converts files "bsedmat","bsexmat","dtmat",
! "vmtxel","eps2_moments","eigenvectors"
! from binary to ascii. If a file is unavailable, just skip it.
!
! input: [filename][ext_in]
! output: [filename][ext_out]
!
! Default for ext_in is "".
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
program bsebinasc
  use global_m
  use bse_convert_m
  implicit none
  character*18 :: ext_in,ext_out
  integer :: nargs, err
  nargs = command_argument_count()
  if(nargs == 1) then
    ext_in = ""
    call get_command_argument(1, ext_out)
  else if(nargs == 2) then
    call get_command_argument(1,ext_in)
    call get_command_argument(2,ext_out)
  else
    call die('Usage: bsebinasc [ext_in] ext_out')
  endif
!-----------------------------------
! Convert bsedmat
  write(6,*) ' Converting bsedmat',TRUNC(ext_in),' -> bsedmat', TRUNC(ext_out)
  call open_file(10,file='bsedmat'//TRUNC(ext_in),form='unformatted', status='old',iostat=err)
  if (err.ne.0) then
    write(0,*) 'WARNING: could not open bsedmat',TRUNC(ext_in), '. Skipping.'
  else
    call open_file(11,file='bsedmat'//TRUNC(ext_out),form='formatted',status='replace')
    call bsemat_binasc(10, 11)
    call close_file(11)
    call close_file(10)
  endif
!---------------------------------------
! Convert bsexmat
  write(6,*) ' Converting bsexmat',TRUNC(ext_in),' -> bsexmat', TRUNC(ext_out)
  call open_file(20,file='bsexmat'//TRUNC(ext_in),form='unformatted', status='old',iostat=err)
  if (err.ne.0) then
    write(0,*) 'WARNING: could not open bsexmat',TRUNC(ext_in), '. Skipping.'
  else
    call open_file(11,file='bsexmat'//TRUNC(ext_out),form='formatted',status='replace')
    call bsemat_binasc(20, 11)
    call close_file(11)
    call close_file(20)
  endif
!---------------------------------
! Convert dtmat
  write(6,*) ' Converting dtmat',TRUNC(ext_in),' -> dtmat', TRUNC(ext_out)
  call open_file(30,file='dtmat'//TRUNC(ext_in),form='unformatted',status='old',iostat=err)
  if (err.ne.0) then
    write(0,*) 'WARNING: could not open dtmat',TRUNC(ext_in), '. Skipping.'
  else
    call open_file(31,file='dtmat'//TRUNC(ext_out),form='formatted',status='replace')
    call dtmat_binasc(30, 31)
    call close_file(31)
    call close_file(30)
  endif
!----------------------------------
! Convert vmtxel
  write(6,*) ' Converting vmtxel',TRUNC(ext_in),' -> vmtxel', TRUNC(ext_out)
  call open_file(40,file='vmtxel'//TRUNC(ext_in),form='unformatted',status='old',iostat=err)
  if (err.ne.0) then
    write(0,*) 'WARNING: could not open vmtxel',TRUNC(ext_in), '. Skipping.'
  else
    call open_file(41,file='vmtxel'//TRUNC(ext_out),form='formatted',status='replace')
    call vmtxel_binasc(40, 41)
    call close_file(40)
    call close_file(41)
  endif
!----------------------------------
! Convert eps2_moments
  write(6,*) ' Converting eps2_moments',TRUNC(ext_in), ' -> eps2_moments',TRUNC(ext_out)
  call open_file(50,file='eps2_moments'//TRUNC(ext_in),form='unformatted', status='old',iostat=err)
  if (err.ne.0) then
    write(0,*) 'WARNING: could not open eps2_moments',TRUNC(ext_in), '. Skipping.'
  else
    call open_file(51,file='eps2_moments'//TRUNC(ext_out),form='formatted',status='replace')
    call eps2_moments_binasc(50, 51)
    call close_file(50)
    call close_file(51)
  endif
!-----------------------------
! Convert eigenvectors
  write(6,*) ' Converting eigenvectors',TRUNC(ext_in),' -> eigenvectors',TRUNC(ext_out)
  call open_file(60,file='eigenvectors'//TRUNC(ext_in),form='unformatted',status='old',iostat=err)
  if (err.ne.0) then
    write(0,*) 'WARNING: could not open eigenvectors',TRUNC(ext_in),'. Skipping.'
  else
    call open_file(61,file='eigenvectors'//TRUNC(ext_out),form='formatted',status='replace')
    call eigenvectors_binasc(60, 61)
    call close_file(60)
    call close_file(61)
  endif
  write(6,*) ' Done '
end program bsebinasc
