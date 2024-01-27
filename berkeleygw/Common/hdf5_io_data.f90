! File hdf5_io_data.f90 automatically created from hdf5_io_data.f90p by mako_preprocess.py.
! Do not edit the resulting file (hdf5_io_data.f90) directly!
!>=========================================================================
!!
!! Module:
!!
!! hdf5_io_data_m Originally by FHJ Last Modified 08/2019 (FHJ)
!!
!! High-level routines to read and write data in HDF5 format.
!! This module uses python Mako template to generate Fortran code.
!! Make sure to run the `mako_preprocess.py` if you edit the source
!! .f90p file, and don`t edit the generated .f90 directly.
!!
!! In the call signatures below:
!! - `loc_id` is the HDF5 file handle,
!! - `dset_name` is the path to the dataset
!! - `error` is an optional argument with the error code, which should
!! typically *not* be passed unless you want to do the error checking
!! yourself and don`t want the code to die in case the HDF5 library
!! fails.
!!
!! Routines:
!! - hdf5_{read,write}_{int,double,complex,logical,scalar}
!! (loc_id, dset_name, buf, error):
!! integer(HID_T), intent(in) :: loc_id
!! character(len=*), intent(in) :: dset_name
!! {data_type}, intent({intent}), dimension(*) :: buf
!! integer, intent(out), optional :: error
!!
!! Reads/writes a single element store in buffer `buf`.
!! When writing, the dataset in `dset_name` may or may not exist, but if
!! it already exist, it must have the correct type and shape.
!!
!!
!! - hdf5_{read,write}_{int,double,complex,logical,scalar}_array
!! (loc_id, dset_name, dims, buf, error):
!! integer(HID_T), intent(in) :: loc_id
!! character(len=*), intent(in) :: dset_name
!! integer, intent(in) :: dims(:)
!! {data_type}, intent({intent}), dimension(*) :: buf
!! integer, intent(out), optional :: error
!!
!! Reads/writes an array of size `dims` from/to buffer `buf`.
!! This routine writes the full array, and not part thereof.
!! When writing, the dataset in `dset_name` may or may not exist, but if
!! it already exist, it must have the correct type and shape.
!!
!!
!! - hdf5_{read,write}_{int,double,complex,logical,scalar}_hyperslab
!! (loc_id, dset_name, rw_count, offset, buf, error)
!! integer(HID_T), intent(in) :: loc_id
!! character(len=*), intent(in) :: dset_name
!! integer, intent(in) :: rw_count(:)
!! integer, intent(in) :: offset(:)
!! {data_type}, intent({intent}), dimension(*) :: buf
!! integer, intent(out), optional :: error
!!
!! Reads/writes part of an array stored in `buf`, where;
!! * `rw_count` is the number of elements read/written from/to
!! the daset for each dimension. This does not need to be the same
!! as the size of the dataset, but it is the overall size of the
!! output/input array `buf`
!! * `offset` is the offset for each dimension, relative to the HDF5
!! dataset when reading/writing from/to the HDF5 file.
!!
!! Unlike the previous subroutines, the dataset must already exist
!! before calling it, since the size of the dataset in file may be
!! different than `rw_count`.
!!
!! Notes:
!! - Logical datasets are stored as integers.
!! * When writing, we map .false.=>0 and .true.=>1
!! * When reading, we map 0=>.false. and nonzero=>.true.
!!
!! - Because there is no native HDF5 support for complex data, we
!! read/write complex dataset as double dataset by appending one
!! extra fast dimension with size=2.
!!
!! - There are interfaces for scalar dataset, which either read/write
!! complex or double data depending on the type of `buf`. However,
!! to make the number of dimensions the same for double or complex
!! scalar data, the scalar routines will **also append an extra
!! dimension with size=1** when reading/writing double data.
!! The code will not add an extra dimension if you call the
!! non-scalar routine to write a double buffer.
!!
!! For example, suppose you want to write `buf(1:N)`:
!! * If `buf` is `real(DP)`:
!! - `write_double_array` writes a real dataset of dimension [N].
!! - `write_scalar_array` writes a real dataset of dimension [1,N].
!! * If `buf` is `complex(DPC)`:
!! - Both `write_complex_array` and `write_scalar_array` write
!! real datasets of dimension [2,N].
!!
!!=========================================================================
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
module hdf5_io_data_m
  use, intrinsic :: iso_c_binding
  use global_m
  use hdf5_io_safe_m
end module hdf5_io_data_m
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
