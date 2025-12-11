# make.inc.  Generated from make.inc.in by configure.

# compilation rules

.SUFFIXES :
.SUFFIXES : .o .c .f90 .h .fh

# most fortran compilers can directly preprocess c-like directives: use
# 	$(MPIF90) $(F90FLAGS) -c $<
# if explicit preprocessing by the C preprocessor is needed, use:
#	$(CPP) $(CPPFLAGS) $< -o $*.F90
#	$(MPIF90) $(F90FLAGS) -c $*.F90 -o $*.o
# remember the tabulator in the first column !!!

.f90.o:
	$(MPIF90) $(F90FLAGS) -c $< -o $(@)

.c.o:
	$(CC) $(CFLAGS)  -c $<

.h.fh:
	$(CPP) $(CPPFLAGS) $< -o $*.fh



# Top QE directory, useful for locating libraries,  linking QE with plugins
# The following syntax should always point to TOPDIR:
TOPDIR = $(dir $(abspath $(filter %make.inc,$(MAKEFILE_LIST))))
# if it doesn't work, uncomment the following line (edit if needed):

# TOPDIR = /mnt/d/WORKS/my_github_repos/ffr-learns-qe/qe-7.5

# DFLAGS  = precompilation options (possible arguments to -D and -U)
#           used by the C compiler and preprocessor
# See include/defs.h.README for a list of options and their meaning
# To use libxc (v>=4.3.0), add -D__LIBXC to DFLAGS
# FDFLAGS = precompilation options used by the Fortran compiler
# With the exception of IBM xlf, FDFLAGS = $(DFLAGS)
# For IBM xlf, FDFLAGS is the same as DFLAGS with separating commas

# MANUAL_DFLAGS  = additional Fortran precompilation option(s), if desired
#                  BEWARE: does not work for IBM xlf, manually edit FDFLAGS
MANUAL_DFLAGS  =
DFLAGS         =  -D__FFTW3 -D__MPI -D__MPI_MODULE
FDFLAGS        = $(DFLAGS) $(MANUAL_DFLAGS)

# IFLAGS = how to locate directories with *.h or *.f90 file to be included
#          typically -I$(TOPDIR)/include -I/some/other/directory/
#          the latter contains .e.g. files needed by FFT libraries
# for libxc add -I/path/to/libxc/include/

IFLAGS         = -I. -I$(TOPDIR)/include

# MOD_FLAG = flag used by f90 compiler to locate modules

MOD_FLAG      = -I

# BASEMOD_FLAGS points to directories containing basic modules,
# while BASEMODS points to the corresponding module libraries.
# More package-specific directories can be added in each Makefile
# and stored into MODFLAGS, and the same for module libraries in QEMODS

BASEMOD_FLAGS= $(MOD_FLAG)$(TOPDIR)/upflib \
               $(MOD_FLAG)$(TOPDIR)/XClib \
               $(MOD_FLAG)$(TOPDIR)/Modules \
               $(MOD_FLAG)$(TOPDIR)/FFTXlib_src \
	       $(MOD_FLAG)$(TOPDIR)/LAXlib \
	       $(MOD_FLAG)$(TOPDIR)/UtilXlib \
	       $(MOD_FLAG)$(TOPDIR)/MBD \
	       $(MOD_FLAG)$(TOPDIR)/KS_Solvers $(FOX_MOD)

# If A depends upon B, A should come before B in the list below
# (most compilers don't care but some don't resolve cross links)
BASEMODS= $(TOPDIR)/Modules/libqemod.a \
          $(TOPDIR)/upflib/libupf.a \
          $(TOPDIR)/XClib/xc_lib.a \
          $(TOPDIR)/FFTXlib_src/libqefft.a \
          $(TOPDIR)/LAXlib/libqela.a \
          $(TOPDIR)/UtilXlib/libutil.a \
          $(TOPDIR)/MBD/libmbd.a

# Compilers: fortran-90, fortran-77, C
# If a parallel compilation is desired, MPIF90 should be a fortran-90
# compiler that produces executables for parallel execution using MPI
# (such as for instance mpif90, mpf90, mpxlf90,...);
# otherwise, an ordinary fortran-90 compiler (f90, g95, xlf90, ifort,...)

MPIF90         = mpif90
F90           = gfortran
CC             = cc

# GPU architecture (Kepler: 35, Pascal: 60, Volta: 70, Turing: 75, Ampere: 80)
GPU_ARCH=

# CUDA runtime (should be compatible with the CUDA driver)
CUDA_RUNTIME=

# CUDA F90 Flags
CUDA_F90FLAGS= $(MOD_FLAG)$(TOPDIR)/external/devxlib/src

# CUDA C Flags
CUDA_CFLAGS=

# C preprocessor and preprocessing flags - for explicit preprocessing,
# if needed (see the compilation rules above)
# preprocessing flags must include DFLAGS and IFLAGS

CPP            = cpp
CPPFLAGS       = -P -traditional -Uvector $(DFLAGS) $(IFLAGS)

# compiler flags: C, F90
# C flags must include DFLAGS and IFLAGS
# F90 flags must include MODFLAGS, IFLAGS, and FDFLAGS with appropriate syntax

CFLAGS         = -O3 $(DFLAGS) $(IFLAGS) $(CUDA_CFLAGS)
F90FLAGS       = $(FFLAGS) -cpp $(FDFLAGS) $(CUDA_F90FLAGS) $(IFLAGS) $(MODFLAGS)

# compiler flags with and without optimization for fortran-77
# the latter is NEEDED to properly compile dlamch.f, used by lapack
FFLAGS         = -O3 -g -fallow-argument-mismatch
FFLAGS_NOOPT   = -O0 -g

# compiler flag needed by some compilers when the main program is not fortran
# Currently used for Yambo

FFLAGS_NOMAIN   = 

# Linker, linker-specific flags (if any)
# Typically LD coincides with F90 or MPIF90, LD_LIBS is empty
# for libxc, set LD_LIBS=-L/path/to/libxc/lib/ -lxcf03 -lxc
# If libxc release is 5.0.0 replace -lxcf03 with -lxcf90

LD             = mpif90
LDFLAGS        = -g
LD_LIBS        = 

# External Libraries (if any) : blas, lapack, fft, MPI

# If you have nothing better, use the local copy
# BLAS_LIBS = $(TOPDIR)/external/lapack/libblas.a

BLAS_LIBS      =  -lblas 

# If you have nothing better, use the local copy
# LAPACK = liblapack
# LAPACK_LIBS = $(TOPDIR)/external/lapack/liblapack.a

LAPACK      = 
LAPACK_LIBS = -L/usr/local/lib -llapack  -lblas 

SCALAPACK_LIBS = 

# nothing needed here if the the internal copy of FFTW is compiled
# (needs -D__FFTW in DFLAGS)

FFT_LIBS       =  -lfftw3 

# HDF5
HDF5_LIBS = 

# FOX
FOX      = 
FOX_MOD  = 
FOX_LIB  = 
FOX_FLAGS = 

# ENVIRON
ENVIRON_LIBS = 

# MPI libraries (should not be needed)

MPI_LIBS       = 

# IBM-specific: MASS libraries, if available and if -D__MASS is defined in FDFLAGS

MASS_LIBS      = 

# CUDA libraries
CUDA_LIBS= -L$(TOPDIR)/external/devxlib/src -ldevXlib
CUDA_EXTLIBS = devxlib

# ar command and flags - for most architectures: AR = ar, ARFLAGS = ruv

AR             = ar
ARFLAGS        = ruv

# ranlib command. If ranlib is not needed (it isn't in most cases) use
# RANLIB = echo

RANLIB         = ranlib

# all internal and external libraries - do not modify

FLIB_TARGETS   = all

LIBXC_LIBS     = 
QELIBS         = $(LIBXC_LIBS) \
                 $(CUDA_LIBS) $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(FOX_LIB) \
		 		 $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(MASS_LIBS) $(HDF5_LIBS) \
			 	 $(ENVIRON_LIBS) $(LD_LIBS)

# wget or curl - useful to download from network
WGET = wget -O

# Install directory - "make install" copies *.x executables there
PREFIX = /usr/local
