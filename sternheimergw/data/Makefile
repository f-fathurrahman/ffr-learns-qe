#include make.inc

QEPATH = /home/efefer/WORKS/STERNHEIMERGW/q-e-qe-6.3

F90 = mpif90
LD = mpif90
F90_OPTS = -c -O3 -g -x f95-cpp-input -D__FFTW3 -D__MPI
LD_OPTS = -g
INCLUDES = \
-I$(QEPATH)/include \
-I$(QEPATH)/FoX/finclude \
-I$(QEPATH)/S3DE/iotk/include \
-I$(QEPATH)/iotk/src \
-I$(QEPATH)/Modules \
-I$(QEPATH)/FFTXlib \
-I$(QEPATH)/LAXlib \
-I$(QEPATH)/UtilXlib \
-I$(QEPATH)/FoX/finclude \
-I$(QEPATH)/PW/src

LIBS = \
$(QEPATH)/PW/src/libpw.a \
$(QEPATH)/GWW/gww/libgww.a \
$(QEPATH)/Modules/libqemod.a \
$(QEPATH)/FFTXlib/libqefft.a \
$(QEPATH)/KS_Solvers/libks_solvers.a \
$(QEPATH)/LAXlib/libqela.a \
$(QEPATH)/UtilXlib/libutil.a \
$(QEPATH)/clib/clib.a \
$(QEPATH)/iotk/src/libiotk.a \
-L$(QEPATH)/FoX/lib  -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys \
-lfftw3 \
-L/home/efefer/mysoftwares/LIBMKL -lmkl_gf_lp64  -lmkl_sequential -lmkl_core


SRC = fft6.f90  gram_schmidt.f90  lapack.f90 norm.f90  parallel.f90  sleep.f90  timing.f90

OBJ = $(SRC:.f90=.o) $(SRC:.f=.o)

#
# Suffix rule for Fortran 90
#
#%.mod :
#	@if [! -f $@ ]; then \
#		rm $(*F).o; \
#		fi
#	$(MAKE) $<

%.o : %.f90
	$(F90) $(F90_OPTS) $(INCLUDES) -c -o $(*F).o $<

#
# Fortran 77 sources
# supress warning
#
.SUFFIXES: .o .f
.f.o:
	$(F77) $(F77_OPTS) -c $<

all: lib

# Targets
lib: $(OBJ)
	ar rcs libsgwdata.a *.o

clean:
	rm -rf *.o *.mod libsgwdata.a *.x

fft6.o : fft6.f90 timing.o 
gram_schmidt.o : gram_schmidt.f90 norm.o 
lapack.o : lapack.f90 
norm.o : norm.f90 
parallel.o : parallel.f90 
sleep.o : sleep.f90 
timing.o : timing.f90