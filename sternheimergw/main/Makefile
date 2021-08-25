#include make.inc

QEPATH = /home/efefer/WORKS/STERNHEIMERGW/q-e-qe-6.3

F90 = mpif90
LD = mpif90
F90_OPTS = -O3 -g -x f95-cpp-input -D__FFTW3 -D__MPI
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
-I$(QEPATH)/PW/src \
-I$(QEPATH)/LR_Modules \
-I../vendor \
-I../data \
-I../util \
-I../algo \
-I../phys

LIBS = \
libsgwmain.a \
../phys/libsgwphys.a \
../algo/libsgwalgo.a \
../data/libsgwdata.a \
../util/libsgwutil.a \
../vendor/libsgwvendor.a \
$(QEPATH)/LR_Modules/liblrmod.a \
$(QEPATH)/PW/src/libpw.a \
$(QEPATH)/Modules/libqemod.a \
$(QEPATH)/FFTXlib/libqefft.a \
$(QEPATH)/KS_Solvers/Davidson/libdavid.a \
$(QEPATH)/KS_Solvers/CG/libcg.a \
$(QEPATH)/LAXlib/libqela.a \
$(QEPATH)/UtilXlib/libutil.a \
$(QEPATH)/dft-d3/libdftd3qe.a \
$(QEPATH)/clib/clib.a \
$(QEPATH)/iotk/src/libiotk.a \
-L$(QEPATH)/FoX/lib  -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys \
-lfftw3 \
-L/home/efefer/mysoftwares/LIBMKL -lmkl_gf_lp64  -lmkl_sequential -lmkl_core

SRC = gw_input.f90  gw_opening.f90  gwq_readin.f90  version.f90

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
	ar rcs libsgwmain.a *.o

main: $(OBJ)
	$(F90) $(F90_OPTS) $(INCLUDES) -c gw.f90 -o gw.o
	$(F90) -g gw.o -o gw.x $(LIBS)

clean:
	rm -rf *.o *.mod libsgwmain.a *.x

gw.o : gw.f90 gw_opening.o 
gw_input.o : gw_input.f90 
gw_opening.o : gw_opening.f90 version.o 
gwq_readin.o : gwq_readin.f90 gw_input.o 
version.o : version.f90