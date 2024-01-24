MAKE = make
AR = ar

LIBXC_HOME = /home/efefer/mysoftwares/libxc-4.3.4dyn

DFLAGS = -D__FFTW -D__LIBXC

F77 = gfortran
F77_OPTS = -O3 -Wall -cpp -fPIC

F90 = gfortran
F90_OPTS = -O3 -I$(LIBXC_HOME)/include $(DFLAGS) -Wall -cpp -fPIC

include OBJECTS_QE.make

include OBJECTS_FFR.make

include OBJECTS_LD1X.make

OBJ = $(OBJ_FFTXLIB) \
$(OBJ_LAXLIB) \
$(OBJ_MODULES) \
$(OBJ_DFTD3) \
$(OBJ_PW) \
$(OBJ_LD1X) \
$(OBJ_FFR) \
$(SRC:.f90=.o) $(SRC:.f=.o)

#
# Suffix rule for Fortran 90
#
%.mod :
	@if [! -f $@ ]; then \
		rm $(*F).o; \
		fi
	$(MAKE) $<

%.o : %.f90
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

#
# Fortran 77 sources
# supress warning
.SUFFIXES: .o .f
.f.o:
	$(F77) $(F77_OPTS) -c $<


%.o : %.c
	gcc -O3 -Wall -fPIC -c -o $(*F).o $<

#
# Targets
#
libdyn:	$(OBJ)
	gfortran -shared -o libqemain.so $(OBJ) -lblas -llapack -Wl,-rpath $(LIBXC_HOME)/lib -L$(LIBXC_HOME)/lib -lxcf03 -lxc

clean:
	rm -rf *.o *.mod libqemain.so *.x

include DEPENDS.make
