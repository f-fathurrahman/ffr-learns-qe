MAKE = make
AR = ar

LIBXC_HOME = /home/efefer/mysoftwares/libxc-4.3.4dyn

DFLAGS = -D__FFTW -D__LIBXC

F77 = gfortran
F77_OPTS = -Wall -cpp -fPIC

F90 = gfortran
F90_OPTS = -I$(LIBXC_HOME)/include $(DFLAGS) -Wall -cpp -fPIC

include OBJECTS_QE.make

include OBJECTS_FFR.make

OBJ = $(OBJ_FFTXLIB) \
$(OBJ_LAXLIB) \
$(OBJ_MODULES) \
$(OBJ_DFTD3) \
$(OBJ_PW) \
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
	gcc -Wall -fPIC -c -o $(*F).o $<

#
# Targets
#
libdyn:	$(OBJ)
	gfortran -shared -o libqemain.so $(OBJ)

clean:
	rm -rf *.o *.mod libqemain.so *.x

include DEPENDS.make