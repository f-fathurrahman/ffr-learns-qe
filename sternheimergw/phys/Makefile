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
-I$(QEPATH)/PW/src \
-I$(QEPATH)/LR_Modules \
-I../vendor \
-I../data \
-I../util \
-I../algo

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

SRC = coulomb.f90       exchange.f90        prepare_q0.f90      sigma.f90 \
coulomb_q0G0.f90  green.f90           prepare_q.f90       sigma_matel.f90 \
coulpade.f90      invert_epsilon.f90  print_matel.f90     sigma_pade.f90 \
do_stern.f90      mix_pot_c.f90       print_matel_im.f90  solve_linter.f90 \
dvqpsi_us.f90     plot_coulomb.f90    sigma_expect.f90

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
	ar rcs libsgwphys.a *.o

clean:
	rm -rf *.o *.mod libsgwphys.a *.x

coulomb.o : coulomb.f90 solve_linter.o 
coulomb_q0G0.o : coulomb_q0G0.f90 solve_linter.o 
coulpade.o : coulpade.f90 
do_stern.o : do_stern.f90 plot_coulomb.o 
dvqpsi_us.o : dvqpsi_us.f90 
exchange.o : exchange.f90 
green.o : green.f90 
invert_epsilon.o : invert_epsilon.f90 
mix_pot_c.o : mix_pot_c.f90 
plot_coulomb.o : plot_coulomb.f90 
prepare_q0.o : prepare_q0.f90 
prepare_q.o : prepare_q.f90 
print_matel.o : print_matel.f90 
print_matel_im.o : print_matel_im.f90 
sigma_expect.o : sigma_expect.f90 
sigma.o : sigma.f90 green.o coulpade.o 
sigma_matel.o : sigma_matel.f90 sigma_expect.o 
sigma_pade.o : sigma_pade.f90 
solve_linter.o : solve_linter.f90
