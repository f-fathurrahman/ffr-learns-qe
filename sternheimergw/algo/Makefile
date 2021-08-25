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
-I../util

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

SRC = allocate_gwq.f90  gauleg_grid.f90    opengwfil.f90      setup_nscf.f90 \
analytic.f90              gmap_sym.f90       pade.f90           setup_nscf_green.f90 \
bicgstab.f90              godby_needs.f90    pade_robust.f90    sigma_grid.f90 \
check_initial_status.f90  gwq_init.f90       pp_output.f90      sigma_io.f90 \
clean_pw_gw.f90           gwq_setup.f90      q_points.f90       star_q.f90 \
close_gwq.f90             initialize_gw.f90  reorder.f90        stern_symm.f90 \
coulomb_vcut.f90          linear_op.f90      run_nscf.f90       stop_gw.f90 \
deallocate_gwq.f90        linear_solver.f90  save_gw_input.f90  truncation.f90 \
freqbins.f90              openfilq.f90       select_solver.f90  unfold_w.f90

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
	ar rcs libsgwalgo.a *.o

clean:
	rm -rf *.o *.mod libsgwalgo.a *.x

allocate_gwq.o : allocate_gwq.f90 
analytic.o : analytic.f90 sigma_grid.o pade_robust.o godby_needs.o freqbins.o 
bicgstab.o : bicgstab.f90 
check_initial_status.o : check_initial_status.f90 
clean_pw_gw.o : clean_pw_gw.f90 
close_gwq.o : close_gwq.f90 
coulomb_vcut.o : coulomb_vcut.f90 
deallocate_gwq.o : deallocate_gwq.f90 
freqbins.o : freqbins.f90 
gauleg_grid.o : gauleg_grid.f90 
gmap_sym.o : gmap_sym.f90 
godby_needs.o : godby_needs.f90 
gwq_init.o : gwq_init.f90 
gwq_setup.o : gwq_setup.f90 
initialize_gw.o : initialize_gw.f90 
linear_op.o : linear_op.f90 
linear_solver.o : linear_solver.f90 
openfilq.o : openfilq.f90 save_gw_input.o 
opengwfil.o : opengwfil.f90 sigma_io.o sigma_grid.o 
pade.o : pade.f90 
pade_robust.o : pade_robust.f90 
pp_output.o : pp_output.f90 
q_points.o : q_points.f90 
reorder.o : reorder.f90 
run_nscf.o : run_nscf.f90 setup_nscf_green.o 
save_gw_input.o : save_gw_input.f90 
select_solver.o : select_solver.f90 linear_solver.o bicgstab.o 
setup_nscf.o : setup_nscf.f90 
setup_nscf_green.o : setup_nscf_green.f90 
sigma_grid.o : sigma_grid.f90 freqbins.o 
sigma_io.o : sigma_io.f90 
star_q.o : star_q.f90 
stern_symm.o : stern_symm.f90 
stop_gw.o : stop_gw.f90 
truncation.o : truncation.f90 coulomb_vcut.o 
unfold_w.o : unfold_w.f90
