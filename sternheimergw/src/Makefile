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
-I$(QEPATH)/LR_Modules

LIBS = \
libsgwmain.a \
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

SRC = \
aaa.f90                   green.f90           print_matel_im.f90 \
allocate_gwq.f90          gwcom.f90           q_points.f90 \
analytic.f90                                  reorder.f90 \
array.f90                 gw_input.f90        run_nscf.f90 \
assert.f90                gw_opening.f90      save_gw_input.f90 \
bicgstab.f90              gwq_init.f90        select_solver.f90 \
check_initial_status.f90  gwq_readin.f90      setup_nscf.f90 \
clean_pw_gw.f90           gwq_setup.f90       setup_nscf_green.f90 \
close_gwq.f90             gw_type.f90         sigma_expect.f90 \
constant.f90              initialize_gw.f90   sigma.f90 \
coulomb.f90               invert_epsilon.f90  sigma_grid.f90 \
coulomb_q0G0.f90          kind_gw.f90         sigma_io.f90 \
coulomb_vcut.f90          linear_op.f90       sigma_matel.f90 \
coulpade.f90              linear_solver.f90   sigma_pade.f90 \
data_lapack.f90           mix_pot_c.f90       sleep.f90 \
deallocate_gwq.f90        norm.f90            solve_linter.f90 \
debug.f90                 openfilq.f90        sort.f90 \
do_stern.f90              opengwfil.f90       star_q.f90 \
dvqpsi_us.f90             pade.f90            stern_symm.f90 \
exchange.f90              pade_robust.f90     stop_gw.f90 \
fft6.f90                  parallel.f90        timing.f90 \
freqbins.f90              plot_coulomb.f90    truncation.f90 \
gauleg_grid.f90           pp_output.f90       unfold_w.f90 \
gmap_sym.f90              prepare_q0.f90      vendor_lapack.f90 \
godby_needs.f90           prepare_q.f90       version.f90 \
gram_schmidt.f90          print_matel.f90

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

all: lib main

# Targets
lib: $(OBJ)
	ar rcs libsgwmain.a *.o

main: $(OBJ)
	$(F90) $(F90_OPTS) $(INCLUDES) -c gw.f90 -o gw.o
	$(F90) -g gw.o -o gw.x $(LIBS)

clean:
	rm -rf *.o *.mod libsgwmain.a *.x

include make.depends