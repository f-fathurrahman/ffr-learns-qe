include ../make.inc

LIBCOMMON = ../Common/libcommon.a
F90_OPTS += -I ./

SRC = \
acc_linalg.f90           evecs.f90           irrbz.f90          random.f90             symmetries.f90 \
acc_mtxel_kernels.f90    fft_parallel.f90    kernel_io.f90      read_cube.f90           \
algos_common.f90         fftw.f90            lapack.f90         read_matrix.f90        tile.f90 \
bessel.f90               find_kpt_match.f90  message.f90        read_rho_vxc.f90       timing.f90 \
bgw_mpi.f90              fixwings.f90        minibzaverage.f90  references.f90         trunc_cell_box.f90 \
blas.f90                 fullbz.f90          misc.f90           scalapack.f90          trunc_cell_box_d.f90 \
cells.f90                genwf.f90           mtxel_optical.f90  scalapack_aux.f90      trunc_cell_wire.f90 \
check_inversion.f90      global.f90          norm.f90           trunc_scell_box_d.f90 \
check_screening.f90      gmap.f90            nrtype.f90         scissors.f90           typedefs.f90 \
checkbz.f90              groupk.f90          sort.f90           vcoul_generator.f90 \
checkgriduniformity.f90  input_utils.f90     os.f90             spinor_symmetries.f90  version.f90 \
createpools.f90          inread_common.f90   peinfo.f90         splines.f90            wfn_rho_vxc_io.f90 \
eqpcor.f90               inversion.f90       subgrp.f90         write_matrix.f90 \
io_utils.f90             push_pop.f90        susymmetries.f90   write_program_header.f90

C_SRC = os_c.c  periodic_table.c


MAIN_PROGRAMS = \
print_version_info.f90 \
test_evecs.f90


OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(C_SRC:.c=.o)

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

%.o : %.c
	$(CC) $(CC_OPTS) -c -o $(*F).o $<

#
# Fortran 77 sources
# supress warning
.SUFFIXES: .o .f
.f.o:
	$(F77) $(F77_OPTS) -c $<


# Targets
lib: $(OBJ)
	ar rcs $(LIBCOMMON) *.o

clean:
	rm -rf *.o *.mod $(LIBCOMMON) *.x


include make.depends