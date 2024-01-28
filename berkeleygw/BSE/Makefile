include ../make.inc

LIBBSE = ../BSE/libbse.a
F90_OPTS += -I ./ -I ../Common

SRC = \
absh.f90 epsmodel.f90  input_q.f90  mtxel_jdos.f90 \
g_sum.f90 inread.f90 mtxel_kernel.f90 \
absp.f90 genwf_co.f90 inread_kernel.f90 \
absp0.f90  bsewrite.f90 genwf_kernel.f90 mtxel_t.f90 \
absp_io.f90 diag.f90            gx_sum.f90 interp.f90  \
absp_lanczos.f90     diagonalize.f90     haydock.f90  intkernel.f90 sortbyq.f90 \
algos_kernel.f90     distrib.f90         input_co.f90 intpts.f90 \
bse_convert.f90      distrib_kernel.f90  input_co_q.f90 intwfn.f90  vmtxel.f90 \
epscopy.f90         input_fi.f90   iterate.f90 \
bse_init.f90         epsdiag.f90         input_kernel.f90  w_sum.f90

C_SRC = 


MAIN_PROGRAMS = absorption.f90 bseascbin.f90 bsebinasc.f90 bsemat_check.f90 \
inteqp.f90 kernel.f90 setup_subsampling_csi.f90 summarize_eigenvectors.f90 \
vmtxel_main.f90


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
	ar rcs $(LIBBSE) *.o

clean:
	rm -rf *.o *.mod $(LIBBSE) *.x


include make.depends