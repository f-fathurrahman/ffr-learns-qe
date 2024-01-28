include ../make.inc

LIBSIGMA = ../Sigma/libsigma.a
F90_OPTS += -I ./ -I ../Common

SRC = \
algos_sigma.f90  genwf_mpi.f90    inread.f90     mtxel_occ.f90  shiftenergy.f90    write_result_dyn.f90 \
ch_converge.f90  input.f90        mtxel.f90      mtxel_vxc.f90  shiftenergy_dyn.f90  wpeff.f90 \
write_result_dyn_hp.f90 \
epscopy.f90      input_outer.f90  mtxel_cor.f90  write_result.f90  write_result_hp.f90

C_SRC = 


MAIN_PROGRAMS = sigma_main.f90 offdiag.f90    sig2wan.f90


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
	ar rcs $(LIBSIGMA) *.o

clean:
	rm -rf *.o *.mod $(LIBSIGMA) *.x


include make.depends