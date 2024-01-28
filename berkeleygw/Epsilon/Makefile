include ../make.inc

LIBEPSILON = ../Epsilon/libepsilon.a
F90_OPTS += -I ./ -I ../Common

SRC = \
algos_epsilon.f90  genwf_mpi.f90  lin_denominator.f90 epsinv.f90 \
chi_convergence.f90 input.f90      mtxel.f90            rqstar.f90 \
chi_summation.f90   input_q.f90    mtxelmultiply.f90 genwf_eps.f90     inread.f90     

C_SRC = 


MAIN_PROGRAMS = \
eps0sym.f90 setup_subsampling_nns.f90 \
epsascbin.f90 epsbinasc.f90 ploteps.f90 \
epsinvomega.f90 epsomega.f90 \
epsmat_merge.f90 epsilon_main.f90 printchi.f90


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
	ar rcs $(LIBEPSILON) *.o

clean:
	rm -rf *.o *.mod $(LIBEPSILON) *.x


include make.depends