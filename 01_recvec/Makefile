# Makefile for simple test programs

include ../qe/make.inc
include ../make.inc.ffr

SRC = \
setup_structure_Si8.f90 \
setup_structure_H.f90 \
setup_fft.f90 \
setup_symmetry.f90 \
setup_kpoints_noshift.f90 \
setup_kpoints_000.f90 \
setup_gvect.f90

OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(SRC:.c=.o)

# Suffix rules
.SUFFIXES: .o .f90

.f90.o:
	$(F90) $(F90_OPTS) $(INC) -c $<

%.mod :
	@if [! -f $@ ]; then rm $(*F).o; fi
	$(MAKE) $<

#
# Here are the targets
#

lib: $(OBJ)
	ar rcs libmain.a *.o

main_2:
	$(F90) $(F90_OPTS) $(INC) main_2.f90 -o main_2.x libmain.a $(LIBS)

main_1:
	$(F90) $(F90_OPTS) $(INC) main_1.f90 -o main_1.x libmain.a $(LIBS)

pwt_gvect:
	$(F90) $(F90_OPTS) $(INC) pwt_gvect.f90 -o pwt_gvect.x ../common/libcommon.a $(LIBS)

pwt_ewald:
	$(F90) $(F90_OPTS) $(INC) pwt_ewald.f90 -o pwt_ewald.x $(LIBS)

clean:
	rm -f *.o *.mod libmain.a x-main *.x


