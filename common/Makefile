# Makefile for simple test programs

include ../SYS/make.inc

SRC = \
t_import_gvect.f90 \
t_import_gvecs.f90 \
t_import_gvecw.f90 \
info_gvect.f90 \
info_gvecs.f90 \
info_structure.f90

OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(SRC:.c=.o)

INC  = -I$(INC_QEMOD) -I$(INC_PW) -I$(INC_QEFFT)

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
	ar rcs libcommon.a *.o

clean:
	rm -f *.o *.mod libcommon.a x-main *.x


