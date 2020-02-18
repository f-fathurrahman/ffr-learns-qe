# Makefile for simple test programs

include ../SYS/make.inc

SRC =

OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(SRC:.c=.o)

QEMODS=../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a
INC  = -I$(INC_QEMOD) -I$(INC_PW) -I$(INC_QEFFT)
LIBS = $(LIB_PW) $(LIB_QEMOD) $(LIB_QEFFT) $(LIB_QELA) $(LIB_CLIB) $(LIB_IOTK) $(LAPACK_LIBS) $(FFT_LIBS) $(BLAS_LIBS)


# Suffix rules
.SUFFIXES: .o .f90

.f90.o:
	$(F90) $(F90_OPTS) $(INC) -c $<

%.mod :
	@if [! -f $@ ]; then \
    rm $(*F).o;  \
  fi
	$(MAKE) $<

#
# Here are the targets
#

test_uspp:
	$(F90) $(F90_OPTS) $(INC) test_uspp.f90 -o test_uspp.x $(LIBS)

t_uspp_param_r:
	$(F90) $(F90_OPTS) $(INC) t_uspp_param_r.f90 -o t_uspp_param_r.x $(LIBS)

clean:
	rm -f *.o *.mod libmain.a *.x


