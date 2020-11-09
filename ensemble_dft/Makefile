# Makefile for simple test programs

include ../qe-6.5/make.inc
include make.inc.ffr

#
# Here are the targets
#

my_cprstart:
	$(MPIF90) -Wall -g $(INC_ALL) my_cprstart.f90 -o my_cprstart.x $(LIBS_ALL)

clean:
	rm -f *.o *.mod libmain.a x-main *.x

