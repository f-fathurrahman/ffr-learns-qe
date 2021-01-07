# Makefile for simple test programs

include ../qe-6.6/make.inc
include make.inc.ffr

#
# Here are the targets
#

minipw:
	$(MPIF90) -Wall -g $(INC_ALL) minipw.f90 -o minipw.x $(LIBS_ALL)

clean:
	rm -f *.o *.mod libmain.a minipw.x

