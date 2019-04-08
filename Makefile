# Makefile for simple test programs

include qe/make.inc
include make.inc.ffr


#
# Here are the targets
#

01_new_gvectors:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 01_new_gvectors.f90 -o 01_new_gvectors.x $(LIBS_ALL)

02_new_vlocal:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 02_new_vlocal.f90 -o 02_new_vlocal.x $(LIBS_ALL)


test_pw_01:
	$(MPIF90) -Wall -g $(INC_ALL) test_pw_01.f90 -o test_pw_01.x $(LIBS_ALL)

pwscf:
	$(MPIF90) -Wall -g $(INC_ALL) pwscf.f90 -o pwscf.x $(LIBS_ALL)

clean:
	rm -f *.o *.mod libmain.a x-main *.x

