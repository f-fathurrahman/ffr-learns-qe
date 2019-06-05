# Makefile for simple test programs

include qe/make.inc
include make.inc.ffr

#
# Here are the targets
#

my_scf:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 my_scf.f90 -o my_scf.x $(LIBS_ALL)

pwt_plot_vltot:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 pwt_plot_vltot.f90 -o pwt_plot_vltot.x $(LIBS_ALL)

01_new_gvectors:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 01_new_gvectors.f90 -o 01_new_gvectors.x $(LIBS_ALL)

02_new_vlocal:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 02_new_vlocal.f90 -o 02_new_vlocal.x $(LIBS_ALL)

03_test_dfft:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 03_test_dfft.f90 -o 03_test_dfft.x $(LIBS_ALL)

04_test_v_of_0:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 04_test_v_of_0.f90 -o 04_test_v_of_0.x $(LIBS_ALL)

05_test_symmetry:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 05_test_symmetry.f90 -o 05_test_symmetry.x $(LIBS_ALL)

test_pw_01:
	$(MPIF90) -Wall -g $(INC_ALL) test_pw_01.f90 -o test_pw_01.x $(LIBS_ALL)

pwscf:
	$(MPIF90) -Wall -g $(INC_ALL) pwscf.f90 -o pwscf.x $(LIBS_ALL)


clean:
	rm -f *.o *.mod libmain.a x-main *.x

