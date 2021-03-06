# Makefile for simple test programs

include qe/make.inc
include make.inc.ffr

#
# Here are the targets
#
test_radial_grid:
	$(MPIF90) -Wall -g $(INC_ALL) test_radial_grid.f90 -o test_radial_grid.x $(LIBS_ALL)

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

06_test_my_symm_base:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 06_test_my_symm_base.f90 -o 06_test_my_symm_base.x $(LIBS_ALL)

07_test_kpoint_grid:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 07_test_kpoint_grid.f90 -o 07_test_kpoint_grid.x $(LIBS_ALL)

08_my_ewald:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 08_my_ewald.f90 -o 08_my_ewald.x $(LIBS_ALL)

09_my_rgen:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 09_my_rgen.f90 -o 09_my_rgen.x $(LIBS_ALL)

10_md_parameters:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 10_md_parameters.f90 -o 10_md_parameters.x $(LIBS_ALL)

11_paw_variables:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 11_paw_variables.f90 -o 11_paw_variables.x $(LIBS_ALL)

12_upf_uspp:
	$(MPIF90) -Wall -g $(INC_ALL) prepare_all.f90 12_upf_uspp.f90 -o 12_upf_uspp.x $(LIBS_ALL)

test_pw_01:
	$(MPIF90) -Wall -g $(INC_ALL) test_pw_01.f90 -o test_pw_01.x $(LIBS_ALL)

pwscf:
	$(MPIF90) -Wall -g $(INC_ALL) pwscf.f90 -o pwscf.x $(LIBS_ALL)


clean:
	rm -f *.o *.mod libmain.a x-main *.x

