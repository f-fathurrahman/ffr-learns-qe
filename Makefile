# Makefile for simple test programs

include SYS/make.inc

#SRC = 
#OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(SRC:.c=.o)

# Suffix rules
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) $(INC_ALL) -c $<

%.mod :
	@if [! -f $@ ]; then \
    rm $(*F).o;  \
  fi
	$(MAKE) $<

#
# Here are the targets
#

lib: $(OBJ)
	ar rcs libmain.a *.o

test_upf_nc:
	$(F90) $(F90_OPTS) $(INC_ALL) test_upf_nc.f90 $(LIB_ALL) -o test_upf_nc.x

test_gth:
	$(F90) $(F90_OPTS) $(INC_ALL) test_gth.f90 $(LIB_ALL) -o test_gth.x

test_uspp:
	$(F90) $(F90_OPTS) $(INC_ALL) test_uspp.f90 $(LIB_ALL) -o test_uspp.x

pw01:
	$(F90) $(F90_OPTS) $(INC_ALL) pw01.f90 $(LIB_ALL) -o pw01.x

pw02:
	$(F90) $(F90_OPTS) $(INC_ALL) pw02.f90 $(LIB_ALL) -o pw02.x

t_mp:
	$(F90) $(F90_OPTS) $(INC_ALL) t_mp.f90 $(LIB_ALL) -o t_mp.x

pwt_electrons_v2:
	$(F90) $(F90_OPTS) $(INC_ALL) cheFSI.f90 pwt_electrons_v2.f90 \
		$(LIB_QE_MATH) -o pwt_electrons_v2.x

pwt_electrons_v1:
	$(F90) $(F90_OPTS) $(INC_ALL) pwt_electrons_v1.f90 $(LIB_QE_MATH) \
		-o pwt_electrons_v1.x

pwt_plot_vltot:
	$(F90) $(F90_OPTS) $(INC_ALL) pwt_plot_vltot.f90 $(LIB_QE_MATH) \
		-o pwt_plot_vltot.x

t_pseudo_v1:
	$(F90) $(F90_OPTS) $(INC_ALL) t_pseudo_v1.f90 $(LIB_QE_MATH) \
		-o t_pseudo_v1.x

free_fcc:
	$(F90) $(F90_OPTS) $(INC_ALL) free_fcc.f90 $(LIB_QE_MATH) \
		-o free_fcc.x

t_gvect_v2:
	$(F90) $(F90_OPTS) $(INC_ALL) t_gvect_v2.f90 $(LIB_QE_MATH) \
	libmain.a	-o t_gvect_v2.x

coul_fcc:
	$(F90) $(F90_OPTS) $(INC_ALL) coul_fcc.f90 $(LIB_QE_MATH) \
		-o coul_fcc.x

t_gvect_v1:
	$(F90) $(F90_OPTS) $(INC_ALL) t_gvect_v1.f90 $(LIB_QE_MATH) \
		-o t_gvect_v1.x

t_struct_v1:
	$(F90) $(F90_OPTS) $(INC_ALL) t_struct_v1.f90 $(LIB_QE_MATH) \
		-o t_struct_v1.x

get_evc:
	$(F90) $(F90_OPTS) $(INC_ALL) get_evc.f90 $(LIB_QE_MATH) -o get_evc.x

pwt_uspp_param_v1:
	$(F90) $(F90_OPTS) $(INC_ALL) pwt_uspp_param_v1.f90 $(LIB_QE_MATH) -o pwt_uspp_param_v1.x

t_uspp_param_r:
	$(F90) $(F90_OPTS) $(INC_ALL) t_uspp_param_r.f90 $(LIB_QE_MATH) -o t_uspp_param_r.x

clean:
	rm -f *.o *.mod libmain.a x-main *.x

