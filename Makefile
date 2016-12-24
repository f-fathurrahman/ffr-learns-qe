# Makefile for simple test programs

include SYS/make.inc.centosffr

SRC = t_import_gvect.f90 t_import_gvecs.f90 t_import_gvecw.f90

OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(SRC:.c=.o)

# Suffix rules
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) -c $<

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

t_mp:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) t_mp.f90 $(LIB_QE_MATH) -o t_mp.x

pwt_electrons_v2:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) cheFSI.f90 pwt_electrons_v2.f90 \
		$(LIB_QE_MATH) -o pwt_electrons_v2.x

pwt_electrons_v1:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) pwt_electrons_v1.f90 $(LIB_QE_MATH) \
		-o pwt_electrons_v1.x

pwt_plot_vltot:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) pwt_plot_vltot.f90 $(LIB_QE_MATH) \
		-o pwt_plot_vltot.x

t_pseudo_v1:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) t_pseudo_v1.f90 $(LIB_QE_MATH) \
		-o t_pseudo_v1.x

free_fcc:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) free_fcc.f90 $(LIB_QE_MATH) \
		-o free_fcc.x

t_gvect_v2:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) t_gvect_v2.f90 $(LIB_QE_MATH) \
	libmain.a	-o t_gvect_v2.x

coul_fcc:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) coul_fcc.f90 $(LIB_QE_MATH) \
		-o coul_fcc.x

t_gvect_v1:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) t_gvect_v1.f90 $(LIB_QE_MATH) \
		-o t_gvect_v1.x

t_struct_v1:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) t_struct_v1.f90 $(LIB_QE_MATH) \
		-o t_struct_v1.x

get_evc:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) -I $(INCLUDE_FFT) \
        get_evc.f90 $(LIB_QE_MATH) -o get_evc.x

pwt_uspp_param_v1:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) pwt_uspp_param_v1.f90 $(LIB_QE_MATH) \
		-o pwt_uspp_param_v1.x

t_uspp_param_r:
	$(F90) $(F90_OPTS) -I $(INCLUDE_MOD) -I $(INCLUDE_PW) t_uspp_param_r.f90 $(LIB_QE_MATH) \
		-o t_uspp_param_r.x

clean:
	rm -f *.o *.mod libmain.a x-main *.x

