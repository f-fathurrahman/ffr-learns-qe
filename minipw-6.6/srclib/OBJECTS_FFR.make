# !!! Use .o extension
OBJ_FFR = \
prepare_all.o my_addusdens.o my_addusdens_g.o \
test_loop_addusdens_g.o my_electrons.o my_electrons_scf.o \
my_compute_magnetization.o my_delta_e.o my_delta_escf.o my_calc_pol.o \
my_qvan2.o my_sum_band.o my_init_us_2.o my_cegterg.o \
my_h_psi.o my_s_psi.o my_vloc_psi_k.o my_c_bands.o \
prepare_h_s_psi.o my_newd.o my_add_vuspsi.o my_usnldiag.o \
my_sum_bec.o my_wfcinit.o my_init_wfc.o my_atomic_wfc.o \
my_potinit.o my_atomic_rho.o my_init_run.o \
my_v_of_rho.o my_v_h.o my_v_xc.o my_v_xc_meta.o my_set_vrs.o \
jl_comm_module.o my_hinit0.o my_setlocal.o \
my_mix_rho.o my_forces.o my_force_us.o my_addusforce.o my_compute_deff.o \
my_vloc_of_g.o my_dvloc_of_g.o \
my_force_lc.o my_force_corr.o my_force_cc.o \
my_stress.o my_stress_ewa.o my_stress_har.o my_stress_loc.o my_stress_cc.o \
my_stress_knl.o my_stress_us.o my_addusstress.o \
my_gradcorr.o \
info_upf.o write_array.o

# fft_interpolate is defined as an interface

# NOTE: my_becmod is deleted. There are many modules and subroutines that use it.
# We modify becmod directly.