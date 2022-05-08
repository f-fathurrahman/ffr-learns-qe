MAKE = make
AR = ar

#DFLAGS = -D__FFTW -D__MPI -D__LIBXC
DFLAGS = -D__FFTW -D__LIBXC

F77 = gfortran # mpif90
#F77_OPTS = -Wall -O3 -cpp
F77_OPTS = -Wall -cpp

F90 = gfortran # mpif90
#F90_OPTS = $(DFLAGS) -Wall -O3 -mtune=native -cpp # -g -fbacktrace -ffpe-trap=zero,overflow,underflow
F90_OPTS = $(DFLAGS) -Wall -cpp # -g -fbacktrace -ffpe-trap=zero,overflow,underflow

LIBXC_HOME = /home/efefer/mysoftwares/libxc-4.3.4/

OBJ_FFTXLIB = \
scatter_mod.o  \
fft_ggen.o  \
fft_fwinv.o  \
fft_scalar.o  \
fft_scalar.ARM_LIB.o  \
fft_scalar.DFTI.o  \
fft_scalar.ESSL.o  \
fft_scalar.FFTW.o  \
fftw_interfaces.o  \
fft_scalar.FFTW3.o  \
fft_scalar.SX6.o  \
fft_parallel.o  \
fft_interfaces.o  \
fft_interpolate.o \
stick_base.o  \
fftw.o  \
fftw_sp.o  \
fftw_dp.o  \
fft_smallbox.o  \
fft_smallbox_type.o  \
fft_support.o  \
fft_error.o  \
fft_stick.o  \
fft_types.o \
tg_gather.o \
fft_helper_subroutines.o \
fft_param.o


OBJ_LAXLIB = la_types.o \
la_error.o \
la_helper.o \
cdiaghg.o \
rdiaghg.o \
dspev_drv.o \
ptoolkit.o \
transto.o \
distools.o \
zhpev_drv.o \
mp_diag.o \
la_param.o


# Some UtilXLib objects are also here
OBJ_MODULES = \
parallel_include.o \
util_param.o \
mp.o \
data_buffer.o \
mp_base.o \
mp_bands.o \
io_global.o \
mp_world.o \
mp_images.o \
mp_pools.o \
command_line_options.o \
mp_exx.o \
mp_global.o \
kind.o \
mp_wave.o \
mp_bands_util.o  mp_base_gpu.o \
parameters.o \
constants.o \
wrappers.o \
io_files.o \
version.o \
environment.o \
parser.o \
wypos.o \
input_parameters.o \
invmat.o \
control_flags.o \
cell_base.o \
bz_form.o \
read_cards.o read_namelists.o \
random_numbers.o \
ions_base.o \
upf_kinds.o \
upf_const.o \
upf_params.o \
radial_grids.o \
atom.o \
pseudo_types.o \
upf_invmat.o \
uspp.o \
fft_base.o \
fft_rho.o \
correlation_lda_lsda.o correlation_gga.o \
exchange_gga.o  exchange_lda_lsda.o \
recvec.o recvec_subs.o \
xc_rVV10.o \
xc_vdW_DF.o \
funct.o \
xc_lda_lsda_drivers.o xc_gga_drivers.o \
metagga.o \
xc_mgga_drivers.o \
upf_utils.o \
xmltools.o \
read_upf_new.o \
read_upf_v1.o \
upf_io.o \
upf_auxtools.o \
upf_to_internal.o \
read_uspp.o \
gth.o \
read_pseudo.o \
open_close_input_file.o \
read_input.o \
error_handler.o \
latgen.o recips.o volume.o \
matches.o \
cryst_to_car.o \
sort.o capital.o \
find_free_unit.o \
test_input_file.o \
int_to_char.o \
generate_k_along_lines.o \
date_and_tim.o \
clocks_handler.o \
copy.o md5_from_file.o md5.o memstat.o c_mkdir.o eval_infix.o cptimer.o \
set_mpi_comm_4_solvers.o \
noncol.o \
becmod.o \
paw_variables.o \
wavefunctions.o \
splinelib.o \
gvecw.o \
check_stop.o \
coulomb_vcut.o \
basic_algebra_routines.o \
constraints_module.o \
mm_dispersion.o \
tsvdw.o \
io_base.o \
fsockets.o \
sockets.o \
fcp_variables.o \
plugin_variables.o \
ws_base.o \
mdiis.o \
plugin_flags.o \
qmmm.o \
run_info.o \
bfgs_module.o \
space_group.o \
wyckoff.o \
plot_io.o

OBJ_PW = \
a2fmod.o \
add_bfield.o \
add_efield.o \
add_vuspsi.o \
add_gatefield.o \
add_paw_to_deeq.o \
add_vhub_to_deeq.o \
addusdens.o \
addusforce.o \
addusstress.o \
allocate_fft.o \
allocate_locpot.o \
allocate_nlpot.o \
allocate_wfc.o \
atomic_rho.o \
atomic_wfc.o \
atomic_wfc_mod.o \
average_pp.o \
acfdt_in_pw.o \
newd.o \
beef.o \
bp_mod.o \
bp_c_phase.o \
bp_calc_btq.o \
bp_qvan3.o \
bp_strings.o \
buffers.o \
c_bands.o \
c_phase_field.o \
orbm_kubo.o \
cdiagh.o \
clean_pw.o \
close_files.o \
compute_becsum.o \
compute_deff.o \
compute_dip.o \
compute_rho.o \
compute_qdipol.o \
compute_qdipol_so.o \
compute_ux.o \
coset.o \
Coul_cut_2D.o \
d_matrix.o \
data_structure.o \
deriv_drhoc.o \
divide_class.o \
divide_class_so.o \
divide_et_impera.o \
dqvan2.o \
drhoc.o \
rotate_wfc.o \
run_driver.o \
dvloc_of_g.o \
dynamics_module.o \
efermig.o \
efermit.o \
electrons.o \
eqvect.o \
esm.o \
ewald.o \
ewald_dipole.o \
extfield.o \
exx_base.o \
exx_band.o \
exx.o \
fcp.o \
find_group.o \
forces_bp_efield.o \
force_cc.o \
force_corr.o \
force_ew.o \
force_hub.o \
force_lc.o \
force_us.o \
forces.o \
g_psi.o \
g_psi_mod.o \
gen_at_dj.o \
gen_at_dy.o \
gen_us_dj.o \
gen_us_dy.o \
get_locals.o \
gk_sort.o \
gradcorr.o \
gweights.o \
g2_kin.o \
hs_psi.o \
hs_1psi.o \
h_epsi_her_apply.o \
h_epsi_her_set.o \
h_psi.o \
h_psi_meta.o \
hinit0.o \
hinit1.o \
init_ns.o \
init_q_aeps.o \
init_run.o \
init_us_0.o \
init_us_b0.o \
init_us_1.o \
init_us_2.o \
init_at_1.o \
init_vloc.o \
input.o \
io_rho_xml.o \
irrek.o \
iweights.o \
intersite_V.o \
init_nsg.o \
nsg_adj.o \
start_k.o \
kpoint_grid.o \
lchk_tauxk.o \
ldaU.o \
make_pointlists.o \
makov_payne.o \
manypw.o \
martyna_tuckerman.o \
memory_report.o \
mix_rho.o \
move_ions.o \
multable.o \
n_plane_waves.o \
new_ns.o \
new_nsb.o \
new_nsg.o \
new_occ.o \
ns_adj.o \
non_scf.o \
offset_atom_wfc.o \
openfil.o \
orthoatwfc.o \
output_tau.o \
para.o \
paw_exx.o \
paw_init.o \
paw_onecenter.o \
paw_symmetry.o \
plugin_print_energies.o \
plugin_scf_energy.o \
plugin_scf_potential.o \
plugin_init_ions.o \
plugin_init_cell.o \
plugin_init_potential.o \
plugin_initbase.o \
plugin_clean.o \
plugin_check.o \
plugin_clock.o \
plugin_summary.o \
plugin_initialization.o \
plugin_ext_forces.o \
plugin_int_forces.o \
plugin_read_input.o \
plus_u_full.o \
potinit.o \
print_clock_pw.o \
print_ks_energies.o \
punch.o \
electrons_base.o \
pwcom.o \
pw2blip.o \
pw2casino.o \
pw2casino_write.o \
qvan2.o \
rdiagh.o \
read_file_new.o \
realus.o \
remove_atomic_rho.o \
report_mag.o \
restart_in_electrons.o \
rho2zeta.o \
ruotaijk.o \
run_pwscf.o \
s_1psi.o \
s_psi.o \
save_in_cbands.o \
save_in_electrons.o \
scale_h.o \
loc_scdm.o \
loc_scdm_k.o \
scf_mod.o \
set_kplusq.o \
set_kup_and_kdw.o \
set_occupations.o \
set_rhoc.o \
set_spin_vars.o \
set_vrs.o \
setlocal.o \
setup.o \
spinor.o \
sph_ind.o \
stop_run.o \
stres_cc.o \
stres_ewa.o \
stres_gradcorr.o \
stres_har.o \
stres_hub.o \
stres_knl.o \
stres_loc.o \
stres_us.o \
stres_nonloc_dft.o \
stres_mgga.o \
stress.o \
struct_fact.o \
sum_band.o \
sumkg.o \
sumkt.o \
summary.o \
symme.o \
symm_base.o \
symmetrize_at.o \
tabd.o \
tetra.o \
transform_becsum_so.o \
transform_becsum_nc.o \
transform_qq_so.o \
trnvecc.o \
update_pot.o \
us_exx.o \
usnldiag.o \
v_of_rho.o \
vcsmd.o \
vcsubs.o \
vhpsi.o \
vloc_of_g.o \
vloc_psi.o \
utils.o \
xdm_dispersion.o \
wfcinit.o \
write_ns.o \
wsweight.o \
weights.o \
ortho_wfc.o \
wgauss.o w1gauss.o \
simpsn.o \
erf.o \
sph_bes.o sph_dbes.o \
ylmr2.o dylmr2.o \
divide.o \
gradutils.o \
thread_util.o \
upf_parallel_include.o \
upf_error.o \
dmxc_drivers.o \
dgcxc_drivers.o \
radial_gradients.o \
rgen.o \
atomic_number.o \
set_hubbard_l.o  set_hubbard_n.o \
trimcheck.o \
atom_weight.o \
set_vdw_corr.o \
david_rci.o \
cegterg.o \
regterg.o \
ccgdiagg.o \
rcgdiagg.o \
ppcg_gamma.o \
ppcg_k.o \
bpcg_gamma.o \
bpcg_k.o \
paro_gamma.o \
pcg_gamma.o \
paro_k.o \
pcg_k.o \
paro_gamma_new.o \
paro_k_new.o \
rotate_HSpsi_gamma.o \
rotate_HSpsi_k.o \
rotate_wfc_gamma.o \
rotate_wfc_k.o \
export_gstart_2_solvers.o \
deviatoric.o \
compute_dipole.o \
remove_tot_torque.o \
linpack.o

OBJ_DFTD3 = \
dftd3_common.o \
dftd3_sizes.o \
dftd3_pars.o \
dftd3_core.o \
dftd3_api.o \
dftd3_qe.o 

OBJ = $(OBJ_FFTXLIB) \
$(OBJ_LAXLIB) \
$(OBJ_MODULES) \
$(OBJ_DFTD3) \
$(OBJ_PW) \
$(SRC:.f90=.o) $(SRC:.f=.o)

#
# Suffix rule for Fortran 90
#
%.mod :
	@if [! -f $@ ]; then \
		rm $(*F).o; \
		fi
	$(MAKE) $<

%.o : %.f90
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

#
# Fortran 77 sources
# supress warning
.SUFFIXES: .o .f
.f.o:
	$(F77) $(F77_OPTS) -c $<


# Targets
lib: $(OBJ)
	ar rcs libmain.a *.o

clean:
	rm -rf *.o *.mod libmain.a *.x


#main: lib
#	mpif90 -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  minipw.f90 -o minipw.x libmain.a /home/efefer/mysoftwares/libxc-4.3.4/lib/libxcf03.a /home/efefer/mysoftwares/libxc-4.3.4/lib/libxc.a -lblas -llapack

#pw: lib
#	mpif90 -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  pwscf.f90 -o pw.x libmain.a /home/efefer/mysoftwares/libxc-4.3.4/lib/libxcf03.a /home/efefer/mysoftwares/libxc-4.3.4/lib/libxc.a -L/home/efefer/mysoftwares/LIBMKL -lmkl_gf_lp64  -lmkl_sequential -lmkl_core

clean_pw.o:
	mpif90 -I -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o clean_pw.o clean_pw.f90

electrons.o:
	mpif90 -I -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o electrons.o electrons.f90

input.o:
	mpif90 -I -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o input.o input.f90

forces.o:
	mpif90 -I -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o forces.o forces.f90

stress.o:
	mpif90 -I -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o stress.o stress.f90

funct.o:
	mpif90 -I $(LIBXC_HOME)/include -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o funct.o funct.f90

xc_lda_lsda_drivers.o:
	mpif90 -I $(LIBXC_HOME)/include -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o xc_lda_lsda_drivers.o xc_lda_lsda_drivers.f90

xc_gga_drivers.o:
	mpif90 -I $(LIBXC_HOME)/include -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o xc_gga_drivers.o xc_gga_drivers.f90

xc_mgga_drivers.o:
	mpif90 -I $(LIBXC_HOME)/include -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o xc_mgga_drivers.o xc_mgga_drivers.f90

dmxc_drivers.o:
	mpif90 -I $(LIBXC_HOME)/include -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o dmxc_drivers.o dmxc_drivers.f90

dgcxc_drivers.o:
	mpif90 -I $(LIBXC_HOME)/include -D__FFTW -D__MPI -D__LIBXC -Wall -O3 -cpp  -c -o dgcxc_drivers.o dgcxc_drivers.f90

fftw.o:
	gcc -O3 -c -o fftw.o fftw.c

fftw_sp.o:
	gcc -O3 -c -o fftw_sp.o fftw_sp.c

fftw_dp.o:
	gcc -O3 -c -o fftw_dp.o fftw_dp.c

fftw_stick.o:
	gcc -O3 -c -o fftw_stick.o fftw_stick.c


include DEPENDS