a2fmod.o : a2fmod.f90 io_files.o io_global.o symm_base.o start_k.o pwcom.o ions_base.o pwcom.o kind.o 
acfdt_in_pw.o : acfdt_in_pw.f90 kind.o 
add_bfield.o : add_bfield.f90 noncol.o mp.o mp_bands.o pwcom.o fft_base.o cell_base.o ions_base.o io_global.o constants.o kind.o 
add_efield.o : add_efield.f90 mp.o fft_types.o fft_base.o mp_bands.o mp_images.o pwcom.o control_flags.o io_global.o pwcom.o extfield.o cell_base.o ions_base.o constants.o kind.o 
add_gatefield.o : add_gatefield.f90 mp.o fft_types.o fft_base.o mp_bands.o mp_images.o pwcom.o control_flags.o io_global.o pwcom.o pwcom.o extfield.o cell_base.o ions_base.o constants.o kind.o 
add_paw_to_deeq.o : add_paw_to_deeq.f90 pwcom.o paw_variables.o uspp.o ions_base.o kind.o 
add_qexsd_step.o : add_qexsd_step.f90 extfield.o fcp_variables.o control_flags.o pwcom.o pwcom.o pwcom.o cell_base.o ions_base.o constants.o kind.o 
addusdens.o : addusdens.f90 mp.o mp_bands.o mp_pools.o uspp.o uspp.o recvec.o fft_interfaces.o cell_base.o ions_base.o kind.o fft_base.o noncol.o control_flags.o realus.o 
addusforce.o : addusforce.f90 fft_interfaces.o mp.o mp_pools.o mp_bands.o uspp.o uspp.o scf_mod.o noncol.o recvec.o fft_base.o cell_base.o realus.o control_flags.o ions_base.o kind.o 
addusstress.o : addusstress.f90 mp.o mp_pools.o fft_interfaces.o uspp.o uspp.o scf_mod.o pwcom.o recvec.o fft_base.o cell_base.o ions_base.o realus.o control_flags.o kind.o 
add_vhub_to_deeq.o : add_vhub_to_deeq.f90 ldaU.o scf_mod.o pwcom.o uspp.o ions_base.o kind.o 
add_vuspsi.o : add_vuspsi.f90 mp.o becmod.o uspp.o uspp.o noncol.o control_flags.o pwcom.o ions_base.o kind.o 
allocate_fft.o : allocate_fft.f90 funct.o wavefunctions.o noncol.o control_flags.o scf_mod.o pwcom.o ions_base.o fft_base.o recvec.o recvec.o io_global.o 
allocate_locpot.o : allocate_locpot.f90 fft_base.o recvec.o pwcom.o ions_base.o 
allocate_nlpot.o : allocate_nlpot.f90 pwcom.o uspp.o uspp.o pwcom.o gvecw.o recvec.o noncol.o pwcom.o pwcom.o pwcom.o ions_base.o control_flags.o 
allocate_wfc.o : allocate_wfc.f90 pwcom.o recvec.o gvecw.o uspp.o wavefunctions.o noncol.o ldaU.o pwcom.o atomic_wfc_mod.o pwcom.o io_global.o 
atom.o : atom.f90 radial_grids.o 
atomic_number.o : atomic_number.f90 upf_utils.o 
atomic_rho.o : atomic_rho.f90 fft_rho.o fft_base.o mp.o mp_bands.o control_flags.o io_global.o uspp.o noncol.o pwcom.o pwcom.o recvec.o cell_base.o ions_base.o atom.o constants.o kind.o 
atomic_wfc.o : atomic_wfc.f90 mp.o mp_bands.o pwcom.o noncol.o uspp.o pwcom.o pwcom.o pwcom.o recvec.o atomic_wfc_mod.o ions_base.o cell_base.o constants.o kind.o 
atomic_wfc_mod.o : atomic_wfc_mod.f90 kind.o 
atom_weight.o : atom_weight.f90 kind.o 
average_pp.o : average_pp.f90 uspp.o atom.o kind.o 
basic_algebra_routines.o : basic_algebra_routines.f90 kind.o 
becmod.o : becmod.f90 mp.o mp_bands.o noncol.o recvec.o control_flags.o kind.o 
beef.o : beef.f90 scf_mod.o pwcom.o control_flags.o funct.o io_global.o kind.o 
bfgs_module.o : bfgs_module.f90 invmat.o basic_algebra_routines.o cell_base.o constants.o io_files.o kind.o 
bp_calc_btq.o : bp_calc_btq.f90 uspp.o constants.o cell_base.o ions_base.o atom.o kind.o 
bpcg_gamma.o : bpcg_gamma.f90 mp.o mp_bands_util.o util_param.o 
bpcg_k.o : bpcg_k.f90 mp.o mp_bands_util.o util_param.o 
bp_c_phase.o : bp_c_phase.f90 mp.o mp_bands.o pwcom.o noncol.o becmod.o bp_mod.o wavefunctions.o pwcom.o pwcom.o pwcom.o uspp.o uspp.o fft_base.o recvec.o constants.o cell_base.o ions_base.o buffers.o io_files.o io_global.o kind.o 
bp_mod.o : bp_mod.f90 cell_base.o fft_base.o mp_images.o mp.o recvec.o becmod.o kind.o 
bp_qvan3.o : bp_qvan3.f90 uspp.o uspp.o ions_base.o kind.o 
bp_strings.o : bp_strings.f90 kind.o 
buffers.o : buffers.f90 io_files.o kind.o 
bz_form.o : bz_form.f90 constants.o io_global.o kind.o 
capital.o : capital.f90 
c_bands.o : c_bands.f90 atomic_wfc_mod.o mp_bands.o becmod.o scf_mod.o g_psi_mod.o noncol.o check_stop.o mp.o mp_pools.o bp_mod.o wavefunctions.o pwcom.o ldaU.o control_flags.o pwcom.o recvec.o uspp.o pwcom.o buffers.o io_files.o io_global.o kind.o 
ccgdiagg.o : ccgdiagg.f90 mp.o mp_bands_util.o util_param.o 
cdiagh.o : cdiagh.f90 mp.o mp_bands.o kind.o 
cdiaghg.o : cdiaghg.f90 laxlib_low.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh zhpev_drv.o mp_diag.o la_types.o la_param.o 
cegterg.o : cegterg.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
cell_base.o : cell_base.f90 control_flags.o invmat.o io_global.o constants.o kind.o 
check_stop.o : check_stop.f90 mp_images.o mp.o io_files.o io_global.o clocks_handler.o kind.o 
clean_pw.o : clean_pw.f90 dftd3_qe.o tsvdw.o control_flags.o Coul_cut_2D.o exx.o bp_mod.o pseudo_types.o realus.o constraints_module.o xdm_dispersion.o mm_dispersion.o radial_grids.o atom.o paw_init.o dynamics_module.o noncol.o pwcom.o fft_types.o fft_base.o extfield.o ldaU.o gth.o uspp.o uspp.o pwcom.o wavefunctions.o symme.o symm_base.o scf_mod.o pwcom.o pwcom.o pwcom.o recvec.o pwcom.o tetra.o pwcom.o ions_base.o pwcom.o atomic_wfc_mod.o 
clocks_handler.o : clocks_handler.f90 parallel_include.o util_param.o 
close_files.o : close_files.f90 bp_mod.o mp.o mp_images.o buffers.o io_files.o pwcom.o control_flags.o ldaU.o 
command_line_options.o : command_line_options.f90 io_global.o mp_world.o mp.o 
compute_becsum.o : compute_becsum.f90 becmod.o paw_variables.o paw_symmetry.o mp.o mp_bands.o mp_pools.o pwcom.o wavefunctions.o uspp.o scf_mod.o buffers.o io_files.o pwcom.o pwcom.o control_flags.o kind.o 
compute_deff.o : compute_deff.f90 noncol.o pwcom.o pwcom.o uspp.o uspp.o ions_base.o kind.o 
compute_dip.o : compute_dip.f90 mp.o mp_bands.o fft_types.o fft_base.o pwcom.o extfield.o pwcom.o cell_base.o kind.o constants.o ions_base.o io_global.o 
compute_dipole.o : compute_dipole.f90 mp.o mp_bands.o fft_types.o fft_base.o cell_base.o kind.o 
compute_qdipol.o : compute_qdipol.f90 uspp.o uspp.o ions_base.o atom.o constants.o kind.o 
compute_qdipol_so.o : compute_qdipol_so.f90 pwcom.o uspp.o pwcom.o ions_base.o kind.o 
compute_rho.o : compute_rho.f90 noncol.o kind.o 
compute_ux.o : compute_ux.f90 noncol.o io_global.o constants.o kind.o 
constants.o : constants.f90 kind.o 
constraints_module.o : constraints_module.f90 cell_base.o ions_base.o input_parameters.o basic_algebra_routines.o io_global.o constants.o kind.o 
control_flags.o : control_flags.f90 parameters.o kind.o 
correlation_gga.o : correlation_gga.f90 kind.o correlation_lda_lsda.o 
correlation_lda_lsda.o : correlation_lda_lsda.f90 kind.o 
coset.o : coset.f90 kind.o 
Coul_cut_2D.o : Coul_cut_2D.f90 pwcom.o uspp.o ions_base.o fft_base.o cell_base.o recvec.o io_global.o constants.o kind.o 
coulomb_vcut.o : coulomb_vcut.f90 
c_phase_field.o : c_phase_field.f90 pwcom.o becmod.o mp_pools.o mp_bands.o mp.o pwcom.o bp_mod.o noncol.o pwcom.o pwcom.o pwcom.o uspp.o uspp.o recvec.o fft_base.o constants.o cell_base.o ions_base.o buffers.o io_files.o io_global.o kind.o 
cryst_to_car.o : cryst_to_car.f90 kind.o 
data_buffer.o : data_buffer.f90 util_param.o 
data_structure.o : data_structure.f90 command_line_options.o symm_base.o realus.o io_global.o gvecw.o recvec.o recvec.o pwcom.o cell_base.o fft_types.o fft_base.o mp_pools.o mp_bands.o mp.o kind.o 
date_and_tim.o : date_and_tim.f90 
david_rci.o : david_rci.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
deriv_drhoc.o : deriv_drhoc.f90 constants.o kind.o 
deviatoric.o : deviatoric.f90 io_global.o kind.o 
dftd3_api.o : dftd3_api.f90 dftd3_core.o dftd3_common.o dftd3_sizes.o 
dftd3_common.o : dftd3_common.f90 
dftd3_core.o : dftd3_core.f90 dftd3_pars.o dftd3_common.o dftd3_sizes.o 
dftd3_pars.o : dftd3_pars.f90 dftd3_common.o dftd3_sizes.o 
dftd3_qe.o : dftd3_qe.f90 cell_base.o ions_base.o io_global.o dftd3_api.o dftd3_core.o dftd3_common.o dftd3_sizes.o 
dftd3_sizes.o : dftd3_sizes.f90 
dftd3_test_code.o : dftd3_test_code.f90 cell_base.o dftd3_api.o 
dgcxc_drivers.o : dgcxc_drivers.f90 xc_gga_drivers.o funct.o kind.o constants.o 
distools.o : distools.f90 
divide_class.o : divide_class.f90 constants.o kind.o 
divide_class_so.o : divide_class_so.f90 constants.o io_global.o noncol.o pwcom.o pwcom.o pwcom.o pwcom.o kind.o 
divide_et_impera.o : divide_et_impera.f90 mp_pools.o kind.o 
divide.o : divide.f90 mp.o 
d_matrix.o : d_matrix.f90 invmat.o random_numbers.o symm_base.o kind.o 
dmxc_drivers.o : dmxc_drivers.f90 constants.o correlation_lda_lsda.o exchange_lda_lsda.o xc_lda_lsda_drivers.o funct.o kind.o 
dqvan2.o : dqvan2.f90 uspp.o uspp.o pwcom.o kind.o 
drhoc.o : drhoc.f90 constants.o kind.o 
dspev_drv.o : dspev_drv.f90 laxlib_kinds.fh la_param.o 
dvloc_of_g.o : dvloc_of_g.f90 esm.o Coul_cut_2D.o constants.o kind.o 
dylmr2.o : dylmr2.f90 kind.o 
dynamics_module.o : dynamics_module.f90 extfield.o pwcom.o symm_base.o random_numbers.o constraints_module.o pwcom.o pwcom.o cell_base.o basic_algebra_routines.o control_flags.o constants.o io_files.o io_global.o ions_base.o kind.o 
efermig.o : efermig.f90 mp_pools.o mp.o constants.o kind.o io_global.o 
efermit.o : efermit.f90 constants.o kind.o io_global.o 
electrons_base.o : electrons_base.f90 io_global.o constants.o kind.o 
electrons.o : electrons.f90 wavefunctions.o buffers.o constants.o plugin_variables.o wrappers.o fcp_variables.o esm.o newd.o xdm_dispersion.o dftd3_qe.o dftd3_api.o mm_dispersion.o mp.o mp_pools.o mp_bands.o io_rho_xml.o pwcom.o noncol.o gvecw.o pwcom.o pwcom.o bp_mod.o atomic_wfc_mod.o cell_base.o loc_scdm_k.o loc_scdm.o ions_base.o paw_symmetry.o paw_onecenter.o paw_variables.o funct.o exx.o uspp.o pwcom.o pwcom.o extfield.o ldaU.o io_files.o control_flags.o scf_mod.o tsvdw.o pwcom.o pwcom.o recvec.o recvec.o fft_base.o io_global.o check_stop.o kind.o 
environment.o : environment.f90 command_line_options.o version.o mp_bands.o mp_pools.o mp_images.o mp_world.o io_global.o io_files.o kind.o 
eqvect.o : eqvect.f90 kind.o 
erf.o : erf.f90 kind.o 
error_handler.o : error_handler.f90 mp.o util_param.o 
esm.o : esm.f90 io_files.o pwcom.o io_global.o fft_scalar.o mp.o mp_bands.o control_flags.o cell_base.o constants.o ions_base.o fft_base.o recvec.o kind.o 
ewald_dipole.o : ewald_dipole.f90 mp.o mp_bands.o pwcom.o ions_base.o cell_base.o constants.o recvec.o kind.o 
ewald.o : ewald.f90 Coul_cut_2D.o martyna_tuckerman.o mp.o mp_bands.o constants.o kind.o 
exchange_gga.o : exchange_gga.f90 kind.o 
exchange_lda_lsda.o : exchange_lda_lsda.f90 constants.o kind.o 
export_gstart_2_solvers.o : export_gstart_2_solvers.f90 mp_bands_util.o 
extfield.o : extfield.f90 kind.o 
exx_band.o : exx_band.f90 command_line_options.o recvec_subs.o fft_base.o mp_bands.o gvecw.o recvec.o pwcom.o cell_base.o parallel_include.o recvec.o mp.o mp_pools.o mp_exx.o buffers.o wavefunctions.o uspp.o pwcom.o io_files.o pwcom.o stick_base.o fft_types.o control_flags.o io_global.o noncol.o kind.o 
exx_base.o : exx_base.f90 mp_exx.o gvecw.o recvec.o constants.o funct.o start_k.o pwcom.o pwcom.o cell_base.o symm_base.o pwcom.o mp.o mp_pools.o mp_pools.o mp_images.o stick_base.o fft_types.o control_flags.o io_global.o noncol.o coulomb_vcut.o kind.o 
exx.o : exx.f90 uspp.o splinelib.o gth.o pwcom.o ions_base.o coulomb_vcut.o pwcom.o constants.o mp_pools.o paw_exx.o paw_variables.o uspp.o fft_interfaces.o scatter_mod.o funct.o pwcom.o buffers.o io_files.o wavefunctions.o exx_base.o us_exx.o becmod.o command_line_options.o exx_band.o realus.o mp.o mp_pools.o pwcom.o mp_bands.o mp_exx.o symm_base.o fft_base.o recvec_subs.o cell_base.o recvec.o gvecw.o stick_base.o fft_types.o control_flags.o io_global.o noncol.o kind.o 
fcp.o : fcp.f90 symm_base.o random_numbers.o ions_base.o pwcom.o pwcom.o cell_base.o mdiis.o fcp_variables.o dynamics_module.o control_flags.o constants.o io_files.o io_global.o kind.o 
fcp_variables.o : fcp_variables.f90 kind.o 
fft_base.o : fft_base.f90 stick_base.o fft_smallbox_type.o fft_types.o parallel_include.o 
fft_error.o : fft_error.f90 fft_param.o 
fft_fwinv.o : fft_fwinv.f90 fft_smallbox_type.o fft_param.o fft_types.o fft_parallel.o fft_smallbox.o fft_scalar.o 
fft_ggen.o : fft_ggen.f90 fft_types.o fft_param.o 
fft_helper_subroutines.o : fft_helper_subroutines.f90 fft_types.o fft_param.o 
fft_interfaces.o : fft_interfaces.f90 fft_smallbox_type.o fft_param.o fft_types.o 
fft_interpolate.o : fft_interpolate.f90 fft_interfaces.o fft_types.o fft_param.o 
fft_parallel.o : fft_parallel.f90 fft_types.o scatter_mod.o fft_scalar.o fft_param.o 
fft_param.o : fft_param.f90 
fft_rho.o : fft_rho.f90 fft_helper_subroutines.o fft_types.o control_flags.o fft_interfaces.o kind.o 
fft_scalar.ARM_LIB.o : fft_scalar.ARM_LIB.f90 fftw_interfaces.o fft_param.o 
fft_scalar.DFTI.o : fft_scalar.DFTI.f90 fft_param.o 
fft_scalar.ESSL.o : fft_scalar.ESSL.f90 fft_param.o 
fft_scalar.o : fft_scalar.f90 fft_scalar.FFTW.o fft_scalar.ARM_LIB.o fft_scalar.SX6.o fft_scalar.ESSL.o fft_scalar.DFTI.o fft_scalar.FFTW3.o fft_param.o 
fft_scalar.FFTW3.o : fft_scalar.FFTW3.f90 fft_param.o 
fft_scalar.FFTW.o : fft_scalar.FFTW.f90 fftw_interfaces.o fft_param.o 
fft_scalar.SX6.o : fft_scalar.SX6.f90 fft_param.o 
fft_smallbox.o : fft_smallbox.f90 fftw_interfaces.o 
fft_smallbox_type.o : fft_smallbox_type.f90 fft_types.o 
fft_support.o : fft_support.f90 fft_param.o 
fft_types.o : fft_types.f90 stick_base.o fft_param.o fft_support.o 
fftw_interfaces.o : fftw_interfaces.f90 
find_free_unit.o : find_free_unit.f90 
find_group.o : find_group.f90 kind.o 
force_cc.o : force_cc.f90 mp.o mp_bands.o wavefunctions.o noncol.o control_flags.o scf_mod.o pwcom.o pwcom.o recvec.o fft_interfaces.o fft_base.o cell_base.o ions_base.o uspp.o atom.o constants.o kind.o 
force_corr.o : force_corr.f90 mp.o mp_bands.o wavefunctions.o control_flags.o scf_mod.o pwcom.o recvec.o fft_interfaces.o fft_base.o cell_base.o ions_base.o uspp.o atom.o constants.o kind.o 
force_ew.o : force_ew.f90 Coul_cut_2D.o mp.o mp_bands.o constants.o kind.o 
force_hub.o : force_hub.f90 io_global.o recvec.o pwcom.o noncol.o mp_bands.o buffers.o pwcom.o wavefunctions.o uspp.o uspp.o becmod.o mp.o mp_pools.o scf_mod.o pwcom.o control_flags.o pwcom.o io_files.o symme.o atomic_wfc_mod.o ldaU.o cell_base.o ions_base.o kind.o 
force_lc.o : force_lc.f90 Coul_cut_2D.o esm.o fft_interfaces.o fft_base.o mp.o mp_bands.o constants.o kind.o 
forces_bp_efield.o : forces_bp_efield.f90 parallel_include.o clocks_handler.o pwcom.o noncol.o becmod.o mp_bands.o mp_world.o mp.o pwcom.o bp_mod.o wavefunctions.o pwcom.o pwcom.o pwcom.o uspp.o uspp.o fft_base.o recvec.o constants.o buffers.o io_files.o io_global.o ions_base.o cell_base.o kind.o 
forces.o : forces.f90 qmmm.o esm.o tsvdw.o xdm_dispersion.o dftd3_qe.o dftd3_api.o mm_dispersion.o martyna_tuckerman.o uspp.o bp_mod.o plugin_flags.o control_flags.o extfield.o ldaU.o scf_mod.o pwcom.o pwcom.o symme.o pwcom.o recvec.o fft_base.o ions_base.o cell_base.o io_global.o kind.o 
force_us.o : force_us.f90 mp.o mp_bands.o mp_pools.o becmod.o buffers.o io_files.o pwcom.o noncol.o wavefunctions.o symme.o pwcom.o pwcom.o uspp.o uspp.o recvec.o pwcom.o ions_base.o cell_base.o control_flags.o kind.o 
fsockets.o : fsockets.f90 
funct.o : funct.f90 xc_rVV10.o xc_vdW_DF.o kind.o io_global.o 
g2_kin.o : g2_kin.f90 pwcom.o gvecw.o recvec.o pwcom.o cell_base.o kind.o 
gen_at_dj.o : gen_at_dj.f90 atomic_wfc_mod.o uspp.o pwcom.o pwcom.o recvec.o pwcom.o cell_base.o ions_base.o atom.o constants.o io_global.o kind.o 
gen_at_dy.o : gen_at_dy.f90 atomic_wfc_mod.o uspp.o pwcom.o pwcom.o recvec.o pwcom.o cell_base.o ions_base.o atom.o constants.o io_global.o kind.o 
generate_k_along_lines.o : generate_k_along_lines.f90 kind.o 
gen_us_dj.o : gen_us_dj.f90 uspp.o splinelib.o gth.o pwcom.o uspp.o pwcom.o recvec.o pwcom.o cell_base.o ions_base.o constants.o kind.o 
gen_us_dy.o : gen_us_dy.f90 uspp.o splinelib.o pwcom.o uspp.o pwcom.o recvec.o pwcom.o cell_base.o ions_base.o constants.o io_global.o kind.o 
get_locals.o : get_locals.f90 noncol.o fft_base.o mp.o mp_bands.o pwcom.o cell_base.o ions_base.o kind.o 
gk_sort.o : gk_sort.f90 pwcom.o constants.o kind.o 
g_psi.o : g_psi.f90 noncol.o g_psi_mod.o kind.o 
g_psi_mod.o : g_psi_mod.f90 kind.o 
gradcorr.o : gradcorr.f90 fft_rho.o fft_interfaces.o fft_base.o pwcom.o xc_gga_drivers.o funct.o cell_base.o pwcom.o recvec.o kind.o constants.o 
gradutils.o : gradutils.f90 fft_helper_subroutines.o fft_types.o fft_interfaces.o cell_base.o recvec.o fft_base.o kind.o 
gth.o : gth.f90 pseudo_types.o upf_params.o upf_const.o upf_kinds.o 
gvecw.o : gvecw.f90 mp.o kind.o 
gweights.o : gweights.f90 kind.o 
h_epsi_her_apply.o : h_epsi_her_apply.f90 mp.o mp_bands.o becmod.o io_global.o pwcom.o constants.o ions_base.o cell_base.o pwcom.o bp_mod.o uspp.o uspp.o recvec.o scf_mod.o pwcom.o ldaU.o pwcom.o pwcom.o kind.o noncol.o 
h_epsi_her_set.o : h_epsi_her_set.f90 becmod.o mp_bands.o mp.o pwcom.o constants.o buffers.o io_files.o ions_base.o cell_base.o pwcom.o bp_mod.o uspp.o uspp.o fft_base.o recvec.o scf_mod.o pwcom.o ldaU.o pwcom.o kind.o pwcom.o noncol.o 
hinit0.o : hinit0.f90 noncol.o io_global.o control_flags.o ldaU.o realus.o pwcom.o recvec.o fft_base.o pwcom.o cell_base.o atomic_wfc_mod.o ions_base.o kind.o 
hinit1.o : hinit1.f90 newd.o paw_symmetry.o paw_onecenter.o paw_variables.o martyna_tuckerman.o realus.o control_flags.o scf_mod.o noncol.o pwcom.o ldaU.o recvec.o fft_base.o cell_base.o ions_base.o 
h_psi.o : h_psi.f90 fft_helper_subroutines.o exx.o fft_base.o realus.o control_flags.o recvec.o ldaU.o uspp.o pwcom.o scf_mod.o pwcom.o becmod.o bp_mod.o mp.o mp_bands.o funct.o noncol.o kind.o 
h_psi_meta.o : h_psi_meta.f90 fft_interfaces.o fft_base.o wavefunctions.o control_flags.o pwcom.o scf_mod.o recvec.o pwcom.o pwcom.o cell_base.o kind.o 
hs_1psi.o : hs_1psi.f90 realus.o noncol.o bp_mod.o control_flags.o kind.o 
hs_psi.o : hs_psi.f90 noncol.o kind.o 
info_upf.o : info_upf.f90 uspp.o ions_base.o 
init_at_1.o : init_at_1.f90 mp.o mp_bands.o uspp.o pwcom.o ions_base.o cell_base.o constants.o atom.o kind.o 
init_ns.o : init_ns.f90 noncol.o uspp.o scf_mod.o ldaU.o pwcom.o ions_base.o kind.o 
init_nsg.o : init_nsg.f90 ldaU.o pwcom.o uspp.o ions_base.o kind.o 
init_q_aeps.o : init_q_aeps.f90 io_global.o control_flags.o uspp.o uspp.o ldaU.o pwcom.o atom.o ions_base.o kind.o 
init_run.o : init_run.f90 uspp.o uspp.o ions_base.o Coul_cut_2D.o tsvdw.o esm.o newd.o recvec_subs.o funct.o fft_base.o bp_mod.o paw_init.o paw_variables.o dynamics_module.o pwcom.o cell_base.o recvec.o recvec.o control_flags.o pwcom.o symme.o pwcom.o 
init_us_0.o : init_us_0.f90 mp.o mp_bands.o uspp.o pwcom.o cell_base.o ions_base.o atom.o constants.o io_global.o recvec.o kind.o 
init_us_1.o : init_us_1.f90 mp.o mp_bands.o paw_variables.o pwcom.o uspp.o uspp.o splinelib.o pwcom.o pwcom.o recvec.o cell_base.o ions_base.o atom.o constants.o kind.o 
init_us_2.o : init_us_2.f90 uspp.o uspp.o splinelib.o gth.o pwcom.o pwcom.o recvec.o constants.o cell_base.o ions_base.o kind.o 
init_us_b0.o : init_us_b0.f90 mp.o mp_bands.o uspp.o pwcom.o ions_base.o atom.o constants.o io_global.o gvecw.o kind.o 
init_vloc.o : init_vloc.f90 Coul_cut_2D.o recvec.o pwcom.o cell_base.o ions_base.o uspp.o kind.o gth.o atom.o 
input.o : input.f90 wyckoff.o recvec.o recvec.o pwcom.o tsvdw.o xdm_dispersion.o dftd3_qe.o dftd3_api.o mm_dispersion.o read_namelists.o constraints_module.o input_parameters.o pwcom.o qmmm.o read_pseudo.o realus.o Coul_cut_2D.o bfgs_module.o symm_base.o pwcom.o noncol.o pwcom.o gvecw.o pwcom.o check_stop.o update_pot.o pwcom.o pwcom.o loc_scdm.o exx.o exx_base.o a2fmod.o esm.o martyna_tuckerman.o ldaU.o start_k.o tetra.o pwcom.o fft_base.o pwcom.o io_files.o extfield.o fcp_variables.o dynamics_module.o pwcom.o run_info.o atomic_wfc_mod.o ions_base.o cell_base.o bp_mod.o io_global.o mp_pools.o constants.o control_flags.o funct.o kind.o 
input_parameters.o : input_parameters.f90 parameters.o kind.o 
intersite_V.o : intersite_V.f90 mp.o mp_images.o io_files.o pwcom.o fft_base.o ldaU.o control_flags.o parameters.o constants.o kind.o cell_base.o ions_base.o io_global.o symm_base.o 
int_to_char.o : int_to_char.f90 
invmat.o : invmat.f90 kind.o 
io_base.o : io_base.f90 io_global.o recvec.o mp.o mp_wave.o kind.o 
io_files.o : io_files.f90 wrappers.o mp_images.o mp.o io_global.o kind.o parameters.o 
io_global.o : io_global.f90 
ions_base.o : ions_base.f90 random_numbers.o cell_base.o io_global.o constants.o parameters.o kind.o 
io_rho_xml.o : io_rho_xml.f90 mp.o mp_images.o mp_bands.o mp_pools.o io_global.o control_flags.o recvec.o cell_base.o scf_mod.o pwcom.o noncol.o funct.o ldaU.o paw_variables.o io_base.o io_files.o kind.o 
irrek.o : irrek.f90 kind.o 
iweights.o : iweights.f90 mp_pools.o mp.o noncol.o kind.o 
jl_comm_module.o : jl_comm_module.f90 becmod.o kind.o 
kind.o : kind.f90 
kpoint_grid.o : kpoint_grid.f90 invmat.o noncol.o io_global.o bp_mod.o kind.o 
la_error.o : la_error.f90 la_param.o 
la_helper.o : la_helper.f90 laxlib_mid.fh laxlib_kinds.fh laxlib_low.fh laxlib_param.fh dspev_drv.o la_types.o la_param.o mp_diag.o 
la_param.o : la_param.f90 
latgen.o : latgen.f90 io_global.o constants.o kind.o 
la_types.o : la_types.f90 laxlib_param.fh 
lchk_tauxk.o : lchk_tauxk.f90 kind.o 
ld1x_add_exchange.o : ld1x_add_exchange.f90 radial_grids.o ld1x_ld1inc.o ld1x_parameters.o constants.o kind.o io_global.o 
ld1x_all_electron.o : ld1x_all_electron.f90 ld1x_ld1inc.o radial_grids.o kind.o 
ld1x_ascheq.o : ld1x_ascheq.f90 radial_grids.o kind.o 
ld1x_ascheqps_drv.o : ld1x_ascheqps_drv.f90 ld1x_ld1inc.o radial_grids.o ld1x_parameters.o kind.o 
ld1x_ascheqps.o : ld1x_ascheqps.f90 radial_grids.o kind.o io_global.o 
ld1x_atomic_paw.o : ld1x_atomic_paw.f90 mp_world.o mp.o io_global.o ld1x_ld1inc.o funct.o ld1x_paw_type.o radial_grids.o constants.o upf_params.o ld1x_parameters.o kind.o 
ld1x_c6_dft.o : ld1x_c6_dft.f90 radial_grids.o ld1x_ld1inc.o constants.o kind.o 
ld1x_c6_tfvw.o : ld1x_c6_tfvw.f90 funct.o radial_grids.o ld1x_ld1inc.o constants.o kind.o 
ld1x_calculate_gipaw_orbitals.o : ld1x_calculate_gipaw_orbitals.f90 radial_grids.o kind.o io_global.o ld1x_ld1inc.o ld1x_parameters.o 
ld1x_cfdsol.o : ld1x_cfdsol.f90 kind.o 
ld1x_chargeps.o : ld1x_chargeps.f90 ld1x_ld1inc.o radial_grids.o ld1x_parameters.o kind.o 
ld1x_compute_chi.o : ld1x_compute_chi.f90 ld1x_ld1inc.o radial_grids.o kind.o io_global.o 
ld1x_compute_chi_tm.o : ld1x_compute_chi_tm.f90 ld1x_ld1inc.o radial_grids.o kind.o 
ld1x_compute_phi.o : ld1x_compute_phi.f90 ld1x_ld1inc.o radial_grids.o constants.o kind.o io_global.o 
ld1x_compute_phi_tm.o : ld1x_compute_phi_tm.f90 ld1x_ld1inc.o radial_grids.o constants.o kind.o io_global.o 
ld1x_compute_phius.o : ld1x_compute_phius.f90 ld1x_ld1inc.o radial_grids.o io_global.o kind.o 
ld1x_compute_potps.o : ld1x_compute_potps.f90 ld1x_ld1inc.o radial_grids.o kind.o 
ld1x_compute_potps_new.o : ld1x_compute_potps_new.f90 ld1x_ld1inc.o radial_grids.o kind.o 
ld1x_compute_q_3bess.o : ld1x_compute_q_3bess.f90 ld1x_ld1inc.o radial_grids.o kind.o 
ld1x_compute_relpert.o : ld1x_compute_relpert.f90 ld1x_ld1inc.o radial_grids.o kind.o 
ld1x_compute_solution.o : ld1x_compute_solution.f90 radial_grids.o kind.o io_global.o 
ld1x_default_conf.o : ld1x_default_conf.f90 kind.o 
ld1x_descreening.o : ld1x_descreening.f90 ld1x_ld1inc.o ld1x_parameters.o radial_grids.o mp_world.o mp.o io_global.o kind.o 
ld1x_dfx_new.o : ld1x_dfx_new.f90 radial_grids.o ld1x_ld1inc.o ld1x_parameters.o kind.o constants.o 
ld1x_dir_outward.o : ld1x_dir_outward.f90 ld1x_ld1inc.o kind.o 
ld1x_dirsol.o : ld1x_dirsol.f90 ld1x_ld1inc.o radial_grids.o kind.o io_global.o 
ld1x_dmixp.o : ld1x_dmixp.f90 kind.o io_global.o 
ld1x_drho0ofvx.o : ld1x_drho0ofvx.f90 radial_grids.o ld1x_ld1inc.o ld1x_parameters.o kind.o constants.o 
ld1x_drhoofv.o : ld1x_drhoofv.f90 radial_grids.o ld1x_ld1inc.o kind.o constants.o 
ld1x_dvex.o : ld1x_dvex.f90 radial_grids.o ld1x_ld1inc.o constants.o kind.o 
ld1x_el_config.o : ld1x_el_config.f90 io_global.o ld1x_parameters.o kind.o 
ld1x_elsd.o : ld1x_elsd.f90 ld1x_ld1inc.o funct.o radial_grids.o constants.o kind.o 
ld1x_elsd_highv.o : ld1x_elsd_highv.f90 ld1x_ld1inc.o radial_grids.o constants.o kind.o 
ld1x_elsdps.o : ld1x_elsdps.f90 funct.o ld1x_ld1inc.o ld1x_parameters.o radial_grids.o constants.o kind.o 
ld1x_elsdps_paw.o : ld1x_elsdps_paw.f90 funct.o ld1x_ld1inc.o ld1x_parameters.o radial_grids.o constants.o kind.o 
ld1x_esic.o : ld1x_esic.f90 ld1x_ld1inc.o radial_grids.o kind.o 
ld1x_export_upf.o : ld1x_export_upf.f90 pseudo_types.o version.o funct.o ld1x_ld1inc.o kind.o constants.o 
ld1x_find_qi.o : ld1x_find_qi.f90 ld1x_ld1inc.o kind.o 
ld1x_gener_pseudo.o : ld1x_gener_pseudo.f90 invmat.o ld1x_atomic_paw.o ld1x_ld1inc.o mp.o io_global.o ld1x_parameters.o radial_grids.o kind.o 
ld1x_grad_log.o : ld1x_grad_log.f90 kind.o 
ld1x_green.o : ld1x_green.f90 radial_grids.o ld1x_ld1inc.o kind.o 
ld1x_import_upf.o : ld1x_import_upf.f90 ld1x_atomic_paw.o ld1x_parameters.o ld1x_paw_type.o pseudo_types.o funct.o ld1x_ld1inc.o radial_grids.o kind.o constants.o 
ld1x_int_0_inf_dr.o : ld1x_int_0_inf_dr.f90 radial_grids.o kind.o 
ld1x_integrate_inward.o : ld1x_integrate_inward.f90 radial_grids.o kind.o 
ld1x_integrate_outward.o : ld1x_integrate_outward.f90 radial_grids.o kind.o 
ld1x_intref.o : ld1x_intref.f90 radial_grids.o kind.o 
ld1x_inward.o : ld1x_inward.f90 radial_grids.o kind.o 
ld1x_kin_e_density.o : ld1x_kin_e_density.f90 
ld1x_kli.o : ld1x_kli.f90 radial_grids.o ld1x_ld1inc.o ld1x_parameters.o kind.o 
ld1x_ld1inc.o : ld1x_ld1inc.f90 ld1x_paw_type.o radial_grids.o ld1x_parameters.o kind.o 
ld1x_ld1_readin.o : ld1x_ld1_readin.f90 ld1x_atomic_paw.o funct.o ld1x_ld1inc.o open_close_input_file.o mp_world.o mp.o io_global.o constants.o upf_params.o ld1x_parameters.o radial_grids.o kind.o 
ld1x_ld1_setup.o : ld1x_ld1_setup.f90 funct.o ld1x_ld1inc.o kind.o 
ld1x_ld1_writeout.o : ld1x_ld1_writeout.f90 open_close_input_file.o ld1x_paw_type.o funct.o ld1x_ld1inc.o mp_world.o mp.o io_global.o radial_grids.o 
ld1x_lderiv.o : ld1x_lderiv.f90 ld1x_ld1inc.o ld1x_parameters.o mp.o io_global.o radial_grids.o kind.o 
ld1x_lderivps.o : ld1x_lderivps.f90 ld1x_ld1inc.o ld1x_parameters.o mp.o io_global.o radial_grids.o kind.o 
ld1x_lschps.o : ld1x_lschps.f90 ld1x_ld1inc.o radial_grids.o kind.o 
ld1x_newd_at.o : ld1x_newd_at.f90 ld1x_ld1inc.o radial_grids.o kind.o 
ld1x_new_potential.o : ld1x_new_potential.f90 ld1x_kli.o ld1x_ld1inc.o funct.o kind.o radial_grids.o constants.o 
ld1x_nodenum.o : ld1x_nodenum.f90 kind.o 
ld1x_normalize.o : ld1x_normalize.f90 ld1x_ld1inc.o io_global.o radial_grids.o ld1x_parameters.o kind.o 
ld1x_occ_spin.o : ld1x_occ_spin.f90 radial_grids.o kind.o 
ld1x_outward.o : ld1x_outward.f90 radial_grids.o kind.o 
ld1x_parameters.o : ld1x_parameters.f90 
ld1x_partial_wave_expansion.o : ld1x_partial_wave_expansion.f90 ld1x_ld1inc.o ld1x_parameters.o mp.o io_global.o radial_grids.o kind.o 
ld1x_paw_type.o : ld1x_paw_type.f90 radial_grids.o kind.o 
ld1x_pseudo_q.o : ld1x_pseudo_q.f90 ld1x_ld1inc.o ld1x_parameters.o io_global.o kind.o 
ld1x_pseudovloc.o : ld1x_pseudovloc.f90 ld1x_ld1inc.o io_global.o radial_grids.o kind.o 
ld1x_read_pseudo_ncpp.o : ld1x_read_pseudo_ncpp.f90 radial_grids.o constants.o kind.o 
ld1x_read_pseudo_rrkj3.o : ld1x_read_pseudo_rrkj3.f90 ld1x_ld1inc.o radial_grids.o funct.o kind.o 
ld1x_run_lda_half.o : ld1x_run_lda_half.f90 ld1x_ld1inc.o ld1x_parameters.o radial_grids.o mp.o io_global.o kind.o 
ld1x_run_pseudo.o : ld1x_run_pseudo.f90 ld1x_atomic_paw.o ld1x_ld1inc.o ld1x_parameters.o radial_grids.o kind.o 
ld1x_run_test.o : ld1x_run_test.f90 ld1x_ld1inc.o ld1x_parameters.o radial_grids.o mp_world.o mp.o io_global.o kind.o 
ld1x_scf.o : ld1x_scf.f90 ld1x_ld1inc.o constants.o radial_grids.o funct.o kind.o 
ld1x_seriebes.o : ld1x_seriebes.f90 kind.o 
ld1x_set_psi_in.o : ld1x_set_psi_in.f90 ld1x_ld1inc.o radial_grids.o kind.o 
ld1x_set_rc_rv.o : ld1x_set_rc_rv.f90 ld1x_ld1inc.o ld1x_parameters.o kind.o 
ld1x_set_rho_core.o : ld1x_set_rho_core.f90 ld1x_ld1inc.o mp_world.o mp.o io_global.o constants.o kind.o 
ld1x_set_sl3.o : ld1x_set_sl3.f90 kind.o 
ld1x_sic_correction.o : ld1x_sic_correction.f90 funct.o ld1x_ld1inc.o constants.o radial_grids.o kind.o 
ld1x_starting_potential.o : ld1x_starting_potential.f90 ld1x_ld1inc.o kind.o 
ld1x_start_potps.o : ld1x_start_potps.f90 ld1x_ld1inc.o io_global.o radial_grids.o kind.o 
ld1x_start_scheq.o : ld1x_start_scheq.f90 radial_grids.o kind.o 
ld1x_test_bessel.o : ld1x_test_bessel.f90 ld1x_ld1inc.o constants.o kind.o io_global.o 
ld1x_trou.o : ld1x_trou.f90 random_numbers.o kind.o 
ld1x_utils.o : ld1x_utils.f90 
ld1x_vdpack.o : ld1x_vdpack.f90 kind.o 
ld1x_vext.o : ld1x_vext.f90 kind.o 
ld1x_v_of_rho_at.o : ld1x_v_of_rho_at.f90 ld1x_ld1inc.o funct.o radial_grids.o constants.o kind.o 
ld1x_vpack.o : ld1x_vpack.f90 kind.o 
ld1x_vxcgc.o : ld1x_vxcgc.f90 metagga.o xc_gga_drivers.o funct.o constants.o xc_lda_lsda_drivers.o kind.o 
ld1x_write_ae_pseudo.o : ld1x_write_ae_pseudo.f90 ld1x_ld1inc.o mp.o io_global.o kind.o 
ld1x_write_cpmd.o : ld1x_write_cpmd.f90 constants.o kind.o funct.o 
ld1x_write_files.o : ld1x_write_files.f90 mp_world.o mp.o radial_grids.o ld1x_ld1inc.o io_global.o kind.o 
ld1x_write_paw_recon.o : ld1x_write_paw_recon.f90 ld1x_ld1inc.o mp_world.o mp.o io_global.o kind.o 
ld1x_write_pseudo.o : ld1x_write_pseudo.f90 funct.o constants.o kind.o 
ld1x_write_results.o : ld1x_write_results.f90 funct.o ld1x_ld1inc.o constants.o mp_world.o mp.o io_global.o kind.o radial_grids.o 
ld1x_write_resultsps.o : ld1x_write_resultsps.f90 funct.o ld1x_ld1inc.o constants.o mp.o io_global.o radial_grids.o kind.o 
ldaU.o : ldaU.f90 control_flags.o ions_base.o atomic_wfc_mod.o parameters.o upf_params.o kind.o 
linpack.o : linpack.f90 kind.o 
loc_scdm.o : loc_scdm.f90 scf_mod.o fft_base.o mp_bands.o mp.o fft_interfaces.o constants.o cell_base.o funct.o control_flags.o buffers.o io_files.o pwcom.o pwcom.o wavefunctions.o pwcom.o noncol.o exx_base.o exx.o io_global.o kind.o 
loc_scdm_k.o : loc_scdm_k.f90 mp.o loc_scdm.o exx_band.o fft_interfaces.o constants.o cell_base.o control_flags.o mp_bands.o pwcom.o pwcom.o noncol.o exx_base.o exx.o io_global.o kind.o 
make_pointlists.o : make_pointlists.f90 noncol.o fft_types.o fft_base.o mp_bands.o ions_base.o io_global.o cell_base.o kind.o 
makov_payne.o : makov_payne.f90 basic_algebra_routines.o control_flags.o mp.o mp_bands.o pwcom.o recvec.o io_files.o constants.o plugin_flags.o pwcom.o scf_mod.o fft_base.o cell_base.o ions_base.o io_global.o kind.o 
manypw.o : manypw.f90 command_line_options.o read_input.o mp_world.o mp.o mp_images.o io_global.o environment.o input_parameters.o 
martyna_tuckerman.o : martyna_tuckerman.f90 ions_base.o recvec.o gvecw.o recvec.o control_flags.o fft_types.o fft_interfaces.o fft_base.o mp.o mp_bands.o cell_base.o ws_base.o constants.o kind.o 
matches.o : matches.f90 
mdiis.o : mdiis.f90 mp.o kind.o 
memory_report.o : memory_report.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp_images.o mp_pools.o mp_bands.o ions_base.o pwcom.o control_flags.o noncol.o pwcom.o uspp.o pwcom.o pwcom.o ldaU.o funct.o atom.o uspp.o pwcom.o pwcom.o gvecw.o recvec.o recvec.o fft_base.o exx_base.o exx.o cell_base.o atomic_wfc_mod.o pwcom.o constants.o kind.o io_global.o 
metagga.o : metagga.f90 constants.o kind.o correlation_gga.o correlation_lda_lsda.o exchange_lda_lsda.o 
mix_rho.o : mix_rho.f90 fft_interfaces.o fft_base.o mp_bands.o mp.o wavefunctions.o pwcom.o cell_base.o constants.o io_files.o ldaU.o io_global.o uspp.o control_flags.o pwcom.o recvec.o recvec.o ions_base.o kind.o scf_mod.o 
mm_dispersion.o : mm_dispersion.f90 mp_images.o mp.o io_global.o constants.o cell_base.o ions_base.o parameters.o kind.o 
move_ions.o : move_ions.f90 pwcom.o fcp_variables.o fcp.o dynamics_module.o basic_algebra_routines.o bfgs_module.o mp.o mp_images.o pwcom.o pwcom.o control_flags.o pwcom.o pwcom.o symm_base.o ions_base.o pwcom.o cell_base.o kind.o io_files.o io_global.o constants.o 
mp_bands.o : mp_bands.f90 parallel_include.o mp.o 
mp_bands_util.o : mp_bands_util.f90 parallel_include.o mp.o 
mp_base.o : mp_base.f90 util_param.o parallel_include.o data_buffer.o 
mp_base_gpu.o : mp_base_gpu.f90 parallel_include.o util_param.o data_buffer.o 
mp_diag.o : mp_diag.f90 la_param.o 
mp_exx.o : mp_exx.f90 io_global.o parallel_include.o mp_bands.o mp.o 
mp.o : mp.f90 parallel_include.o util_param.o 
mp_global.o : mp_global.f90 mp.o parallel_include.o command_line_options.o mp_pools.o mp_exx.o mp_bands.o mp_bands.o mp_pools.o mp_images.o mp_world.o 
mp_images.o : mp_images.f90 parallel_include.o io_global.o mp.o 
mp_pools.o : mp_pools.f90 parallel_include.o mp.o 
mp_wave.o : mp_wave.f90 parallel_include.o kind.o 
mp_world.o : mp_world.f90 parallel_include.o io_global.o mp.o 
multable.o : multable.f90 
my_addusdens.o : my_addusdens.f90 kind.o fft_base.o noncol.o control_flags.o 
my_addusdens_g.o : my_addusdens_g.f90 mp.o mp_bands.o mp_pools.o uspp.o uspp.o noncol.o recvec.o fft_interfaces.o fft_base.o cell_base.o ions_base.o kind.o 
my_add_vuspsi.o : my_add_vuspsi.f90 becmod.o uspp.o uspp.o noncol.o control_flags.o pwcom.o ions_base.o kind.o 
my_atomic_rho.o : my_atomic_rho.f90 fft_rho.o fft_base.o mp.o mp_bands.o control_flags.o io_global.o uspp.o noncol.o pwcom.o pwcom.o recvec.o cell_base.o ions_base.o atom.o constants.o kind.o 
my_calc_pol.o : my_calc_pol.f90 bp_mod.o constants.o cell_base.o io_global.o kind.o 
my_c_bands.o : my_c_bands.f90 atomic_wfc_mod.o mp_bands.o becmod.o scf_mod.o g_psi_mod.o noncol.o check_stop.o mp.o mp_pools.o bp_mod.o wavefunctions.o pwcom.o ldaU.o control_flags.o pwcom.o recvec.o uspp.o pwcom.o buffers.o io_files.o io_global.o kind.o 
my_cegterg.o : my_cegterg.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
my_compute_magnetization.o : my_compute_magnetization.f90 mp.o pwcom.o cell_base.o pwcom.o mp_bands.o noncol.o pwcom.o fft_base.o scf_mod.o kind.o 
my_delta_e.o : my_delta_e.f90 mp.o paw_variables.o mp_bands.o cell_base.o pwcom.o fft_base.o scf_mod.o funct.o kind.o 
my_delta_escf.o : my_delta_escf.f90 mp.o paw_variables.o mp_bands.o cell_base.o pwcom.o fft_base.o scf_mod.o kind.o funct.o 
my_electrons.o : my_electrons.f90 control_flags.o funct.o io_global.o kind.o 
my_electrons_scf.o : my_electrons_scf.f90 constants.o plugin_variables.o wrappers.o fcp_variables.o esm.o my_newd.o paw_symmetry.o paw_onecenter.o paw_variables.o tsvdw.o xdm_dispersion.o dftd3_qe.o dftd3_api.o mm_dispersion.o mp.o mp_pools.o io_rho_xml.o pwcom.o noncol.o extfield.o ldaU.o io_files.o control_flags.o scf_mod.o pwcom.o gvecw.o pwcom.o pwcom.o pwcom.o pwcom.o pwcom.o recvec.o recvec.o fft_base.o bp_mod.o ions_base.o cell_base.o io_global.o check_stop.o kind.o 
my_h_psi.o : my_h_psi.f90 fft_helper_subroutines.o exx.o fft_base.o realus.o control_flags.o recvec.o ldaU.o uspp.o pwcom.o scf_mod.o pwcom.o becmod.o bp_mod.o mp.o mp_bands.o funct.o noncol.o kind.o 
my_init_run.o : my_init_run.f90 uspp.o uspp.o ions_base.o Coul_cut_2D.o tsvdw.o esm.o my_newd.o recvec_subs.o funct.o fft_base.o bp_mod.o paw_init.o paw_variables.o dynamics_module.o pwcom.o cell_base.o recvec.o recvec.o control_flags.o pwcom.o symme.o pwcom.o 
my_init_us_2.o : my_init_us_2.f90 uspp.o uspp.o splinelib.o gth.o pwcom.o pwcom.o recvec.o constants.o cell_base.o ions_base.o kind.o 
my_init_wfc.o : my_init_wfc.f90 funct.o mp.o mp_bands.o random_numbers.o wavefunctions.o noncol.o uspp.o pwcom.o pwcom.o recvec.o atomic_wfc_mod.o constants.o becmod.o bp_mod.o kind.o 
my_newd.o : my_newd.f90 ldaU.o realus.o pwcom.o uspp.o mp.o mp_pools.o mp_bands.o noncol.o wavefunctions.o control_flags.o uspp.o scf_mod.o pwcom.o recvec.o fft_interfaces.o fft_base.o cell_base.o ions_base.o kind.o 
my_potinit.o : my_potinit.f90 paw_onecenter.o paw_init.o paw_variables.o uspp.o fft_rho.o io_base.o io_rho_xml.o mp_bands.o mp.o pwcom.o io_files.o noncol.o ldaU.o pwcom.o funct.o scf_mod.o control_flags.o recvec.o recvec.o fft_base.o pwcom.o pwcom.o atomic_wfc_mod.o ions_base.o cell_base.o io_global.o constants.o kind.o 
my_qvan2.o : my_qvan2.f90 uspp.o uspp.o pwcom.o kind.o 
my_set_vrs.o : my_set_vrs.f90 fft_interfaces.o fft_base.o funct.o kind.o 
my_s_psi.o : my_s_psi.f90 fft_base.o wavefunctions.o realus.o control_flags.o ions_base.o uspp.o pwcom.o uspp.o becmod.o mp.o mp_bands.o funct.o noncol.o kind.o 
my_sum_band.o : my_sum_band.f90 fft_helper_subroutines.o becmod.o paw_variables.o paw_symmetry.o funct.o mp.o mp_bands.o mp_pools.o pwcom.o pwcom.o noncol.o wavefunctions.o uspp.o uspp.o buffers.o io_files.o symme.o scf_mod.o pwcom.o ldaU.o pwcom.o recvec.o recvec.o fft_interfaces.o fft_base.o ions_base.o cell_base.o control_flags.o pwcom.o kind.o 
my_sum_bec.o : my_sum_bec.f90 mp.o mp_bands.o us_exx.o realus.o wavefunctions.o noncol.o pwcom.o pwcom.o uspp.o uspp.o ions_base.o control_flags.o becmod.o kind.o 
my_usnldiag.o : my_usnldiag.f90 noncol.o pwcom.o uspp.o uspp.o pwcom.o pwcom.o ions_base.o kind.o 
my_v_h.o : my_v_h.f90 Coul_cut_2D.o esm.o martyna_tuckerman.o mp.o mp_bands.o control_flags.o cell_base.o pwcom.o recvec.o fft_interfaces.o fft_base.o kind.o constants.o 
my_vloc_psi_k.o : my_vloc_psi_k.f90 wavefunctions.o fft_helper_subroutines.o fft_interfaces.o fft_base.o mp_bands.o pwcom.o pwcom.o kind.o parallel_include.o 
my_v_of_rho.o : my_v_of_rho.f90 tsvdw.o control_flags.o cell_base.o scf_mod.o funct.o ldaU.o ions_base.o noncol.o recvec.o fft_base.o kind.o 
my_v_xc.o : my_v_xc.f90 mp.o mp_bands.o scf_mod.o xc_lda_lsda_drivers.o funct.o pwcom.o cell_base.o pwcom.o recvec.o fft_base.o io_global.o constants.o kind.o 
my_v_xc_meta.o : my_v_xc_meta.f90 mp_bands.o mp.o scf_mod.o xc_mgga_drivers.o funct.o cell_base.o pwcom.o recvec.o fft_base.o io_global.o constants.o kind.o 
my_wfcinit.o : my_wfcinit.f90 mp_images.o mp.o pwcom.o wavefunctions.o uspp.o buffers.o io_files.o pwcom.o ldaU.o pwcom.o control_flags.o pwcom.o bp_mod.o atomic_wfc_mod.o io_global.o 
newd.o : newd.f90 ldaU.o realus.o pwcom.o uspp.o mp.o mp_pools.o mp_bands.o noncol.o wavefunctions.o control_flags.o uspp.o scf_mod.o pwcom.o recvec.o fft_interfaces.o fft_base.o cell_base.o ions_base.o kind.o 
new_nsb.o : new_nsb.f90 becmod.o mp.o mp_pools.o buffers.o io_files.o wavefunctions.o control_flags.o pwcom.o pwcom.o symm_base.o ldaU.o pwcom.o ions_base.o kind.o io_global.o 
new_ns.o : new_ns.f90 mp_bands.o recvec.o noncol.o uspp.o uspp.o becmod.o mp.o mp_pools.o buffers.o io_files.o wavefunctions.o control_flags.o pwcom.o pwcom.o symm_base.o ldaU.o pwcom.o ions_base.o kind.o io_global.o 
new_nsg.o : new_nsg.f90 becmod.o mp.o mp_global.o buffers.o io_files.o wavefunctions.o control_flags.o pwcom.o pwcom.o symm_base.o ldaU.o pwcom.o ions_base.o kind.o io_global.o 
new_occ.o : new_occ.f90 mp.o mp_bands.o buffers.o io_files.o recvec.o noncol.o wavefunctions.o control_flags.o pwcom.o pwcom.o pwcom.o atomic_wfc_mod.o constants.o kind.o io_global.o 
noncol.o : noncol.f90 parameters.o kind.o 
non_scf.o : non_scf.f90 wavefunctions.o pwcom.o pwcom.o pwcom.o buffers.o io_files.o io_global.o pwcom.o control_flags.o check_stop.o bp_mod.o kind.o 
n_plane_waves.o : n_plane_waves.f90 mp_pools.o mp.o kind.o 
ns_adj.o : ns_adj.f90 io_global.o noncol.o pwcom.o scf_mod.o ldaU.o ions_base.o kind.o 
nsg_adj.o : nsg_adj.f90 io_global.o noncol.o pwcom.o scf_mod.o ldaU.o ions_base.o kind.o 
offset_atom_wfc.o : offset_atom_wfc.f90 ldaU.o ions_base.o noncol.o uspp.o 
open_close_input_file.o : open_close_input_file.f90 io_global.o 
openfil.o : openfil.f90 bp_mod.o noncol.o io_files.o ldaU.o pwcom.o pwcom.o atomic_wfc_mod.o control_flags.o buffers.o kind.o 
orbm_kubo.o : orbm_kubo.f90 mp_world.o mp_pools.o recvec.o scf_mod.o uspp.o becmod.o bp_mod.o constants.o mp.o cell_base.o pwcom.o start_k.o recvec.o fft_base.o pwcom.o pwcom.o noncol.o buffers.o io_files.o io_global.o kind.o pwcom.o 
orthoatwfc.o : orthoatwfc.f90 pwcom.o mp.o mp_bands.o noncol.o control_flags.o becmod.o uspp.o pwcom.o ldaU.o pwcom.o atomic_wfc_mod.o ions_base.o io_files.o io_global.o buffers.o kind.o 
ortho_wfc.o : ortho_wfc.f90 mp.o mp_bands.o io_global.o kind.o 
output_tau.o : output_tau.f90 ions_base.o cell_base.o constants.o kind.o io_global.o 
para.o : para.f90 parallel_include.o mp_images.o mp.o mp_pools.o kind.o 
parallel_include.o : parallel_include.f90 
parameters.o : parameters.f90 
paro_gamma.o : paro_gamma.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
paro_gamma_new.o : paro_gamma_new.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
paro_k.o : paro_k.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
paro_k_new.o : paro_k_new.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
parser.o : parser.f90 mp_images.o mp.o kind.o io_global.o 
paw_exx.o : paw_exx.f90 pwcom.o paw_onecenter.o atom.o constants.o io_global.o paw_variables.o uspp.o uspp.o ions_base.o kind.o 
paw_init.o : paw_init.f90 constants.o mp.o mp_images.o funct.o pwcom.o radial_grids.o atom.o noncol.o atomic_wfc_mod.o random_numbers.o paw_symmetry.o scf_mod.o uspp.o paw_variables.o uspp.o ions_base.o pwcom.o kind.o 
paw_onecenter.o : paw_onecenter.f90 io_global.o pwcom.o uspp.o radial_grids.o xc_gga_drivers.o xc_lda_lsda_drivers.o funct.o constants.o noncol.o uspp.o pwcom.o ions_base.o atom.o mp.o mp_images.o paw_variables.o kind.o 
paw_symmetry.o : paw_symmetry.f90 constants.o io_global.o uspp.o symm_base.o ions_base.o uspp.o pwcom.o noncol.o cell_base.o pwcom.o mp.o mp_images.o kind.o 
paw_variables.o : paw_variables.f90 kind.o 
pcg_gamma.o : pcg_gamma.f90 mp.o mp_bands_util.o util_param.o 
pcg_k.o : pcg_k.f90 mp.o mp_bands_util.o util_param.o 
plot_io.o : plot_io.f90 kind.o io_global.o 
plugin_check.o : plugin_check.f90 plugin_flags.o 
plugin_clean.o : plugin_clean.f90 plugin_flags.o 
plugin_clock.o : plugin_clock.f90 io_global.o plugin_flags.o 
plugin_ext_forces.o : plugin_ext_forces.f90 plugin_flags.o kind.o io_global.o mp_images.o mp.o 
plugin_flags.o : plugin_flags.f90 parameters.o kind.o 
plugin_initbase.o : plugin_initbase.f90 mp_bands.o fft_base.o plugin_flags.o 
plugin_init_cell.o : plugin_init_cell.f90 plugin_flags.o fft_base.o kind.o 
plugin_initialization.o : plugin_initialization.f90 plugin_flags.o io_files.o kind.o io_global.o 
plugin_init_ions.o : plugin_init_ions.f90 plugin_flags.o fft_base.o kind.o 
plugin_init_potential.o : plugin_init_potential.f90 fft_base.o kind.o plugin_flags.o 
plugin_int_forces.o : plugin_int_forces.f90 pwcom.o martyna_tuckerman.o fft_interfaces.o fft_base.o recvec.o cell_base.o plugin_flags.o control_flags.o pwcom.o ions_base.o io_global.o kind.o 
plugin_print_energies.o : plugin_print_energies.f90 plugin_flags.o io_files.o kind.o io_global.o 
plugin_read_input.o : plugin_read_input.f90 plugin_flags.o 
plugin_scf_energy.o : plugin_scf_energy.f90 plugin_flags.o scf_mod.o pwcom.o fft_base.o io_files.o kind.o io_global.o 
plugin_scf_potential.o : plugin_scf_potential.f90 plugin_flags.o scf_mod.o pwcom.o fft_base.o kind.o io_global.o 
plugin_summary.o : plugin_summary.f90 plugin_flags.o 
plugin_variables.o : plugin_variables.f90 parameters.o kind.o 
plus_u_full.o : plus_u_full.f90 pwcom.o noncol.o uspp.o pwcom.o pwcom.o pwcom.o recvec.o atomic_wfc_mod.o ions_base.o cell_base.o symm_base.o ldaU.o random_numbers.o invmat.o constants.o kind.o 
potinit.o : potinit.f90 paw_onecenter.o paw_init.o paw_variables.o uspp.o fft_rho.o io_base.o io_rho_xml.o mp_bands.o mp.o pwcom.o io_files.o noncol.o ldaU.o pwcom.o funct.o scf_mod.o control_flags.o recvec.o recvec.o fft_base.o pwcom.o pwcom.o atomic_wfc_mod.o ions_base.o cell_base.o io_global.o constants.o kind.o 
ppcg_gamma.o : ppcg_gamma.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp_bands_util.o mp.o util_param.o 
ppcg_k.o : ppcg_k.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp_bands_util.o mp.o util_param.o 
prepare_all.o : prepare_all.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh check_stop.o read_input.o environment.o command_line_options.o mp_bands.o mp_pools.o mp_world.o mp_global.o 
prepare_h_s_psi.o : prepare_h_s_psi.f90 buffers.o uspp.o bp_mod.o mp_bands.o pwcom.o becmod.o io_files.o wavefunctions.o pwcom.o pwcom.o 
print_clock_pw.o : print_clock_pw.f90 bp_mod.o funct.o ldaU.o noncol.o realus.o uspp.o paw_variables.o control_flags.o io_global.o 
print_ks_energies.o : print_ks_energies.f90 pwcom.o mp_pools.o mp.o mp_bands.o pwcom.o pwcom.o pwcom.o pwcom.o pwcom.o io_global.o constants.o control_flags.o kind.o 
pseudo_types.o : pseudo_types.f90 upf_kinds.o 
ptoolkit.o : ptoolkit.f90 laxlib_low.fh laxlib_param.fh laxlib_kinds.fh mp_diag.o la_param.o la_types.o 
punch.o : punch.f90 xdm_dispersion.o wavefunctions.o a2fmod.o io_rho_xml.o ions_base.o pwcom.o scf_mod.o pwcom.o wrappers.o pwcom.o control_flags.o io_files.o io_global.o 
pw2blip.o : pw2blip.f90 fft_scalar.o fft_support.o cell_base.o constants.o control_flags.o mp.o mp_pools.o io_global.o kind.o 
pw2casino.o : pw2casino.f90 plugin_flags.o io_files.o pwcom.o noncol.o mp_pools.o mp_bands.o mp_images.o kind.o 
pw2casino_write.o : pw2casino_write.f90 exx.o becmod.o pw2blip.o buffers.o mp.o mp_bands.o mp_pools.o funct.o wavefunctions.o io_files.o io_global.o uspp.o uspp.o control_flags.o gvecw.o pwcom.o pwcom.o ldaU.o scf_mod.o pwcom.o pwcom.o recvec.o fft_interfaces.o fft_base.o pwcom.o constants.o run_info.o cell_base.o ions_base.o kind.o 
pwcom.o : pwcom.f90 upf_params.o parameters.o kind.o 
qmmm.o : qmmm.f90 constraints_module.o fft_types.o constants.o ions_base.o cell_base.o input_parameters.o parallel_include.o kind.o mp.o mp_pools.o mp_world.o io_global.o 
qvan2.o : qvan2.f90 uspp.o uspp.o pwcom.o kind.o 
radial_gradients.o : radial_gradients.f90 kind.o 
radial_grids.o : radial_grids.f90 upf_const.o upf_kinds.o 
random_numbers.o : random_numbers.f90 kind.o 
rcgdiagg.o : rcgdiagg.f90 mp.o mp_bands_util.o util_param.o 
rdiagh.o : rdiagh.f90 mp.o mp_bands.o kind.o 
rdiaghg.o : rdiaghg.f90 laxlib_mid.fh laxlib_low.fh laxlib_param.fh laxlib_kinds.fh dspev_drv.o mp_diag.o la_types.o la_param.o 
read_cards.o : read_cards.f90 cell_base.o bz_form.o wrappers.o input_parameters.o parser.o wypos.o io_global.o kind.o 
read_file_new.o : read_file_new.f90 constants.o realus.o symm_base.o cell_base.o pwcom.o pwcom.o pwcom.o io_rho_xml.o scf_mod.o fft_base.o fft_rho.o recvec.o recvec_subs.o ions_base.o Coul_cut_2D.o esm.o ldaU.o funct.o newd.o paw_onecenter.o paw_init.o paw_variables.o uspp.o read_pseudo.o uspp.o kind.o pwcom.o gvecw.o recvec.o wavefunctions.o pwcom.o noncol.o pwcom.o io_files.o buffers.o control_flags.o io_global.o 
read_input.o : read_input.f90 open_close_input_file.o mp_images.o mp.o io_global.o read_cards.o read_namelists.o input_parameters.o kind.o 
read_namelists.o : read_namelists.f90 constants.o mp_images.o mp.o io_global.o input_parameters.o kind.o 
read_ncpp.o : read_ncpp.f90 upf_const.o pseudo_types.o upf_params.o upf_kinds.o 
read_pseudo.o : read_pseudo.f90 gth.o read_uspp.o upf_to_internal.o upf_auxtools.o read_upf_new.o read_upf_v1.o wrappers.o radial_grids.o funct.o pseudo_types.o io_global.o mp_images.o mp.o kind.o uspp.o uspp.o atom.o ions_base.o io_files.o 
read_upf_new.o : read_upf_new.f90 upf_utils.o pseudo_types.o upf_kinds.o xmltools.o 
read_upf_v1.o : read_upf_v1.f90 upf_utils.o pseudo_types.o upf_kinds.o 
read_uspp.o : read_uspp.f90 pseudo_types.o upf_const.o upf_invmat.o upf_io.o upf_params.o upf_kinds.o 
realus.o : realus.f90 fft_helper_subroutines.o becmod.o pwcom.o recvec.o pwcom.o wavefunctions.o fft_interfaces.o scf_mod.o noncol.o pwcom.o control_flags.o splinelib.o mp.o mp_bands.o fft_types.o atom.o uspp.o uspp.o cell_base.o ions_base.o constants.o funct.o fft_base.o io_global.o kind.o 
recips.o : recips.f90 kind.o 
recvec.o : recvec.f90 constants.o mp.o kind.o 
recvec_subs.o : recvec_subs.f90 constants.o mp.o fft_ggen.o fft_types.o kind.o 
regterg.o : regterg.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
remove_atomic_rho.o : remove_atomic_rho.f90 cell_base.o mp_bands.o mp_pools.o io_base.o scf_mod.o pwcom.o recvec.o control_flags.o io_files.o io_global.o kind.o 
remove_tot_torque.o : remove_tot_torque.f90 kind.o 
report_mag.o : report_mag.f90 pwcom.o noncol.o scf_mod.o constants.o io_global.o ions_base.o kind.o 
restart_in_electrons.o : restart_in_electrons.f90 pwcom.o pwcom.o io_files.o io_global.o kind.o 
rgen.o : rgen.f90 kind.o 
rho2zeta.o : rho2zeta.f90 kind.o io_global.o constants.o 
rotate_HSpsi_gamma.o : rotate_HSpsi_gamma.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
rotate_HSpsi_k.o : rotate_HSpsi_k.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
rotate_wfc.o : rotate_wfc.f90 control_flags.o kind.o 
rotate_wfc_gamma.o : rotate_wfc_gamma.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
rotate_wfc_k.o : rotate_wfc_k.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh mp.o mp_bands_util.o util_param.o 
run_driver.o : run_driver.f90 command_line_options.o beef.o io_files.o update_pot.o fsockets.o pwcom.o pwcom.o cell_base.o ions_base.o pwcom.o control_flags.o mp_images.o mp.o check_stop.o upf_params.o parameters.o io_global.o 
run_info.o : run_info.f90 
run_pwscf.o : run_pwscf.f90 noncol.o pwcom.o ions_base.o constants.o kind.o newd.o exx_band.o exx.o exx_base.o fft_types.o atomic_wfc_mod.o funct.o qmmm.o fft_base.o pwcom.o scf_mod.o update_pot.o mp_images.o check_stop.o pwcom.o command_line_options.o pwcom.o control_flags.o cell_base.o upf_params.o parameters.o io_global.o 
ruotaijk.o : ruotaijk.f90 kind.o 
s_1psi.o : s_1psi.f90 pwcom.o realus.o noncol.o control_flags.o becmod.o uspp.o kind.o 
save_in_cbands.o : save_in_cbands.f90 pwcom.o pwcom.o io_files.o io_global.o kind.o 
save_in_electrons.o : save_in_electrons.f90 control_flags.o pwcom.o pwcom.o io_files.o io_global.o kind.o 
scale_h.o : scale_h.f90 mp_bands.o mp.o funct.o exx.o exx_base.o start_k.o control_flags.o pwcom.o pwcom.o recvec.o constants.o pwcom.o cell_base.o io_global.o kind.o 
scatter_mod.o : scatter_mod.f90 fft_param.o fft_types.o 
scf_mod.o : scf_mod.f90 mp.o mp_bands.o paw_onecenter.o cell_base.o constants.o wavefunctions.o control_flags.o extfield.o uspp.o paw_variables.o recvec.o recvec.o fft_interfaces.o fft_base.o funct.o buffers.o ions_base.o ldaU.o pwcom.o kind.o 
set_hubbard_l.o : set_hubbard_l.f90 io_global.o 
set_hubbard_n.o : set_hubbard_n.f90 io_global.o 
set_kplusq.o : set_kplusq.f90 kind.o 
set_kup_and_kdw.o : set_kup_and_kdw.f90 kind.o 
setlocal.o : setlocal.f90 Coul_cut_2D.o qmmm.o esm.o martyna_tuckerman.o mp.o mp_bands.o control_flags.o fft_interfaces.o fft_base.o pwcom.o scf_mod.o recvec.o extfield.o cell_base.o ions_base.o constants.o kind.o io_global.o 
set_mpi_comm_4_solvers.o : set_mpi_comm_4_solvers.f90 mp_bands_util.o 
set_occupations.o : set_occupations.f90 kind.o 
set_rhoc.o : set_rhoc.f90 scf_mod.o mp.o mp_bands.o pwcom.o recvec.o fft_rho.o fft_base.o cell_base.o ions_base.o uspp.o atom.o io_global.o kind.o 
set_spin_vars.o : set_spin_vars.f90 
setup.o : setup.f90 laxlib_hi.fh laxlib_mid.fh laxlib_param.fh laxlib_kinds.fh laxlib_low.fh laxlib.fh extfield.o fcp_variables.o paw_variables.o funct.o exx_base.o exx.o noncol.o pwcom.o pwcom.o mp.o mp_bands.o mp_pools.o mp_images.o pwcom.o bp_mod.o ldaU.o uspp.o uspp.o pwcom.o control_flags.o pwcom.o symm_base.o tetra.o start_k.o electrons_base.o pwcom.o pwcom.o recvec.o gvecw.o recvec.o atomic_wfc_mod.o ions_base.o cell_base.o io_files.o io_global.o parameters.o constants.o kind.o 
set_vdw_corr.o : set_vdw_corr.f90 io_global.o 
set_vrs.o : set_vrs.f90 fft_interfaces.o fft_base.o funct.o kind.o 
simpsn.o : simpsn.f90 upf_kinds.o 
sort.o : sort.f90 kind.o 
space_group.o : space_group.f90 kind.o 
sph_bes.o : sph_bes.f90 constants.o kind.o 
sph_dbes.o : sph_dbes.f90 constants.o kind.o 
sph_ind.o : sph_ind.f90 kind.o 
spinor.o : spinor.f90 kind.o 
splinelib.o : splinelib.f90 upf_kinds.o 
s_psi.o : s_psi.f90 fft_base.o wavefunctions.o realus.o control_flags.o ions_base.o uspp.o pwcom.o uspp.o becmod.o mp.o mp_bands.o funct.o noncol.o kind.o 
start_k.o : start_k.f90 cell_base.o kind.o 
stick_base.o : stick_base.f90 fft_param.o 
stop_run.o : stop_run.f90 io_files.o environment.o mp_global.o io_global.o 
stres_cc.o : stres_cc.f90 mp.o mp_bands.o wavefunctions.o control_flags.o pwcom.o scf_mod.o pwcom.o pwcom.o recvec.o fft_interfaces.o fft_base.o cell_base.o ions_base.o uspp.o atom.o kind.o 
stres_ewa.o : stres_ewa.f90 Coul_cut_2D.o mp.o mp_bands.o constants.o kind.o 
stres_gradcorr.o : stres_gradcorr.f90 fft_rho.o fft_types.o mp.o mp_bands.o pwcom.o xc_mgga_drivers.o xc_gga_drivers.o funct.o kind.o 
stres_har.o : stres_har.f90 Coul_cut_2D.o mp.o mp_bands.o wavefunctions.o control_flags.o scf_mod.o recvec.o fft_interfaces.o fft_base.o pwcom.o cell_base.o constants.o kind.o 
stres_hub.o : stres_hub.f90 uspp.o recvec.o pwcom.o noncol.o mp_bands.o control_flags.o mp.o mp_pools.o io_global.o symme.o scf_mod.o buffers.o io_files.o atomic_wfc_mod.o pwcom.o uspp.o pwcom.o becmod.o ldaU.o pwcom.o cell_base.o ions_base.o wavefunctions.o kind.o 
stres_knl.o : stres_knl.f90 mp.o mp_bands.o mp_pools.o wavefunctions.o noncol.o control_flags.o pwcom.o symme.o buffers.o io_files.o pwcom.o gvecw.o recvec.o cell_base.o constants.o kind.o 
stres_loc.o : stres_loc.f90 Coul_cut_2D.o mp.o mp_bands.o uspp.o wavefunctions.o control_flags.o pwcom.o scf_mod.o recvec.o fft_interfaces.o fft_base.o cell_base.o ions_base.o gth.o atom.o kind.o 
stres_mgga.o : stres_mgga.f90 mp_bands.o mp_pools.o mp.o fft_base.o fft_interfaces.o pwcom.o pwcom.o io_files.o buffers.o pwcom.o funct.o wavefunctions.o scf_mod.o recvec.o cell_base.o noncol.o control_flags.o kind.o 
stres_nonloc_dft.o : stres_nonloc_dft.f90 io_global.o xc_rVV10.o xc_vdW_DF.o fft_base.o mp.o funct.o kind.o 
stress.o : stress.f90 esm.o tsvdw.o exx.o xdm_dispersion.o dftd3_qe.o dftd3_api.o mm_dispersion.o uspp.o bp_mod.o symme.o funct.o control_flags.o scf_mod.o pwcom.o ldaU.o fft_base.o recvec.o pwcom.o constants.o ions_base.o cell_base.o kind.o io_global.o 
stres_us.o : stres_us.f90 mp.o becmod.o mp_bands.o mp_pools.o noncol.o pwcom.o wavefunctions.o uspp.o uspp.o control_flags.o pwcom.o pwcom.o pwcom.o constants.o ions_base.o kind.o 
struct_fact.o : struct_fact.f90 recvec.o constants.o kind.o 
sum_band.o : sum_band.f90 us_exx.o realus.o fft_helper_subroutines.o becmod.o paw_variables.o paw_symmetry.o funct.o mp.o mp_bands.o mp_pools.o pwcom.o pwcom.o noncol.o wavefunctions.o uspp.o uspp.o buffers.o io_files.o symme.o scf_mod.o pwcom.o ldaU.o pwcom.o recvec.o recvec.o fft_interfaces.o fft_base.o ions_base.o cell_base.o control_flags.o pwcom.o kind.o 
sumkg.o : sumkg.f90 mp.o mp_pools.o kind.o 
sumkt.o : sumkt.f90 kind.o 
summary.o : summary.f90 pwcom.o pwcom.o pwcom.o symm_base.o atom.o io_files.o fcp.o fcp_variables.o exx.o realus.o martyna_tuckerman.o esm.o mp.o mp_bands.o gvecw.o pwcom.o uspp.o pwcom.o bp_mod.o funct.o pwcom.o noncol.o control_flags.o pwcom.o ldaU.o pwcom.o pwcom.o fft_base.o recvec.o recvec.o pwcom.o ions_base.o cell_base.o constants.o run_info.o kind.o io_global.o 
symm_base.o : symm_base.f90 invmat.o cell_base.o io_global.o kind.o 
symme.o : symme.f90 constants.o mp_bands.o parallel_include.o recvec.o symm_base.o cell_base.o kind.o 
symmetrize_at.o : symmetrize_at.f90 kind.o pwcom.o io_global.o 
tabd.o : tabd.f90 kind.o 
test0.o : test0.f90 fft_param.o fft_support.o fft_scalar.o fft_parallel.o fft_interfaces.o fft_types.o 
test.o : test.f90 laxlib_low.fh laxlib_hi.fh laxlib_param.fh laxlib_kinds.fh dspev_drv.o la_param.o la_types.o 
test_input_file.o : test_input_file.f90 
test_loop_addusdens_g.o : test_loop_addusdens_g.f90 mp.o mp_bands.o mp_pools.o uspp.o uspp.o noncol.o recvec.o fft_interfaces.o ions_base.o kind.o 
tetra.o : tetra.f90 constants.o pwcom.o pwcom.o pwcom.o pwcom.o io_global.o kind.o 
tg_gather.o : tg_gather.f90 fft_types.o fft_param.o 
thread_util.o : thread_util.f90 util_param.o 
transform_becsum_nc.o : transform_becsum_nc.f90 pwcom.o noncol.o pwcom.o uspp.o ions_base.o kind.o 
transform_becsum_so.o : transform_becsum_so.f90 pwcom.o noncol.o uspp.o pwcom.o uspp.o ions_base.o kind.o 
transform_qq_so.o : transform_qq_so.f90 pwcom.o uspp.o ions_base.o kind.o 
transto.o : transto.f90 laxlib_kinds.fh 
trimcheck.o : trimcheck.f90 
trnvecc.o : trnvecc.f90 kind.o 
tsvdw.o : tsvdw.f90 uspp.o parallel_include.o mp.o mp_images.o mp_bands.o kind.o ions_base.o io_global.o funct.o fft_base.o constants.o cell_base.o 
update_pot.o : update_pot.f90 constants.o becmod.o wavefunctions.o uspp.o buffers.o pwcom.o fft_rho.o io_base.o mp_bands.o mp_pools.o paw_onecenter.o paw_variables.o pwcom.o noncol.o pwcom.o extfield.o pwcom.o ldaU.o scf_mod.o pwcom.o control_flags.o fft_interfaces.o fft_base.o cell_base.o pwcom.o recvec.o mp_images.o mp.o ions_base.o io_files.o io_global.o kind.o 
upf_auxtools.o : upf_auxtools.f90 upf_io.o upf_const.o pseudo_types.o upf_kinds.o 
upf_const.o : upf_const.f90 upf_kinds.o 
upf_erf.o : upf_erf.f90 upf_kinds.o 
upf_error.o : upf_error.f90 upf_parallel_include.o 
upf_invmat.o : upf_invmat.f90 upf_kinds.o 
upf_io.o : upf_io.f90 
upf_kinds.o : upf_kinds.f90 
upf_parallel_include.o : upf_parallel_include.f90 
upf_params.o : upf_params.f90 
upf_to_internal.o : upf_to_internal.f90 upf_kinds.o radial_grids.o pseudo_types.o 
upf_utils.o : upf_utils.f90 
us_exx.o : us_exx.f90 fft_interfaces.o gvecw.o symm_base.o mp_pools.o mp.o io_global.o mp_bands.o funct.o pwcom.o noncol.o realus.o pwcom.o fft_types.o control_flags.o uspp.o constants.o recvec.o uspp.o ions_base.o cell_base.o becmod.o kind.o 
usnldiag.o : usnldiag.f90 noncol.o pwcom.o uspp.o uspp.o pwcom.o pwcom.o ions_base.o kind.o 
uspp.o : uspp.f90 upf_const.o upf_invmat.o pseudo_types.o upf_params.o upf_kinds.o 
util_param.o : util_param.f90 parallel_include.o 
utils.o : utils.f90 noncol.o pwcom.o becmod.o io_global.o kind.o 
vcsmd.o : vcsmd.f90 constraints_module.o io_files.o pwcom.o parameters.o control_flags.o pwcom.o pwcom.o dynamics_module.o pwcom.o ions_base.o cell_base.o constants.o io_global.o kind.o 
vcsubs.o : vcsubs.f90 io_global.o constants.o kind.o 
version.o : version.f90 version.h 
vhpsi.o : vhpsi.f90 mp_bands.o noncol.o mp.o control_flags.o ions_base.o scf_mod.o pwcom.o ldaU.o becmod.o kind.o 
vloc_of_g.o : vloc_of_g.f90 Coul_cut_2D.o esm.o constants.o kind.o 
vloc_psi.o : vloc_psi.f90 noncol.o pwcom.o pwcom.o pwcom.o pwcom.o fft_helper_subroutines.o wavefunctions.o fft_interfaces.o fft_base.o mp_bands.o kind.o parallel_include.o 
v_of_rho.o : v_of_rho.f90 Coul_cut_2D.o esm.o martyna_tuckerman.o fft_interfaces.o xc_lda_lsda_drivers.o pwcom.o mp_bands.o mp.o xc_mgga_drivers.o pwcom.o io_global.o constants.o tsvdw.o control_flags.o cell_base.o scf_mod.o funct.o ldaU.o ions_base.o noncol.o recvec.o fft_base.o kind.o 
volume.o : volume.f90 kind.o 
w1gauss.o : w1gauss.f90 constants.o kind.o 
wavefunctions.o : wavefunctions.f90 kind.o 
weights.o : weights.f90 io_global.o mp.o mp_pools.o mp_images.o pwcom.o pwcom.o tetra.o pwcom.o pwcom.o pwcom.o kind.o 
wfcinit.o : wfcinit.f90 funct.o mp_bands.o random_numbers.o noncol.o recvec.o constants.o becmod.o kind.o mp_images.o mp.o pwcom.o wavefunctions.o uspp.o buffers.o io_files.o pwcom.o ldaU.o pwcom.o control_flags.o pwcom.o bp_mod.o atomic_wfc_mod.o io_global.o 
wgauss.o : wgauss.f90 constants.o kind.o 
wrappers.o : wrappers.f90 io_global.o kind.o 
write_ns.o : write_ns.f90 noncol.o ldaU.o scf_mod.o io_global.o pwcom.o ions_base.o constants.o kind.o 
ws_base.o : ws_base.f90 invmat.o kind.o 
wsweight.o : wsweight.f90 kind.o 
wyckoff.o : wyckoff.f90 space_group.o kind.o 
wypos.o : wypos.f90 kind.o 
xc_gga_drivers.o : xc_gga_drivers.f90 correlation_gga.o exchange_gga.o funct.o kind.o 
xc_lda_lsda_drivers.o : xc_lda_lsda_drivers.f90 correlation_lda_lsda.o exchange_lda_lsda.o funct.o kind.o 
xc_mgga_drivers.o : xc_mgga_drivers.f90 metagga.o funct.o kind.o 
xc_rVV10.o : xc_rVV10.f90 cell_base.o recvec.o control_flags.o fft_interfaces.o fft_base.o io_global.o mp_bands.o mp.o constants.o kind.o 
xc_vdW_DF.o : xc_vdW_DF.f90 cell_base.o recvec.o correlation_lda_lsda.o control_flags.o fft_interfaces.o fft_base.o io_global.o mp_bands.o mp_images.o mp.o constants.o kind.o 
xdm_dispersion.o : xdm_dispersion.f90 funct.o paw_onecenter.o paw_variables.o io_files.o mp_bands.o mp.o mp_images.o pwcom.o cell_base.o fft_types.o fft_base.o io_global.o scf_mod.o control_flags.o splinelib.o atom.o uspp.o ions_base.o constants.o kind.o 
xmltools.o : xmltools.f90 upf_kinds.o 
ylmr2.o : ylmr2.f90 constants.o kind.o 
zhpev_drv.o : zhpev_drv.f90 laxlib_kinds.fh la_param.o 
