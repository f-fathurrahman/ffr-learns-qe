# Makefile for phonon (PH)

include ../qe-6.6/make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../qe-6.6/PW/src \
          $(MOD_FLAG)../qe-6.6/dft-d3 \
          $(MOD_FLAG)../qe-6.6/LR_Modules
PHOBJS = \
acfdtest.o \
add_dkmds.o \
add_for_charges.o \
add_zstar_ue.o \
add_zstar_ue_us.o \
addcore.o \
adddvscf_ph_mag.o \
addnlcc.o \
addnlcc_zstar_eu_us.o \
addusddens.o \
addusddense.o \
addusdynmat.o \
adddvhubscf.o \
allocate_part.o \
allocate_pert.o \
allocate_phq.o \
apply_trev.o \
yambo.o \
ahc.o \
dvscf_interpolate.o \
bcast_ph_input.o \
c_bands_ph.o \
cch_psi_all.o \
ccg_psi.o \
check_if_partial_dyn.o \
check_initial_status.o \
check_restart_recover.o \
clean_pw_ph.o \
clinear.o \
close_phq.o \
compute_alphasum.o \
compute_becalp.o \
compute_becsum_ph.o \
compute_drhous.o \
compute_drhous_nc.o \
compute_dvloc.o \
compute_nldyn.o \
compute_weight.o \
deallocate_part.o \
deallocate_phq.o \
d2ionq.o \
d2ionq_disp.o \
d2nsq_bare.o \
dnsq_bare.o \
dnsq_orth.o \
dnsq_scf.o \
davcio_drho.o \
dfile_autoname.o \
dfile_star.o \
dielec.o \
do_phonon.o \
do_q2r.o \
drho.o \
drhodv.o \
drhodvloc.o \
drhodvnl.o \
drhodvus.o \
dvanqq.o \
dvkb3.o \
dvpsi_e.o \
dvqpsi_us.o \
dvqpsi_us_only.o \
dyndia.o \
dynmat0.o \
dynmat_us.o \
dynmatcc.o \
dynmatrix.o \
dynmat_hub_bare.o \
dynmat_hub_scf.o \
dwfc.o \
dvqhub_barepsi_us.o \
dvqhub_barepsi_us2.o \
delta_sphi.o \
dpsi_orth.o \
ef_shift.o \
ep_matrix_element_wannier.o \
elph.o \
el_ph_collect.o \
elphon.o \
elph_scdft_mod.o \
elph_tetra_mod.o \
find_equiv_sites.o \
find_mode_sym.o \
find_irrep.o \
generate_dynamical_matrix_c.o \
generate_effective_charges_c.o \
gmressolve_all.o \
incdrhous.o \
incdrhous_nc.o \
initialize_ph.o \
init_representations.o \
io_dyn_mat.o \
io_dyn_mat_old.o \
io_pattern.o \
non_scf_ph.o \
obsolete.o \
openfilq.o \
phcom.o \
ph_restart.o \
phescf.o \
phq_init.o \
phq_readin.o \
phq_recover.o \
phq_setup.o \
phq_summary.o \
phqscf.o \
polariz.o \
print_clock_ph.o \
prepare_q.o \
prepare_sym_analysis.o \
psidspsi.o \
psymdvscf.o \
psyme.o \
psym_dmag.o \
psym_dmage.o \
punch_plot_e.o \
q_points.o \
q_points_wannier.o \
q2qstar_ph.o \
random_matrix.o \
read_wfc_rspace_and_fwfft.o \
rigid.o \
rotate_dvscf_star.o \
rotate_and_add_dyn.o \
rotate_pattern_add.o \
run_nscf.o \
save_ph_input.o \
set_asr_c.o \
set_drhoc.o \
set_int12_nc.o \
set_irr.o \
set_irr_nosym.o \
set_irr_sym.o \
setlocq.o \
solve_e.o \
solve_e_fpol.o \
solve_linter.o \
stop_ph.o \
summarize.o \
sym_and_write_zue.o \
sym_def.o \
sym_dmag.o \
sym_dmage.o \
symdvscf.o \
symdyn_munu.o \
symdynph_gq.o \
syme.o \
symm.o \
symmorphic_or_nzb.o \
swfc.o \
syme_dns.o \
sym_dns.o \
sym_dns_wrapper.o \
transform_int_so.o \
transform_int_nc.o \
transform_alphasum_nc.o \
transform_alphasum_so.o \
tra_write_matrix.o \
trntnsc.o \
write_epsilon_and_zeu.o \
write_qplot_data.o \
write_matrix.o \
write_modes.o \
write_rec.o \
write_hub.o \
zstar_eu.o \
zstar_eu_us.o \
my_phqscf.o my_setup_nscf.o my_cgsolve_all.o

# Put my_* files at the end (no need to edit make.depend)


RAMANOBJS =    \
psyme2.o   \
chi_test.o \
dielec_test.o \
dvpsi_e2.o \
el_opt.o   \
hdiag.o    \
pcgreen.o  \
ramanm.o   \
syme2.o    \
solve_e2.o \
solve_e_nscf.o \
write_ramtns.o \
d2mxc.o    \
raman_mat.o\
raman.o    \
dhdrhopsi.o\
set_dvscf.o

PHAUXOBJS = \
write_eigenvectors.o

LRMODS = ../qe-6.6/LR_Modules/liblrmod.a
PWOBJS = ../qe-6.6/PW/src/libpw.a
QEMODS = ../qe-6.6/Modules/libqemod.a ../qe-6.6/KS_Solvers/libks_solvers.a ../qe-6.6/upflib/libupf.a \
         ../qe-6.6/FFTXlib/libqefft.a ../qe-6.6/LAXlib/libqela.a ../qe-6.6/UtilXlib/libutil.a ../qe-6.6/dft-d3/libdftd3qe.a

all : libs-ph ph.x dynmat.x matdyn.x q2r.x lambda.x fqha.x q2qstar.x \
	alpha2f.x epa.x dvscf_q2r.x postahc.x

libs-ph : libph.a libphaux.a

ph.x : phonon.o libph.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ phonon.o libph.a \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

dynmat.x : dynmat.o libphaux.a  libph.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ dynmat.o libphaux.a libph.a \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

matdyn.x : matdyn.o libphaux.a libph.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ matdyn.o libphaux.a libph.a \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

q2r.x : q2r.o libphaux.a libph.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ q2r.o libphaux.a libph.a \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)


q2qstar.x : q2qstar.o libph.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ q2qstar.o libph.a \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

lambda.x : lambda.o $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ lambda.o \
                 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

alpha2f.x : alpha2f.o libph.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ alpha2f.o libph.a \
                 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

epa.x : epa.o libph.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ epa.o libph.a \
                 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

dvscf_q2r.x : dvscf_q2r.o libphaux.a libph.a $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ dvscf_q2r.o libphaux.a libph.a \
         $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

fqha.x : fqha.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ fqha.o \
                 $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

postahc.x : postahc.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ postahc.o \
                 $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

libph.a : $(PHOBJS) $(RAMANOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

libphaux.a : $(PHAUXOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ph.x dynmat.x matdyn.x \
	q2r.x q2trans.x q2trans_fd.x lambda.x fqha.x \
	alpha2f.x disp_width.x

include make.depend
# DO NOT DELETE
