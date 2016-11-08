import numpy as np
import healpy as hp

import masks

#tag to add before filenames
file_tag = 'fermi-test-BS-1'

#what is the location of this global_var file? (so we can put a copy in the logs directory)
gv_file  = '/Users/bsafdi/Dropbox/pulsars_local_new/python/pulsars/global_var.py'

#routine flags and dependencies
#routines in make_maps.py
do_convert_maps             = 1
do_load_total_maps          = 1     #needs do_convert_maps (always)
#routines in diff_fermi.py
do_make_fermi_map           = 1     #(run if energy range changed)
do_load_fermi_map           = 1     #needs do_make_fermi_map (once)
do_make_fermi_plots         = 0     #needs do_load_fermi_map (always)
#routines in likelihood.py
do_make_PDF_plot            = 0     #needs do_load_galprop_map or do_load_harm_sim_map, do_load_simulation (always)
do_time_likelihood          = 0     #needs do_load_galprop_map or do_load_harm_sim_map, do_load_simulation (always)
#routines in mle.py
do_make_asimov              = 0     #needs do_load_galprop_map (always)
do_scan                     = 0     #needs do_load_galprop_map, and do_load_simulation or do_make_asimov (always)
#routines in tri.py
do_make_chain_plots         = 0
do_make_flux_fraction_plot  = 0     #needs do_make_chain_plots, do_load_galprop_map, do_load_simulation (always)
do_make_iso_fraction_plot   = 0     #needs do_make_chain_plots, do_load_galprop_map, do_load_simulation (always)
do_make_triangle_plot       = 0     #needs do_make_chain_plots (always)
do_make_walker_plots        = 0     #needs do_make_chain_plots (always)

#plot mode for simulation skymap plots
skymap_plot_mode            = 1     #0 = full sky mollview, 1 = inner galaxy cartview               

#flag to select diffuse model to use for simulation
which_diff_model_sim        = 0     #0 = Fermi, 1 = Galprop
#flag to select diffuse model to use for PDF
which_diff_model_PDF        = 0     #0 = Fermi, 1 = Galprop, 2 = harmonic simulated, 3 = wavelet-denoised simulated

#specify directory for logs
logs_dir  = 'logs/'

#defs for hists.py
#specify directory for plots
plots_dir = 'plots/' + file_tag + '/'
sim_plots_dir = 'plots/' + file_tag + '/sim_plots/'

#defs for make_maps.py
#healpy resolution parameters
NSIDE   = 128
NPIX    = hp.nside2npix(NSIDE)
#specify directories for original CTBCORE and converted maps
CTBCORE_maps_dir    = '/home/slee/Storage/CTBCORE-maps/' #change directory as needed
#CTBCORE_maps_dir    = '/group/hepheno/samuelkl/CTBCORE-maps/' #change directory as needed
data_dir            = CTBCORE_maps_dir + 'allsky/p10_ultraclean_Q2/specbin/' 
converted_dir       = CTBCORE_maps_dir + 'converted_maps/'
#energy bins for original CTBCORE filenames
en_bins_str     = ['000.3', '000.4', '000.5', '000.6', '000.8', '000.9', '001.2', '001.5', \
                   '001.9', '002.4', '003.0', '003.8', '004.8', '006.0', '007.5', '009.5', '011.9']
num_en_bins     = len(en_bins_str) - 1               
#initialize array for exact energy bins (retrieved from header files)
en_bins         = [0.3]
#initialize array for mean exposures in each energy bin (calculated during map conversion)
mean_exposure   = []                        #mean exposure in each en_bin, averaged over sky
avg_exposure    = 4.E10                     #overall mean exposure (averaged over sky and energy)
                                            #put in by hand for now, use for scaling Fermi diffuse and IGPS fluxes to counts
#initialize total maps
total_map               = np.zeros(NPIX)    #total map (summed over energy bins) with point-source removal
total_map_nopsc         = np.zeros(NPIX)    #total map with no point-source removal
total_map_smth          = np.zeros(NPIX)    #total smoothed map with point-source removal
total_map_nopsc_smth    = np.zeros(NPIX)    #total smoothed map with no point source removal

#defs for masks.py
#set mask for all maps
#mask_for_all_maps = masks.mask_lat_band(60, 120, NSIDE)
#mask_for_all_maps = masks.mask_not_lat_band(82.5, 97.5, NSIDE)
#mask_for_all_maps = np.logical_not(masks.mask_ring(1, 10, 90, 0, NSIDE)) + masks.mask_lat_band(87.5, 92.5, NSIDE)*np.logical_not(masks.mask_ring(0, 3, 90, 0, NSIDE))
mask_for_all_maps = np.logical_not(masks.mask_ring(5, 10, 90, 0, NSIDE)) + masks.mask_lat_band(87.5, 92.5, NSIDE)
#mask_for_all_maps = np.zeros(NPIX)
npixROI           = np.sum(np.logical_not(mask_for_all_maps))

#initialize diffuse model maps
diff_model_map          = np.zeros(NPIX)
diff_model_map_PSF      = np.zeros(NPIX)
diff_model_min          = 0.0
diff_model_min_PSF      = 0.0
diff_model_map_mean_before_zeroing = 0.0

#defs for diff_galprop.py
#specify directory and filenames of Galprop output
galprop_dir         = CTBCORE_maps_dir + 'results_54_05620004/'
ics_file            = 'ics_isotropic_healpix_54_05620004'
brems_file          = 'bremss_healpix_54_05620004'
pi0_file            = 'pi0_decay_healpix_54_05620004'
galprop_min_en_bin  = 8     #integrate galprop output above given bin, 0 = 0.3 GeV, 5 = 0.95 GeV, 8 = 1.89 GeV
#initialize Galprop diffuse-model maps
gmt         = np.zeros(NPIX)
gmt_um      = np.zeros(NPIX)
gmt_um_mean = 0.0
gmt_mean    = 0.0
gmt_min     = 0.0

#defs for diff_fermi.py
#specify directory and filenames of Fermi healpix maps
fermi_dir       = CTBCORE_maps_dir + 'fermi_healpix/'
fermi_en_min    = 11        #5 = 0.28 GeV, 9 = 0.98 GeV, 11 = 1.83 GeV
fermi_en_max    = 17        #17 = 11.98 GeV
fermi_en        = np.zeros(fermi_en_max - fermi_en_min)
#initialize Fermi diffuse-model maps
fmt_um          = np.zeros(NPIX)
fmt_um_mean     = 0.0
fmt_um_min      = 0.0
fmt             = np.zeros(NPIX)
fmt_mean        = 0.0
fmt_min         = 0.0
fermi_mean_exp  = 0.0
#for galprop 0.95-11.9 GeV use ~ 11., 1.83-11.9 GeV use ~3. for nside=128                  

theta_true = [A_true, n1_true, n2_true, Sb_true, x_iso_true]

#defs for simulate_nfw.py
r_s             = 20.   #NFW scale radius in kpc
gamma_nfw       = 1.26  #gNFW power, use values from Daylan et al. fits
sigma_v         = 1.7   #in 10^-26 cm^3 / s
rho_0           = 0.3   #in GeV / cm^3
m_chi           = 35.   #in GeV
photons_per_ann = 1.67  #1.67 for 1.9-11.9 GeV, 3.63 for 1.-11.9 GeV

#defs for diff_harm_sim.py
#choose point-source threshold and multipole truncation for harmonic diffuse-model maps
point_source_threshold = 1000.0
lmax = 256
#initialize harmonic diffuse-model maps
harm_sim_map        = np.zeros(NPIX)
harm_sim_map_um     = np.zeros(NPIX)
hmt_um              = np.zeros(NPIX)

#defs for diff_wave_sim.py
#choose parameters for wavelet-denoised diffuse-model maps
wave_l_c = 256
wave_j_max = 4
wave_kappa = 0.
wave_n_max = 30
#initialize wavelet-denoised diffuse-model maps
wave_sim_map        = np.zeros(NPIX)
wave_sim_map_um     = np.zeros(NPIX)

#defs for special.py
#choose cutoff for pre-calculated factorial arrays
sp_max      = 2000

#defs for likelihood.pyx
f_ary               = [1.]
df_rho_div_f_ary    = [1.]
l10A_min            = 1.
l10A_max            = 6.
n1_min              = 5.
n1_max              = 10.
n2_min              = 0.5
n2_max              = 1.75
Sb_min              = 3.
Sb_max              = 50.

#defs for mle.py
test_prior_ranges = 1   #use flat likelihood function to ensure prior ranges cover all PS flux fractions (should be 0 for regular runs)
nsims             = 150
#specify directory for scan results
chains_dir = 'chains/'
#choose parameters for scan
PDF_bin     = 1.
k_max       = 70       #cutoff for photon counts (k_max < sp_max) -- also used for making simulation and PDF plots
ndim        = 5         #number of parameters (fixed to 5)
#nwalkers    = 200       #number of walkers
#iterations  = 2000      #number of iterations
#threads     = 15        #number of threads -- 15 for hepheno, change for local
nwalkers    = 1000       #number of walkers
iterations  = 500       #number of iterations
threads     = 4         #number of threads -- 15 for hepheno, change for local

asimov_counts_hist = np.zeros(k_max)

#defs for tri.py
#specify parameter-name text for plot labels
param_names = ["$log_{10} A$", "$n_1$", "$n_2$", "$S_b$", "$x_{iso}$"]
nburnin     = 100       #number of iterations to discard for burn-in
