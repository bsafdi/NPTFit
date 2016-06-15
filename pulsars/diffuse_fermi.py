import os
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits

import masks

import logging
logger = logging.getLogger(__name__)

from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
from scipy import integrate

def make_fermi_healpix(diff_dir, nside, fermi_en_min, fermi_en_max, nside_up_factor=2, model='gll_iem_v02_P6_V11_DIFFUSE'):
    '''Input diff_dir (where Fermi diffuse-model maps/files are located), nside, and indices of min & max Fermi energies.
    Converts Fermi diffuse-model intensity maps (in GLAT-GLON) to healpix and writes to diff_dir/healpix/model-en.hpx.
    Cubic interpolation over GLAT-GLON grid sampled at nside_up and downgraded to nside.
    Healpix intensity maps given in units 1/cm^2/s/sr/MeV.'''
    diff_healpix_dir = diff_dir + model + '/healpix/' 
    if not os.path.exists(diff_healpix_dir):
        os.mkdir(diff_healpix_dir)
    logging.info('Loading Fermi diffuse model...')
    fermi_fits = fits.open(diff_dir + model + '/' + model + '.fit')
    lng_lat_maps = fermi_fits[0].data
    dlng = 360./np.shape(lng_lat_maps)[2]
    dlat = 180./np.shape(lng_lat_maps)[1]
    lng_array = np.arange(-180, 180, dlng) + dlng/2
    lat_array = np.arange(-90, 90, dlat) + dlat/2
    nside_up = nside_up_factor*nside
    npix_up = 12*nside_up**2
    hp_pix_in_lat_lng_up = np.rad2deg(hp.pix2ang(nside_up, range(npix_up)))
    fermi_fits.close()
    for en in range(fermi_en_min, fermi_en_max+1):
        logging.info('%s', 'Converting Fermi diffuse model at Fermi energy index ' + str(en))
        diff_interp = RectBivariateSpline(lng_array, lat_array, lng_lat_maps[en].T, kx=1, ky=1)
        # diff_healpix_128 = diff_interp(mod(hp_pix_in_lat_lng[1] - 180, 360) - 180, mod(-hp_pix_in_lat_lng[0], 180) - 90, grid=False)
        diff_healpix_up = diff_interp(np.mod(hp_pix_in_lat_lng_up[1] - 180, 360) - 180, np.mod(-hp_pix_in_lat_lng_up[0], 180) - 90, grid=False)
        diff_healpix = hp.ud_grade(diff_healpix_up, nside)
        diff_healpix_filename = diff_healpix_dir + model + '-' + str(en) + '.hpx'
        hp.write_map(diff_healpix_filename, diff_healpix)
    logging.info('%s', 'Fermi diffuse model converted to healpix maps and saved to ' + diff_healpix_dir)

def get_fermi_en_int(diff_dir, nside, fermi_en_min, fermi_en_max, model='gll_iem_v02_P6_V11_DIFFUSE'):
    '''Input diff_dir (where Fermi diffuse-model maps/files are located), nside, and indices of min & max Fermi energies.
    Return arrays fermi_en and fermi_intensity_maps.
    fermi_en gives energies (in GeV).
    fermi_intensity_maps gives intensity maps at those energies in units of 1/cm^2/s/sr/GeV.'''
    logging.info('Loading Fermi healpix maps...')
    diff_healpix_dir = diff_dir + model + '/healpix/' 
    fermi_fits = fits.open(diff_dir + model + '/' + model + '.fit')
    fermi_en = np.array([0.001 * en for en_list in fermi_fits[1].data[fermi_en_min:fermi_en_max+1] for en in en_list])
    fermi_fits.close()
    #load/convert healpix Fermi diffuse-model intensity maps
    npix = 12*nside**2
    fermi_intensity_maps = np.zeros((fermi_en_max - fermi_en_min + 1, npix))
    for en in range(fermi_en_min, fermi_en_max+1):
        input_map_filename = diff_healpix_dir + model + '-' + str(en) + '.hpx'
        if not os.path.exists(input_map_filename):
            logging.info('Fermi healpix maps not created yet.  Converting...')
            make_fermi_healpix(diff_dir, nside, fermi_en_min, fermi_en_max, model=model)
        input_map = hp.read_map(input_map_filename, verbose=False)
        fermi_intensity_maps[en - fermi_en_min] = hp.ud_grade(1000.*input_map, nside)
    logging.info('Fermi healpix maps loaded.')
    return fermi_en, fermi_intensity_maps
    
def rebin_fermi_int(intensity_maps, fermi_en, en_bins):    
    '''Input healpix maps of differential intensity from Fermi diffuse model, 
    corresponding energies fermi_en, and (CTBCORE) energy-bin edges en_bins.
    Return two arrays of binned intensity maps (integrated over bins given by en_bins and fermi_en, respectively), 
    using linear interpolation in log space.'''
    if min(fermi_en) > min(en_bins):
        logging.info('Error: Lowest edge of en_bins out of range of fermi_en!  Adjust fermi_en_min.')
    elif max(en_bins) > max(fermi_en):
        logging.info('Error: Highest edge of en_bins out of range of fermi_en!  Adjust fermi_en_max.')
    else:
        #assume intensity = A * (fermi_en/GeV)^-alpha in each energy bin
        #power-law slope in energy bins defined by en_bins
        alpha = [-np.log(intensity_maps[en+1]/intensity_maps[en])/np.log(fermi_en[en+1]/fermi_en[en]) 
                 for en in range(len(fermi_en)-1)]
        #power-law normalization in energy bins defined by en_bins
        A = [intensity_maps[en]*np.power(fermi_en[en], alpha[en]) for en in range(len(fermi_en)-1)]
        #calculate binned intensity maps in energy bins defined by fermi_en (will use below)
        fermi_en_binned_intensity_maps = np.array([A[en]/(alpha[en] - 1.) * (np.power(fermi_en[en], 1. - alpha[en]) - 
                                                                             np.power(fermi_en[en+1], 1. - alpha[en]))
                                                   for en in range(len(fermi_en)-1)])
        #make binned intensity maps in energy bins defined by en_bins
        if len(np.shape(intensity_maps)) > 1:
            en_bins_binned_intensity_maps = np.zeros((len(en_bins) - 1, len(intensity_maps[0])))
        else:
            en_bins_binned_intensity_maps = np.zeros(len(en_bins) - 1)
#        print 'Energy-bin edges:', en_bins
#        print 'Fermi energies:', fermi_en
        for en_bin in range(len(en_bins) - 1):
            en_lo = en_bins[en_bin]
            en_hi = en_bins[en_bin + 1]
#            print ''
#            print 'Calculating binned intensity map for bin:', en_lo, en_hi
            #get indices of range in fermi_en that completely contains [en_lo, en_hi]
            fermi_en_lo_idx = max(np.sum(en_lo >= fermi_en) - 1, 0)
            fermi_en_hi_idx = min(np.sum(fermi_en <= en_hi), len(fermi_en))
            #if [en_lo, en_hi] completely contained in single fermi_en bin, integrate
            if fermi_en_hi_idx == fermi_en_lo_idx + 1:
                A_lo = A[fermi_en_lo_idx]
                alpha_lo = alpha[fermi_en_lo_idx]
                en_bins_binned_intensity_maps[en_bin] += A_lo/(alpha_lo - 1.) * (np.power(en_lo, 1. - alpha_lo) - 
                                                                                 np.power(en_hi, 1. - alpha_lo))
#                print 'Integrated intensity map in range:', en_lo, en_hi  
            else:  
                fermi_en_lo = fermi_en[fermi_en_lo_idx + 1]
                fermi_en_hi = fermi_en[fermi_en_hi_idx - 1]
                #add intensity map integrated over en_lo to fermi_en_lo
                if en_lo < fermi_en[-2]:
                    A_lo = A[fermi_en_lo_idx]
                    alpha_lo = alpha[fermi_en_lo_idx]
                    en_bins_binned_intensity_maps[en_bin] += A_lo/(alpha_lo - 1.) * (np.power(en_lo, 1. - alpha_lo) - 
                                                                                     np.power(fermi_en_lo, 1. - alpha_lo))
#                    print 'Added binned intensity map in range:', en_lo, fermi_en_lo
                #add binned intensity maps of fermi_en bins completely contained in [en_lo, en_hi]
                if fermi_en_lo_idx + 2 < fermi_en_hi_idx:
                    en_bins_binned_intensity_maps[en_bin] += np.sum(fermi_en_binned_intensity_maps[fermi_en_lo_idx + 1:fermi_en_hi_idx - 1],
                                                                    axis=0)
#                    print 'Added binned intensity maps in range:', fermi_en[fermi_en_lo_idx + 1], fermi_en[fermi_en_hi_idx - 1]
                #add intensity map integrated over fermi_en_hi to en_hi
                A_hi = A[fermi_en_hi_idx - 1]
                alpha_hi = alpha[fermi_en_hi_idx - 1]
                en_bins_binned_intensity_maps[en_bin] += A_hi/(alpha_hi - 1.) * (np.power(fermi_en_hi, 1. - alpha_hi) - 
                                                                                 np.power(en_hi, 1. - alpha_hi))
#                print 'Added binned intensity map in range:', fermi_en_hi, en_hi
        logging.info('Fermi diffuse model rebinned in energy.')
        return en_bins_binned_intensity_maps, fermi_en_binned_intensity_maps
        
def get_fermi_iso(diff_dir, model='gll_iem_v02_P6_V11_DIFFUSE'):    
    '''Input diff_dir (where Fermi diffuse-model maps/files are located) and indicate diffuse model.
    Return fermi_en_iso (energies in Fermi isotropic model) and isotropic differential intensity.'''
    if model == 'gll_iem_v02_P6_V11_DIFFUSE':
        fermi_iso_file = diff_dir + model +'/' + 'isotropic_iem_front_v02_P6_V11_DIFFUSE.txt'
    elif model == 'gll_iem_v05_rev1':
        fermi_iso_file = diff_dir + model +'/' + 'iso_clean_front_v05.txt'
    else:
        logging.info('Fermi diffuse model ' + model + ' not known!')
        return
    fermi_iso_en, fermi_iso_int = np.loadtxt(fermi_iso_file).T[:2]
    fermi_iso_en = 0.001 * fermi_iso_en
    fermi_iso_int = 1000. * fermi_iso_int
    return fermi_iso_en, fermi_iso_int
    
def rebin_fermi_iso(diff_dir, en_bins, model='gll_iem_v02_P6_V11_DIFFUSE'):    
    '''Input diff_dir (where Fermi diffuse-model maps/files are located), (CTBCORE) energy-bin edges en_bins and indicate diffuse model.
    Return fermi_en_iso (energies in Fermi isotropic model) and two arrays of binned intensity maps 
    (integrated over bins given by en_bins and fermi_iso_en, respectively),
    using linear interpolation in log space.'''
    if model == 'gll_iem_v02_P6_V11_DIFFUSE':
        fermi_iso_file = diff_dir + model +'/' + 'isotropic_iem_front_v02_P6_V11_DIFFUSE.txt'
    elif model == 'gll_iem_v05_rev1':
        fermi_iso_file = diff_dir + model +'/' + 'iso_clean_front_v05.txt'
    else:
        logging.info('Fermi diffuse model ' + model + ' not known!')
        return
    fermi_iso_en, fermi_iso_int = np.loadtxt(fermi_iso_file).T[:2]
    fermi_iso_en = 0.001 * fermi_iso_en
    fermi_iso_int = 1000. * fermi_iso_int
    en_bins_binned_intensity, fermi_iso_en_binned_intensity = rebin_fermi_int(fermi_iso_int, fermi_iso_en, en_bins)
    return en_bins_binned_intensity, fermi_iso_en_binned_intensity

def integrate_spectrum_fermi(spectrum, fermi_en):    
    '''Input healpix maps of differential intensity from Fermi diffuse model and corresponding energies fermi_en.
    Return binned intensity maps (integrated over energies bins defined by fermi_en, using linear interpolation in log space).
    Can also apply to intensity maps.'''
    alpha = [-np.log(spectrum[en+1]/spectrum[en])/np.log(fermi_en[en+1]/fermi_en[en]) for en in range(len(fermi_en)-1)]
    A = [spectrum[en]*np.power(fermi_en[en], alpha[en]) for en in range(len(fermi_en)-1)]
    return np.array([A[en]/(1. - alpha[en])*(np.power(fermi_en[en+1], 1.-alpha[en]) - np.power(fermi_en[en], 1.-alpha[en]))
                    for en in range(len(fermi_en)-1)])             

def load_fermi_map_tot(diff_dir, sim_dir, fermi_en_min, fermi_en_max, nside,
                       avg_exposure, mask_for_all_maps, model='gll_iem_v02_P6_V11_DIFFUSE', do_make_fermi_map=1):
    #create or load total Fermi diffuse-model map
    #load Fermi energies from model.fit
    fermi_fits = fits.open(diff_dir + model + '/' + model + '.fit')
    fermi_en = np.array([0.001 * en for en_list in fermi_fits[1].data[fermi_en_min:fermi_en_max+1] for en in en_list])
    fermi_en_range_str = '{0:.2f}-{1:.2f}'.format(fermi_en[0], fermi_en[-1])
    fermi_fits.close()
    #create total Fermi diffuse-model map and save to file
    if do_make_fermi_map:
        logging.info('Creating Fermi diffuse model map...')
        fermi_healpix_maps = get_fermi_en_int(diff_dir, nside, fermi_en_min, fermi_en_max, model)[1]
        fermi_healpix_maps_integrated = integrate_spectrum_fermi(fermi_healpix_maps, fermi_en).transpose()
        xbg_um = np.sum(fermi_healpix_maps_integrated * avg_exposure * 4.*np.pi/len(fermi_healpix_maps[0]), axis=1)
        hp.write_map(sim_dir + 'fermi_map_tot-' + fermi_en_range_str + '.fits', xbg_um)
        logging.info('%s', 'Fermi diffuse model map saved to ' + sim_dir + 'fermi_map_tot-' +  fermi_en_range_str + '.fits')
    #load total Fermi diffuse-model map from file
    else:
        xbg_um = hp.read_map(sim_dir + 'fermi_map_tot-' + fermi_en_range_str + '.fits', verbose=False)
    x_iso_true = np.mean(xbg_um)
    xbg = hp.ma(xbg_um - x_iso_true)
    xbg.mask = mask_for_all_maps
    logger.info('Fermi diffuse model loaded.')
    logger.info('Fermi diffuse model has overall mean of {0:.2f} photons per pixel from {1:.2f}-{2:.2f} GeV.'.format(np.mean(xbg_um),
                                                                                                                     fermi_en[0]/1000.,
                                                                                                                     fermi_en[-1]/1000.))
    logger.info('In the unmasked region containing {0:.2f} pixels, Fermi diffuse model has mean of {1:.2f} photons per pixel from {2:.2f}-{3:.2f} GeV.'.format(np.sum(np.logical_not(mask_for_all_maps)), np.mean(xbg.compressed()) + x_iso_true, fermi_en[0], fermi_en[-1]))
    return fermi_en, xbg_um, xbg, x_iso_true

    
#def plot_fermi_map_tot():
#    fermi_en_range_str = '{0:.2f}-{1:.2f} GeV'.format(gv.fermi_en[0]/1000., gv.fermi_en[-1]/1000.)
#    #plot histogram of total Fermi diffuse-model map and save to file
#    hists.plot_hist(xbg.compressed(), title='Fermi diffuse model: counts histogram (integrated over ' + fermi_en_range_str + ')',
#                    plot_file='fermi-fmt-hist.png',
#                    xmin=0, xmax=gv.k_max)
#    plt.clf()
#    #plot skymap of total Fermi diffuse-model map and save to file
#    hp.mollview(xbg, norm='hist', title='Fermi diffuse skymap (integrated over ' + fermi_en_range_str + ')')
#    plt.savefig(gv.plots_dir + 'fermi-fmt-map.png')
#    logging.info('%s', 'Saved Fermi diffuse model histogram and skymap to ' + gv.plots_dir + 'fermi-fmt-hist.png and ' \
#            + gv.plots_dir + 'fermi-fmt-map.png')
#    #plt.show()
#    plt.close()
