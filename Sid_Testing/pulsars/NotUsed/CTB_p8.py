import os
import numpy as np
import healpy as hp
from astropy.io import fits

import logging
logger = logging.getLogger(__name__)

#energy bins for original CTBCORE filenames
en_bins_str     = ['000.3', '000.4', '000.5', '000.6', '000.8', '000.9', '001.2', '001.5', \
                   '001.9', '002.4', '003.0', '003.8', '004.8', '006.0', '007.5', '009.5', '011.9', '015.0','018.9','023.8','030.0','037.8','047.5','059.9','075.4','094.9','119.4','150.4','189.3','238.3','300.0']
num_en_bins     = len(en_bins_str) - 1               

def read_CTB_count_map(input_filename, nside):
    #convert original CTBCORE intensity (counts/cm^2/s) to counts map at nside
    #open original .fits file
    f = fits.open(input_filename)
    #[0] field contains intensity, [1] field contains exposure -- use to get counts and convert to nside
    count_map = hp.ud_grade(np.around(f[0].data*f[1].data*4*np.pi/float(len(f[0].data))), nside, power=-2)
    f.close()
    return count_map
    
def fits_filename(CTB_dir, en_bin, nopsc, smoothed):
    #function to help build .fits filenames depending on energy bin, point-source subtraction, and smoothing
    if smoothed:
        fwhm = '120'
    else:
        fwhm = '000'
    if nopsc:
        fn = 'fermi-allsky-' + en_bins_str[en_bin] + '-' + en_bins_str[en_bin+1] \
                        + 'GeV-fwhm' + fwhm + '-0512-bestpsf-nopsc' #+ fwhm + '-0256-front-nopsc'
    else:
        fn = 'fermi-allsky-' + en_bins_str[en_bin] + '-' + en_bins_str[en_bin+1] \
                        + 'GeV-fwhm' + fwhm + '-0512-bestpsf-nopsc' #+ fwhm + '-0256-front'
    return CTB_dir + fn + '.fits'
   
def get_CTB(CTB_dir, nside, CTB_en_min=0, CTB_en_max=num_en_bins):
    '''Input CTB_dir (where CTBCORE specbin maps are located), nside, and indices of lower (upper) edge of min (max) CTBCORE energy bin.
    Return arrays en_bins, count_maps, exposure_maps, psc_masks.
    en_bins gives edges of energy bins (in GeV).
    count_maps, exposure_maps, psc_masks gives counts, exposures (in cm^2 s), and point-source masks in energy bins at nside.'''
    if CTB_en_min > 0:
        map_nopsc_filename = fits_filename(CTB_dir, CTB_en_min-1, nopsc=True, smoothed=False)
        map_nopsc = fits.open(map_nopsc_filename)
        en_bins = [map_nopsc[0].header["EMAX"]]
        map_nopsc.close()
    else:
        en_bins = [0.3]
    count_maps = []
    exposure_maps = []
    psc_masks = []
    for en_bin in range(CTB_en_min, CTB_en_max):
        map_nopsc_filename = fits_filename(CTB_dir, en_bin, nopsc=True, smoothed=False)
        map_psc_filename = fits_filename(CTB_dir, en_bin, nopsc=False, smoothed=False)
        map_nopsc = fits.open(map_nopsc_filename)
        map_psc = fits.open(map_psc_filename)
        en_bin = map_nopsc[0].header["EMAX"]
        en_bins.append(en_bin)
        count_maps.append(read_CTB_count_map(map_nopsc_filename, nside))
        exposure_maps.append(hp.ud_grade(map_nopsc[1].data, nside))
        psc_masks.append(hp.ud_grade(map_psc[1].data == 0, nside, power=-2))
        map_nopsc.close()
        map_psc.close()
    logging.info('CTBCORE energy-bin edges, count maps, exposure maps, and point-source masks retrieved.')
    return np.array(en_bins), np.array(count_maps), np.array(exposure_maps), np.array(psc_masks)
    
def get_CTB_total_count_maps(CTB_dir, nside, min_en_bin=0, max_en_bin=num_en_bins):
    '''Input CTB_dir (where CTBCORE specbin maps are located), nside, and indices of lower (upper) edge of min (max) CTBCORE energy bin.
    Return total count maps (summed over energy bins) for psc-masked, unmasked, psc-masked+smoothed, and unmasked+smoothed CTBCORE maps.
    NOTE: SHOULD REMASK THESE MAPS WITH PSC-MASK OF MIN_EN_BIN TO REMOVE COUNTS WITHIN THAT MASK COMING FROM HIGHER ENERGIES (WHERE PSF IS SMALLER).'''
    npix = 12*nside**2
    total_map               = np.zeros(npix)    #total map (summed over energy bins) with point-source removal
    total_map_nopsc         = np.zeros(npix)    #total map with no point-source removal
    total_map_smth          = np.zeros(npix)    #total smoothed map with point-source removal
    total_map_nopsc_smth    = np.zeros(npix)    #total smoothed map with no point source removal
    #load maps in all energy bins and sum to get integrated map   
    for en_bin in range(min_en_bin, max_en_bin):
        total_map += read_CTB_count_map(fits_filename(CTB_dir, en_bin, nopsc=False, smoothed=False), nside)
        total_map_nopsc += read_CTB_count_map(fits_filename(CTB_dir, en_bin, nopsc=True, smoothed=False), nside)
        total_map_smth += read_CTB_count_map(fits_filename(CTB_dir, en_bin, nopsc=False, smoothed=True), nside)
        total_map_nopsc_smth += read_CTB_count_map(fits_filename(CTB_dir, en_bin, nopsc=True, smoothed=True), nside)
    logging.info('Total maps (summed over all energy bins) loaded.')
    return total_map, total_map_nopsc, total_map_smth, total_map_nopsc_smth    
