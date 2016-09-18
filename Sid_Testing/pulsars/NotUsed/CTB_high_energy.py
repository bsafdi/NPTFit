import os
import numpy as np
import healpy as hp
from astropy.io import fits

import logging

logger = logging.getLogger(__name__)


class CTB:
	""" A class to load the data for different configurations. Based on an earlier version of the file CTB by Nick Rodd"""

	def __init__(self, high_energy=False, new_en_set=[], dataset = 'PASS8_July15_'):
		
		self.dataset = dataset

		if new_en_set == []:
			if high_energy:
				self.en_bins_str = self.getting_energy_bins(high_energy = True)
				self.dataset = 'PASS8_May16_'
				self.lowest_en = 50
			else:
				self.en_bins_str = self.getting_energy_bins()
				self.lowest_en = 0.3
		else:
			self.en_bins_str = new_en_set 
			self.lowest_en = new_en_set[0]

		self.num_en_bins = len(self.en_bins_str) - 1

	@staticmethod
	def getting_energy_bins(high_energy = False):
		
		if not high_energy:
			# energy bins for original CTBCORE filenames
			en_bins_str = ['000.3', '000.4', '000.5', '000.6', '000.8', '000.9', '001.2', '001.5', \
						   '001.9', '002.4', '003.0', '003.8', '004.8', '006.0', '007.5', '009.5', '011.9', '015.0', '018.9',
						   '023.8', '030.0', '037.8', '047.5', '059.9', '075.4', '094.9', '119.4', '150.4', '189.3', '238.3',
						   '300.0']
		else:
			en_bins_str = ['050.0' , '072.3', '104.6', '151.2', '218.7', '316.2', '457.3', '661.3', '956.3', '01383.0', '02000.0']
		return en_bins_str


#['050.0-072.3' , '072.3-104.6', '104.6-151.2', '151.2-218.7', '218.7-316.2', '316.2-457.3', '457.3-661.3', '661.3-956.3' ]

	@staticmethod
	def read_CTB_count_map(input_filename, nside):
		# convert original CTBCORE intensity (counts/cm^2/s) to counts map at nside
		# open original .fits file
		f = fits.open(input_filename)
		# [0] field contains intensity, [1] field contains exposure -- use to get counts and convert to nside
		count_map = hp.ud_grade(np.around(f[0].data * f[1].data * 4 * np.pi / float(len(f[0].data))), nside, power=-2)
		f.close()
		return count_map

	def fits_filename(self,CTB_dir, en_bin, nopsc=False, smoothed=True, is_p8=True, is_ps_model=False, is_ps_mask=False, eventclass=5, eventtype=3):
		# function to help build .fits filenames depending on energy bin, point-source subtraction, and smoothing
		en_bins_str = self.en_bins_str
		if smoothed:
			fwhm = '120'
		else:
			fwhm = '000'
		if is_p8:
			if is_ps_model:
				fn = 'fermi-allsky-' + en_bins_str[en_bin] + '-' + en_bins_str[en_bin + 1] \
					 + 'GeV-fwhm' + fwhm + '-0512-bestpsf-pscmdl'
			elif is_ps_mask:
				fn = 'fermi-allsky-' + en_bins_str[en_bin] + '-' + en_bins_str[en_bin + 1] \
					 + 'GeV-fwhm' + fwhm + '-0512-bestpsf-mask'
			else:
				fn = 'fermi-allsky-' + en_bins_str[en_bin] + '-' + en_bins_str[en_bin + 1] \
					 + 'GeV-fwhm' + fwhm + '-0512-bestpsf-nopsc'  # + fwhm + '-0256-front-nopsc'
		elif nopsc:
			fn = 'fermi-allsky-' + en_bins_str[en_bin] + '-' + en_bins_str[en_bin + 1] \
				 + 'GeV-fwhm' + fwhm + '-0256-front-nopsc'
		else:
			fn = 'fermi-allsky-' + en_bins_str[en_bin] + '-' + en_bins_str[en_bin + 1] \
				 + 'GeV-fwhm' + fwhm + '-0256-front'
		return CTB_dir + fn + '.fits'



	def fits_filename_newstyle(self, CTB_dir, en_bin, nopsc, smoothed, is_p8, is_ps_model=False, is_ps_mask=False, eventclass=5, eventtype=3):
		# function to help build .fits filenames depending on energy bin, point-source subtraction, and smoothing
		# First establish which directory - depends on eventclass. Only setup for ultracleanveto (5) and source (2) at the moment
		if eventclass == 2:
			subdir = self.dataset + 'Source/specbin/'
		if eventclass == 5:
			subdir = self.dataset + 'UltracleanVeto/specbin/'
		if eventclass == 4:
			subdir = self.dataset + 'Ultraclean/specbin/'
		# Next determine which datatype - for ultracleanveto, only option 3 (Q4/bestpsf) is currently implemented, whilst for source 0, 3, 4 and 5 are available
		if eventtype == 0:
			datatype = ''
		if eventtype == 3:
			datatype = '-bestpsf'
		if eventtype == 4:
			datatype = '-psftop2'
		if eventtype == 5:
			datatype = '-psftop3'
		if smoothed:
			fwhm = '120'
		else:
			fwhm = '000'
		# Removed is_p8 keyword as all data we use now is pass8
		if is_ps_model:
			suffix = '-pscmdl'
		elif is_ps_mask:
			suffix = '-mask'
		else:
			suffix = '-nopsc'
		fn = 'fermi-allsky-' + self.en_bins_str[en_bin] + '-' + self.en_bins_str[en_bin + 1] \
			 + 'GeV-fwhm' + fwhm + '-0512' + datatype + suffix
		return CTB_dir + subdir + fn + '.fits'


	def get_CTB(self, CTB_dir, nside, CTB_en_min=0, CTB_en_max=29, is_p8=False, eventclass=5, eventtype=3, newstyle=0):
		'''Input CTB_dir (where CTBCORE specbin maps are located), nside, and indices of lower (upper) edge of min (max) CTBCORE energy bin.
		Return arrays en_bins, count_maps, exposure_maps, psc_masks.
		en_bins gives edges of energy bins (in GeV).
		count_maps, exposure_maps, psc_masks gives counts, exposures (in cm^2 s), and point-source masks in energy bins at nside.'''

		if CTB_en_min > 0:
			if newstyle:
				map_nopsc_filename = self.fits_filename_newstyle(CTB_dir, CTB_en_min - 1, nopsc=True, smoothed=False,
															is_p8=is_p8, eventclass=eventclass, eventtype=eventtype)
			else:
				map_nopsc_filename = self.fits_filename(CTB_dir, CTB_en_min - 1, nopsc=True, smoothed=False, is_p8=is_p8, eventclass=eventclass, eventtype=eventtype)
			map_nopsc = fits.open(map_nopsc_filename)
			en_bins = [map_nopsc[0].header["EMAX"]]
			map_nopsc.close()
		else:
			en_bins = [self.lowest_en]
		count_maps = []
		exposure_maps = []
		psc_masks = []
		for en_bin in range(CTB_en_min, CTB_en_max):
			if newstyle:
				map_nopsc_filename = self.fits_filename_newstyle(CTB_dir, en_bin, nopsc=True, smoothed=False, is_p8=is_p8, eventclass=eventclass, eventtype=eventtype)
				
				map_psc_filename = self.fits_filename_newstyle(CTB_dir, en_bin, nopsc=False, smoothed=False, is_p8=is_p8, eventclass=eventclass, eventtype=eventtype)
			else:
				map_nopsc_filename = self.fits_filename(CTB_dir, en_bin, nopsc=True, smoothed=False, is_p8=is_p8,
												   eventclass=eventclass, eventtype=eventtype)
				map_psc_filename = self.fits_filename(CTB_dir, en_bin, nopsc=False, smoothed=False, is_p8=is_p8,
												 eventclass=eventclass, eventtype=eventtype)
			map_nopsc = fits.open(map_nopsc_filename)
			map_psc = fits.open(map_psc_filename)
			en_bin = map_nopsc[0].header["EMAX"]
			en_bins.append(en_bin)
			count_maps.append(self.read_CTB_count_map(map_nopsc_filename, nside))
			exposure_maps.append(hp.ud_grade(map_nopsc[1].data, nside))
			psc_masks.append(hp.ud_grade(map_psc[1].data == 0, nside, power=-2))
			map_nopsc.close()
			map_psc.close()
		logging.info('CTBCORE energy-bin edges, count maps, exposure maps, and point-source masks retrieved.')
		return np.array(en_bins), np.array(count_maps), np.array(exposure_maps), np.array(psc_masks)


	def get_ps_models(self, CTB_dir, nside, CTB_en_min=0, CTB_en_max=29, eventclass=5, eventtype=3, newstyle=0):
		ps_models = []
		for en_bin in range(CTB_en_min, CTB_en_max):
			if newstyle:
				f = fits.open(
					self.fits_filename_newstyle(CTB_dir, en_bin, 0, 0, is_p8=True, is_ps_model=True, eventclass=eventclass,
										   eventtype=eventtype))
				f_mask = fits.open(
					self.fits_filename_newstyle(CTB_dir, en_bin, 0, 0, is_p8=True, eventclass=eventclass, eventtype=eventtype))
			else:
				f = fits.open(self.fits_filename(CTB_dir, en_bin, 0, 0, is_p8=True, is_ps_model=True, eventclass=eventclass,
											eventtype=eventtype))
				f_mask = fits.open(
					self.fits_filename(CTB_dir, en_bin, 0, 0, is_p8=True, eventclass=eventclass, eventtype=eventtype))
			count_map = hp.ud_grade(np.around(f[0].data * f_mask[1].data * 4 * np.pi / float(len(f[0].data))), nside,
									power=-2)
			ps_models.append(count_map)
			f.close()
			f_mask.close()
		return np.array(ps_models)


	def get_ps_masks(self, CTB_dir, nside, CTB_en_min=0, CTB_en_max=29, eventclass=5, eventtype=3, newstyle=0):
		ps_masks = []
		for en_bin in range(CTB_en_min, CTB_en_max):
			if newstyle:
				f = fits.open(
					self.fits_filename_newstyle(CTB_dir, en_bin, 0, 0, is_p8=True, is_ps_mask=True, eventclass=eventclass,
										   eventtype=eventtype))
			else:
				f = fits.open(self.fits_filename(CTB_dir, en_bin, 0, 0, is_p8=True, is_ps_mask=True, eventclass=eventclass,
											eventtype=eventtype))
			ps_masks.append(hp.ud_grade(f[0].data, nside, power=-2))
			f.close()
		return np.array(ps_masks)


	def get_CTB_total_count_maps(self, CTB_dir, nside, min_en_bin=0, max_en_bin=29, is_p8=False, eventclass=5,
								 eventtype=3, newstyle=0):
		'''Input CTB_dir (where CTBCORE specbin maps are located), nside, and indices of lower (upper) edge of min (max) CTBCORE energy bin.
		Return total count maps (summed over energy bins) for psc-masked, unmasked, psc-masked+smoothed, and unmasked+smoothed CTBCORE maps.
		NOTE: SHOULD REMASK THESE MAPS WITH PSC-MASK OF MIN_EN_BIN TO REMOVE COUNTS WITHIN THAT MASK COMING FROM HIGHER ENERGIES (WHERE PSF IS SMALLER).'''
		npix = 12 * nside ** 2
		total_map = np.zeros(npix)  # total map (summed over energy bins) with point-source removal
		total_map_nopsc = np.zeros(npix)  # total map with no point-source removal
		total_map_smth = np.zeros(npix)  # total smoothed map with point-source removal
		total_map_nopsc_smth = np.zeros(npix)  # total smoothed map with no point source removal
		# load maps in all energy bins and sum to get integrated map
		for en_bin in range(min_en_bin, max_en_bin):
			if newstyle:
				total_map += read_CTB_count_map(
					self.fits_filename_newstyle(CTB_dir, en_bin, nopsc=False, smoothed=False, is_p8=is_p8, eventclass=eventclass,
										   eventtype=eventtype), nside)
				total_map_nopsc += read_CTB_count_map(
					self.fits_filename_newstyle(CTB_dir, en_bin, nopsc=True, smoothed=False, is_p8=is_p8, eventclass=eventclass,
										   eventtype=eventtype), nside)
				total_map_smth += read_CTB_count_map(
					self.fits_filename_newstyle(CTB_dir, en_bin, nopsc=False, smoothed=True, is_p8=is_p8, eventclass=eventclass,
										   eventtype=eventtype), nside)
				total_map_nopsc_smth += read_CTB_count_map(
					self.fits_filename_newstyle(CTB_dir, en_bin, nopsc=True, smoothed=True, is_p8=is_p8, eventclass=eventclass,
										   eventtype=eventtype), nside)
			else:
				total_map += read_CTB_count_map(
					self.fits_filename(CTB_dir, en_bin, nopsc=False, smoothed=False, is_p8=is_p8, eventclass=eventclass,
								  eventtype=eventtype), nside)
				total_map_nopsc += read_CTB_count_map(
					self.fits_filename(CTB_dir, en_bin, nopsc=True, smoothed=False, is_p8=is_p8, eventclass=eventclass,
								  eventtype=eventtype), nside)
				total_map_smth += read_CTB_count_map(
					self.fits_filename(CTB_dir, en_bin, nopsc=False, smoothed=True, is_p8=is_p8, eventclass=eventclass,
								  eventtype=eventtype), nside)
				total_map_nopsc_smth += read_CTB_count_map(
					self.fits_filename(CTB_dir, en_bin, nopsc=True, smoothed=True, is_p8=is_p8, eventclass=eventclass,
								  eventtype=eventtype), nside)
		logging.info('Total maps (summed over all energy bins) loaded.')
		return total_map, total_map_nopsc, total_map_smth, total_map_nopsc_smth
