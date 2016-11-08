import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib
from scipy.interpolate import interp1d as interp1d


class smooth:
	def __init__(self,map_to_smooth):
		self.map_to_smooth = map_to_smooth
		self.npix = len(self.map_to_smooth)
		self.nside = hp.npix2nide(self.npix)
		
		pass


	def smear_map(self,psf,mult_sigma_for_smooth): #psf in sigma degrees
        self.non_zero_vals = np.where(self.counts_map > 0)[0]
        self.swp_inst = swp.smooth_gaussian_psf(psf,self.counts_map, mult_sigma_for_smooth = mult_sigma_for_smooth)
        self.counts_map_psf = np.sum([self.swp_inst.smooth_the_pixel(nzv) for nzv in self.non_zero_vals],axis=0)