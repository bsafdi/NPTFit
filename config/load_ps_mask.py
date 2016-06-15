import sys, os

# current_dir = os.getcwd()
# change_path = ".."
# os.chdir(change_path)

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


class load_ps_mask:
	def __init__(self,mask_dir,en,mask_type = '0.99',nside=128):
		self.en = en
		self.mask_dir = mask_dir
		self.mask_type = mask_type
		self.nside = nside

		self.make_string()

		self.load_map()


	def make_string(self):
		if self.en < 1.:
			self.round = '%.8f'
		elif self.en < 10.:
			self.round = '%.7f'
		elif self.en < 100.:
			self.round = '%.6f'
		else:
			self.round = '%.5f'
		self.en_string = str( self.round % self.en)
		self.the_string = 'pscmask_3FGL-energy' + self.en_string + '_' + self.mask_type + '_Q2_ULTRACLEAN_FRONT.fits'

	def load_map(self):
		self.mask_fits = fits.open(self.mask_dir + self.the_string  )
		self.the_ps_mask_raw_nside = self.mask_fits[0].data

		self.the_ps_mask = np.logical_not(np.logical_not(hp.ud_grade(self.the_ps_mask_raw_nside,self.nside,dtype=float)))






