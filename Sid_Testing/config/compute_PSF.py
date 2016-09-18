import numpy as np
import pulsars.psf as psf
import os

class compute_PSF:
    def __init__(self, sigma_PSF_deg, psf_dir, nside=128, num_f_bins=10):
        self.sigma_PSF_deg = sigma_PSF_deg
        self.psf_dir = psf_dir
        self.nside = nside
        self.num_f_bins = num_f_bins

        self.f_ary_file = self.psf_dir + 'f_ary-' + str(self.nside) + '-' + str(np.round(sigma_PSF_deg, 3)) + '-' + str(
            self.num_f_bins) + '.dat'

        self.make_or_load_psf()

    def make_or_load_psf(self):
        print 'The name of the PSF file is ', self.f_ary_file
        if not os.path.exists(self.f_ary_file):
            print 'we have to make PSF ...'
            self.f_ary, self.df_rho_div_f_ary = psf.make_PSF(self.psf_dir, self.nside, self.sigma_PSF_deg)
            print 'finished making PSF ...'
        else:
            print 'just loading PSF ...'
            self.f_ary, self.df_rho_div_f_ary = psf.load_PSF(self.psf_dir, self.nside, self.sigma_PSF_deg)
