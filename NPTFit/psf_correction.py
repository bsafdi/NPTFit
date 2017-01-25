###############################################################################
# psf_correction.py
###############################################################################
#
# Add f_ary and df_rho_div_f_ary to self, both of which are needed to
# compute the psf correction to the NPTF calculation.
#
# If the arrays already exist they are loaded, otherwise the computation is
# passed onto psf_calculation and returned.
#
# Code defaults to a gaussian psf unless specified otherwise.
#
# NB: this module fundamentally assumes that the analysis is being performed
# on a spherical region pixelised using HEALPix. If this is not the case, the
# PSF correction must be computed differently.
#
###############################################################################

from __future__ import print_function
from __future__ import absolute_import

import os
import numpy as np
from . import psf_compute


class PSFCorrection:
    def __init__(self, psf_dir=None, num_f_bins=10, n_psf=50000,
                 n_pts_per_psf=1000, f_trunc=0.01, nside=128,
                 psf_sigma_deg=None, delay_compute=False):
        """ Load or calculate the NPTF PSF correction

            Adds f_ary and df_rho_div_ary to self. Defaults to a Gaussian PSF,
            for user defined PSFs see below

            :param psf_dir: directory where psf files are stored
            :param num_f_bins: number of bins to calculate rho(f) in
            :param n_psf: number of psfs to place when calculating rho(f)
            :param n_pts_per_psf: number of points to place per psf in calculation
            :param f_trunc: minimum flux fraction to keep track of
            :param nside: nside of the map rho(f) is used on
            :param psf_sigma_deg: 1 sigma containment of the gaussian psf
            :param delay_compute: set to true to define a custom psf

            Custom PSF details:

            To implement a custom (non-gaussian) psf create an instance of
            PSFCorrection with the above parameters (ignore psf_sigma_deg)
            setting delay_compute = True

            Then redefine the following attributes:

            :param psf_r_func: The psf as a function of r, distance from
                   the center of the point
            :param sample_psf_max: maximum distance to sample the psf from
                   the center, should be larger for psfs with long tails
            :param psf_samples: number of samples to make from the psf
                   (linearly spaced) from 0 to sample_psf_max,
                   should be large enough to adequately represent the full psf
            :param psf_tag: label the psf is saved with

            Finally execute the calculation using make_or_load_psf_corr()
        """

        # User must either specify a sigma, or use a custom PSF
        if not delay_compute:
            assert(psf_sigma_deg is not None), \
                   "Must either specify a sigma or use a custom PSF"
        else:
            # If custom PSF give sigma temporary value
            if psf_sigma_deg is None:
                psf_sigma_deg = 0.1812

        self.psf_dir = psf_dir
        self.num_f_bins = num_f_bins
        self.n_psf = n_psf
        self.n_pts_per_psf = n_pts_per_psf
        self.f_trunc = f_trunc
        self.nside = nside

        if self.psf_dir is None:
            self.psf_dir = os.getcwd() + '/psf_dir/'
        self.make_dirs([self.psf_dir])

        # Convert psf from degrees to radians
        # Only used if the psf is a gaussian
        psf_sigma = psf_sigma_deg*np.pi/180.

        # If using a custom PSF, set delay_compute = True and manually
        # set the following 4 variables for that instance of PSFCorrection
        self.psf_r_func = lambda r: np.exp(-r**2 / (2.*psf_sigma**2))
        self.sample_psf_max = 5. * psf_sigma
        self.psf_samples = 10000
        self.psf_tag = 'gauss_'+str(self.nside) + '_' + \
                       str(np.round(psf_sigma_deg, 3)) + '_' + \
                       str(self.num_f_bins) + '_' + \
                       str(self.n_psf) + '_' + \
                       str(self.n_pts_per_psf) + '_' + \
                       str(self.f_trunc)

        # If set delay evaluation so that a user defined PSF can be entered
        # User then has to manually call this to implement the calculation
        if not delay_compute:
            self.make_or_load_psf_corr()

    def make_or_load_psf_corr(self):
        """ Function to load or calculate f_ary and df_rho_div_f_ary and
            append them to self
        """

        self.psf_corr_file = self.psf_dir + self.psf_tag + '.npy'
        if not os.path.exists(self.psf_corr_file):
            self.f_ary, self.df_rho_div_f_ary \
                = psf_compute.psf_corr(self.nside, self.num_f_bins, self.n_psf,
                                       self.n_pts_per_psf, self.f_trunc,
                                       self.psf_r_func, self.sample_psf_max,
                                       self.psf_samples)
            tosave = np.array([self.f_ary, self.df_rho_div_f_ary])
            # Check again if exists, for multicore runs
            if not os.path.exists(self.psf_corr_file):
                np.save(self.psf_corr_file, tosave)
                print("File saved as:", self.psf_corr_file)
        else:
            print("Loading the psf correction from:", self.psf_corr_file)
            loadpsf = np.load(self.psf_corr_file)
            self.f_ary = loadpsf[0]
            self.df_rho_div_f_ary = loadpsf[1]

    @staticmethod
    def make_dirs(dirs):
        """ Creates directories if they do not already exist
        """

        for d in dirs:
            if not os.path.exists(d):
                try:
                    os.mkdir(d)
                except OSError as e:
                    if e.errno != 17:
                        raise
