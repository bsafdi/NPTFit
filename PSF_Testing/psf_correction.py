#------------------------------------------------------------------------------
# NPTF PSF correction
#------------------------------------------------------------------------------
##
## Add f_ary and df_rho_div_f_ary to self, both of which are needed to
## compute the psf correction to the NPTF calculation.
##
## If the arrays already exist they are loaded, otherwise the computation is
## passed onto psf_calculation and returned.
##
## Code defaults to a gaussian psf unless specified otherwise.
##
#------------------------------------------------------------------------------

import os
import numpy as np
import psf_compute

class psf_correction:
    def __init__(self, num_f_bins=10, n_psf=50000, n_pts_per_psf=1000, 
                 f_trunc=0.01, psf_sigma_deg=0.1812, delay_compute = False):
        """ Load or calculate the NPTF PSF correction
            Adds f_ary and df_rho_div_ary to self
            Defaults to a gaussian psf, for user defined psfs see below

            Args:
                num_f_bins:    number of bins to calculate rho(f) in
                n_psf:         number of psfs to place when calculating rho(f)
                n_pts_per_psf: number of points to place per psf in calculation
                f_trunc:       minimum flux fraction to keep track of
                psf_sigma_deg: 1 sigma containment of the gaussian psf
                delay_compute: set to true to define a custom psf
            
            Custom PSF details:
                To implement a custom (non-gaussian) psf create an instance of
                psf_correction with the above parameters (ignore psf_sigma_deg)
                setting delay_compute = True

                Then redefine the following attributes:
                    psf_r_func:     The psf as a function of r, distance from
                                    the center of the point
                    sample_psf_max: maximum distance to sample the psf from
                                    the center, should be larger for psfs
                                    with long tails
                    psf_samples:    number of samples to make from the psf
                                    (linearly spaced) from 0 to sample_psf_max,
                                    should be large enough to adequately
                                    represent the full psf
                    psf_tag:        label the psf is saved with

                Finally execute the calculation using make_or_load_psf_corr()
            """
        
        self.num_f_bins = num_f_bins
        self.n_psf = n_psf
        self.n_pts_per_psf = n_pts_per_psf
        self.f_trunc = f_trunc

        # TODO: remove this hack
        self.psf_dir = '/zfs/nrodd/NPTFit/PSF_Testing/psf_dir/'
        self.nside = 128

        # Convert psf from degrees to radians
        # Only used if the psf is a gaussian
        psf_sigma = psf_sigma_deg*np.pi/180.

        # If using a custom PSF, set delay_compute = True and manually 
        # set the following 4 variables for that instance of psf_correction 
        self.psf_r_func = lambda r: np.exp(-r**2 / (2.*psf_sigma**2))
        self.sample_psf_max = 5. * psf_sigma
        self.psf_samples = 10000
        self.psf_tag = 'gauss_'+str(self.nside)+'_'+str(np.round(psf_sigma_deg, 3))

        
        # If set delay evaluation so that a user defined PSF can be entered
        # User then has to manually call this to implement the calculation
        self.psf_corr_file = self.psf_dir + self.psf_tag + '.npy'
        if not delay_compute:
            self.make_or_load_psf_corr()

    def make_or_load_psf_corr(self):
        """ Function to load or calculate f_ary and df_rho_div_f_ary and 
            append them to self
            """
        if not os.path.exists(self.psf_corr_file):
            print "Computing the psf correction ..."
            self.f_ary, self.df_rho_div_f_ary \
            = psf_compute.psf_corr(self.nside, self.num_f_bins, self.n_psf, 
                                   self.n_pts_per_psf, self.f_trunc, 
                                   self.psf_r_func, self.sample_psf_max, 
                                   self.psf_samples)
            print "... done"
            tosave = np.array([self.f_ary, self.df_rho_div_f_ary])
            np.save(self.psf_corr_file,tosave)
            print "File saved as:",self.psf_corr_file
        else:
            print "Loading the psf correction from:",self.psf_corr_file
            loadpsf = np.load(self.psf_corr_file)
            self.f_ary = loadpsf[0]
            self.df_rho_div_f_ary = loadpsf[1]
