###############################################################################
# nptfit.py
###############################################################################
#
# Head file for performing bayesian scans using (non-)Poissonian template fits.
#
###############################################################################

from __future__ import absolute_import

import numpy as np

from .nptf_scan import NPTFScan  # Sets up and performs the scan


class NPTF(NPTFScan):
    def __init__(self, tag='Untagged', work_dir=None):
        # Initialise NPTFScan
        NPTFScan.__init__(self, tag=tag, work_dir=work_dir)

    def configure_for_scan(self, f_ary=np.array([1.0]), 
                           df_rho_div_f_ary=np.array([1.0]), nexp=1):
        """ Set up a non-poissonian scan (or poissonian if no NP templates)

            Takes PSF parameters as input -- determined in PSFCorrection

            :param f_ary: Photon leakage probabilities characterizing PSF,
            sum(f_ary) = 1.0
            :param df_rho_div_f_ary: df*rho(f)/f for integrating over f as a sum
            Defaults to these being a delta function, which is the case if
            there is no PSF correction
            :param nexp: Number of exposure regions to calculate the NPTF within
        """

        self.f_ary = f_ary
        self.df_rho_div_f_ary = df_rho_div_f_ary
        self.nexp = nexp

        # If no non-Poissonian models, use Poissonian likelihood
        if len(self.non_poiss_models) == 0:
            # Set number of exposure regions to 1 and define
            # the Poissonian likelihood
            self.nexp = 1
            self.ll = self.log_like_ptf
        else:
            # Define the non-Poissonian likelihood
            self.ll = self.log_like_nptf

        # Prepare likelihood, data and templates for the scan
        self.configure_for_scan_internal()
