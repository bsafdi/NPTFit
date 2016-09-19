import sys, os

import healpy as hp
import numpy as np

from config.bayesian_scan import bayesian_scan

sys.path.append("/group/hepheno/smsharma/NPTFit/Testing")
from NPLL import log_like

class bayesian_scan_NPTF(bayesian_scan):
    """ The file that sets up and executes the scan """

    def __init__(self, *args, **kwargs):

        bayesian_scan.__init__(self, *args, **kwargs)

    def initiate_ps_from_f_ary(self, f_ary, df_rho_div_f_ary, nbreak=1):
        self.nbreak = nbreak
        self.configure_for_edep_scan()
        self.f_ary = f_ary
        self.df_rho_div_f_ary = df_rho_div_f_ary

        self.num_non_poiss_models = len(self.non_poiss_models)

        self.NPT_dist_compressed_ary = [self.templates_dict_nested[self.non_poiss_models.keys()[i]]['summed_templates'] for i in range(self.num_non_poiss_models)] 

        self.NPTF_inst_array = [
                NPTF_ps(self.NPT_dist_compressed_ary, self.f_ary[en], self.df_rho_div_f_ary[en], self.CTB_masked_compressed[en], self.psf_dir, 
                                    nside=self.nside, nbreak=self.nbreak) for en in range(self.number_energy_bins)]

        self.ll = self.log_like_NPTF_edep

    def initiate_poissonian_edep(self):
        self.configure_for_edep_scan()
        self.ll = self.log_like_poisson_edep

class NPTF_ps():
    def __init__(self, PS_dist_compressed, f_ary, df_rho_div_f_ary, CTB_map_compressed,  *args, **kwargs):
        self.PS_dist_compressed = PS_dist_compressed
        self.f_ary = f_ary
        self.df_rho_div_f_ary = df_rho_div_f_ary
        self.CTB_map_compressed = CTB_map_compressed

    def make_ll_PDF_PSF(self):
        # Setup basic arrays - setup so can adjust how good the fit is by multiplying xbg and PS_dist

        ll = log_like(self.xbg_PSF_compressed, self.theta_PS, self.f_ary,
                                                       self.df_rho_div_f_ary, self.PS_dist_compressed,
                                                       np.array(self.CTB_map_compressed, dtype='int32'))
        return ll
