import sys, os

import healpy as hp
import numpy as np

from bayesian_scan import bayesian_scan
from NPLL import log_like

class bayesian_scan_NPTF(bayesian_scan):
    """ The file that sets up and executes the scan """

    def __init__(self, *args, **kwargs):

        bayesian_scan.__init__(self, *args, **kwargs)

    def initiate_ps_from_f_ary(self, f_ary, df_rho_div_f_ary, nbreak=1, nexp=2):
        self.nbreak = nbreak
        self.f_ary = f_ary
        self.df_rho_div_f_ary = df_rho_div_f_ary
        self.nexp = nexp

        # self.exposure_pixels_compressed()
        self.configure_for_edep_scan()

        self.NPT_dist_compressed_ary = [self.templates_dict_nested[self.non_poiss_models.keys()[i]]['template_masked_compressed'] for i in range(self.n_non_poiss_models)] 
        
        self.NPT_dist_compressed_exp_ary = [self.templates_dict_nested[self.non_poiss_models.keys()[i]]['template_masked_compressed_expreg'] for i in range(self.n_non_poiss_models)]

        self.NPTF_inst_array = [
                NPTF_ps(self.NPT_dist_compressed_exp_ary, self.f_ary, self.df_rho_div_f_ary, np.array(self.masked_compressed_data_expreg), 
                                    nside=self.nside, nbreak=self.nbreak, nexp = self.nexp, exposure_mean_list=self.exposure_means_list, exposure_mean=self.exposure_mean)]

        self.ll = self.log_like_NPTF_edep

    def initiate_poissonian_edep(self):
        self.configure_for_edep_scan()
        self.ll = self.log_like_poisson_edep

class NPTF_ps():
    def __init__(self, PS_dist_compressed, f_ary, df_rho_div_f_ary, compressed_data_array,  
                                    nside, nbreak, nexp, exposure_mean_list, exposure_mean):
        self.PS_dist_compressed = PS_dist_compressed
        self.f_ary = f_ary
        self.df_rho_div_f_ary = df_rho_div_f_ary
        self.compressed_data_array = compressed_data_array
        self.nside = nside
        self.nbreak = nbreak
        self.nexp = nexp
        self.exposure_means_list = exposure_mean_list
        self.exposure_mean = exposure_mean

    def make_ll_PDF_PSF(self):
        # Setup basic arrays - setup so can adjust how good the fit is by multiplying xbg and PS_dist

        ll = 0.0

        for i in range(self.nexp):

            theta_PS = [self.theta_PS[0] * self.exposure_mean / self.exposure_means_list[i]] + list(
                self.theta_PS[1:self.nbreak+2]) + \
                       [self.theta_PS[3+i] * self.exposure_means_list[i] / self.exposure_mean for i in range(self.nbreak)]

            if len(np.array(self.compressed_data_array.T[i], dtype='int32')) == 0: continue

            ll += log_like(self.xbg_PSF_compressed[i], self.theta_PS, self.f_ary,
                                                       self.df_rho_div_f_ary, [self.PS_dist_compressed[0][i]],
                                                       np.array(self.compressed_data_array.T[i], dtype='int32'))
        return ll
