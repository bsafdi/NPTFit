import sys, os

# current_dir = os.getcwd()
# change_path = ".."
# os.chdir(change_path)

# This file sets up and executres the scan
import healpy as hp
import numpy as np
import config.set_dirs as sd
import pulsars.masks as masks
import pulsars.CTB as CTB
import config.make_templates as mt
import config.compute_PSF as CPSF

import pymultinest
import triangle
import matplotlib.pyplot as plt

import likelihood as llh

import logging

from config.config_file import config
from config.NPTF import NPTF
import config.NPTF_models as NPTF_models
from config.bayesian_scan import bayesian_scan

from config.analyze_results import analyze_results as ar


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
                NPTF_models.NPTF_ps(self.NPT_dist_compressed_ary, self.psf_dir, self.CTB_masked_compressed[en],
                                      k_max=self.k_max_array[en],  nside=self.nside, nbreak=self.nbreak) for en in range(self.number_energy_bins)]

        self.ll = self.log_like_NPTF_edep

    def initiate_poissonian_edep(self):
        self.configure_for_edep_scan()
        self.ll = self.log_like_poisson_edep
