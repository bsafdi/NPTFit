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
        # self.use_exposure = False

    def initiate_1_ps_edep_new_from_f_ary(self, f_ary, df_rho_div_f_ary, nbreak=1):
        self.initiate_1_ps_edep(0, f_ary=f_ary, df_rho_div_f_ary=df_rho_div_f_ary, nbreak=nbreak)

    def initiate_1_ps_edep(self, sigma_psf_deg, nbreak=1, **kwargs):
        if self.use_exposure:
            self.initiate_1_ps_edep_exposure(sigma_psf_deg, nbreak, **kwargs)
        else:

            self.nbreak = nbreak
            self.configure_for_edep_scan()

            if 'f_ary' and 'df_rho_div_f_ary' in kwargs.keys():
                f_ary = kwargs['f_ary']
                df_rho_div_f_ary = kwargs['df_rho_div_f_ary']
                kwarg_array = [{'f_ary': f_ary[en], 'df_rho_div_f_ary': df_rho_div_f_ary[en]} for en in
                               range(self.number_energy_bins)]
            else:
                kwarg_array = [{'sigma_psf_deg': sigma_psf_deg[en]} for en in range(self.number_energy_bins)]

            if 'PS_1_dist' in kwargs.keys():
                PS_1_dist = kwargs['PS_1_dist']  # summed, but not masked or compressed
                self.PS_1_dist_compressed = self.return_masked([PS_1_dist], compressed=True)[0]
            elif 'PS_1_dist_compressed' in kwargs.keys():
                self.PS_1_dist_compressed = kwargs['PS_1_dist_compressed']
            else:
                self.PS_1_dist_compressed = self.templates_dict_nested[self.non_poiss_models.keys()[0]]['summed_templates']


            self.NPTF_inst_array = [
                    NPTF_models.NPTF_1_ps(self.PS_1_dist_compressed, self.nexp, self.psf_dir, self.CTB_masked_compressed[en],
                                          k_max=self.k_max_array[en], Sc=self.Sc, nside=self.nside, nbreak=self.nbreak,
                                          **kwarg_array[en]) for en in range(self.number_energy_bins)]

            # self.NPTF_inst.set_ps()

            self.ll = self.log_like_NPTF_edep

    def initiate_1_ps_edep_exposure(self, sigma_psf_deg, nbreak=1, **kwargs):
        self.nbreak = nbreak
        self.exposure_pixels_compressed()
        self.configure_for_edep_scan()
        self.use_exposure = True

        if 'f_ary' and 'df_rho_div_f_ary' in kwargs.keys():
            f_ary = kwargs['f_ary']
            df_rho_div_f_ary = kwargs['df_rho_div_f_ary']
            kwarg_array = [{'f_ary': f_ary[en], 'df_rho_div_f_ary': df_rho_div_f_ary[en]} for en in
                           range(self.number_energy_bins)]
        else:
            kwarg_array = [{'sigma_psf_deg': sigma_psf_deg[en]} for en in range(self.number_energy_bins)]

        if 'PS_1_dist' in kwargs.keys():
            PS_1_dist = kwargs['PS_1_dist']  # summed, but not masked or compressed
            self.PS_1_dist_compressed_exposure = self.return_masked_exposure([PS_1_dist], compressed=True)[0]
            # TODO: Check the statement right before
        elif 'PS_1_dist_compressed' in kwargs.keys():
            self.PS_1_dist_compressed_exposure = kwargs['PS_1_dist_compressed']
        else:
            # The hackiest worst way to flip the energy and the number of regions
            self.PS_1_dist_compressed_exposure = self.templates_dict_nested[self.non_poiss_models.keys()[0]]['flux_summed_templates_exposure']


        self.NPTF_inst_array = [
                NPTF_models.NPTF_1_ps(self.PS_1_dist_compressed_exposure, self.nexp, self.psf_dir,
                                      np.array(self.CTB_masked_compressed_exposure)[::,en],
                                      k_max=self.k_max_array[en], Sc=self.Sc, nside=self.nside, nbreak=self.nbreak,
                                      use_exposure=True, exposure=self.exposure_means_list,
                                      exposure_mean=self.exposure_mean, **kwarg_array[en])
                for en in range(self.number_energy_bins)]


        # self.NPTF_inst.set_ps()

        self.ll = self.log_like_NPTF_edep

    def initiate_2_ps_edep_new_from_f_ary(self, f_ary, df_rho_div_f_ary, nbreak=1):
        self.initiate_2_ps_edep(0, f_ary=f_ary, df_rho_div_f_ary=df_rho_div_f_ary, nbreak=nbreak)

    def initiate_2_ps_edep(self, sigma_psf_deg, nbreak=1, **kwargs):
        self.nbreak = nbreak
        self.configure_for_edep_scan()

        if 'f_ary' and 'df_rho_div_f_ary' in kwargs.keys():
            f_ary = kwargs['f_ary']
            df_rho_div_f_ary = kwargs['df_rho_div_f_ary']
            kwarg_array = [{'f_ary': f_ary[en], 'df_rho_div_f_ary': df_rho_div_f_ary[en]} for en in
                           range(self.number_energy_bins)]
        else:
            kwarg_array = [{'sigma_psf_deg': sigma_psf_deg[en]} for en in range(self.number_energy_bins)]

        self.PS_1_dist_compressed = self.templates_dict_nested[self.non_poiss_models.keys()[0]]['summed_templates']
        if len(self.fixed_ps_model_dict.keys()) > 0:
            PS_2_dist_uncompressed = self.fixed_ps_model_dict.values()[0]['template']
            PS_2_dist_compressed = hp.ma(PS_2_dist_uncompressed)
            PS_2_dist_compressed.mask = self.mask_total
            self.PS_2_dist_compressed = PS_2_dist_compressed.compressed()
        else:
            self.PS_2_dist_compressed = self.templates_dict_nested[self.non_poiss_models.keys()[1]]['summed_templates']

        self.NPTF_inst_array = [
            NPTF_models.NPTF_2_ps(self.PS_1_dist_compressed, self.PS_2_dist_compressed, self.psf_dir,
                                  self.CTB_masked_compressed[en], k_max=self.k_max_array[en], Sc=self.Sc,
                                  nside=self.nside, **kwarg_array[en]) for en in range(self.number_energy_bins)]

        if len(self.fixed_ps_model_dict.keys()) > 0:
            print 'we are using a fixed PS template with medians ', self.fixed_ps_model_dict.values()[0]['medians']
            for NPTF_inst, ebin in map(None, self.NPTF_inst_array, range(self.number_energy_bins)):
                NPTF_inst.init_1_PS(
                    self.single_model_parameters_edep(self.fixed_ps_model_dict.values()[0]['medians'], energy_bin=ebin))

        # self.NPTF_inst.set_ps()

        self.ll = self.log_like_NPTF_edep

    def initiate_3_ps_edep_new_from_f_ary(self, f_ary, df_rho_div_f_ary, nbreak=1):
        self.initiate_3_ps_edep(0, f_ary=f_ary, df_rho_div_f_ary=df_rho_div_f_ary, nbreak=nbreak)

    def initiate_3_ps_edep(self, sigma_psf_deg, nbreak=1, **kwargs):
        self.nbreak = nbreak
        self.configure_for_edep_scan()

        if 'f_ary' and 'df_rho_div_f_ary' in kwargs.keys():
            f_ary = kwargs['f_ary']
            df_rho_div_f_ary = kwargs['df_rho_div_f_ary']
            kwarg_array = [{'f_ary': f_ary[en], 'df_rho_div_f_ary': df_rho_div_f_ary[en]} for en in
                           range(self.number_energy_bins)]
        else:
            kwarg_array = [{'sigma_psf_deg': sigma_psf_deg[en]} for en in range(self.number_energy_bins)]

        self.PS_1_dist_compressed = self.templates_dict_nested[self.non_poiss_models.keys()[0]]['summed_templates']

        if len(self.fixed_ps_model_dict.keys()) == 0:
            self.PS_2_dist_compressed = self.templates_dict_nested[self.non_poiss_models.keys()[1]]['summed_templates']
            self.PS_3_dist_compressed = self.templates_dict_nested[self.non_poiss_models.keys()[2]]['summed_templates']
        else:
            PS_3_dist_uncompressed = self.fixed_ps_model_dict.values()[0]['template']
            PS_3_dist_compressed = hp.ma(PS_3_dist_uncompressed)
            PS_3_dist_compressed.mask = self.mask_total
            self.PS_3_dist_compressed = PS_3_dist_compressed.compressed()
        if len(self.fixed_ps_model_dict.keys()) == 1:
            self.PS_2_dist_compressed = self.templates_dict_nested[self.non_poiss_models.keys()[1]]['summed_templates']
        if len(self.fixed_ps_model_dict.keys()) > 1:
            PS_2_dist_uncompressed = self.fixed_ps_model_dict.values()[1]['template']
            PS_2_dist_compressed = hp.ma(PS_2_dist_uncompressed)
            PS_2_dist_compressed.mask = self.mask_total
            self.PS_2_dist_compressed = PS_2_dist_compressed.compressed()

        self.NPTF_inst_array = [
            NPTF_models.NPTF_3_ps(self.PS_1_dist_compressed, self.PS_2_dist_compressed, self.PS_3_dist_compressed,
                                  self.psf_dir, self.CTB_masked_compressed[en], k_max=self.k_max_array[en], Sc=self.Sc,
                                  nside=self.nside, **kwarg_array[en]) for en in range(self.number_energy_bins)]

        if len(self.fixed_ps_model_dict.keys()) == 1:
            print 'we are using a single fixed PS template--2 floating templates--with medians ', \
            self.fixed_ps_model_dict.values()[0]['medians']
            for NPTF_inst, ebin in map(None, self.NPTF_inst_array, range(self.number_energy_bins)):
                NPTF_inst.init_1_PS(
                    self.single_model_parameters_edep(self.fixed_ps_model_dict.values()[0]['medians'], energy_bin=ebin))
        elif len(self.fixed_ps_model_dict.keys()) == 2:
            print 'we are using two fixed PS templates--1 floating templates--with medians ', \
            self.fixed_ps_model_dict.values()[0]['medians'], self.fixed_ps_model_dict.values()[1]['medians']
            for NPTF_inst, ebin in map(None, self.NPTF_inst_array, range(self.number_energy_bins)):
                NPTF_inst.init_2_PS(
                    self.single_model_parameters_edep(self.fixed_ps_model_dict.values()[0]['medians'], energy_bin=ebin),
                    self.single_model_parameters_edep(self.fixed_ps_model_dict.values()[1]['medians'], energy_bin=ebin))

        # self.NPTF_inst.set_ps()

        self.ll = self.log_like_NPTF_edep

    def initiate_poissonian_edep(self):
        self.configure_for_edep_scan()
        self.ll = self.log_like_poisson_edep
