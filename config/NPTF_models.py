import sys, os

# current_dir = os.getcwd()
# change_path = ".."
# os.chdir(change_path)

# This file sets up and executes the scan
import healpy as hp
import numpy as np
import math
import config.set_dirs as sd
import pulsars.masks as masks
import pulsars.CTB as CTB
import pulsars.likelihood_psf_pdep_new as llpsfn
import pulsars.likelihood_psf_pdep_2break as llpsfb
import pulsars.likelihood_psf_pdep_3break as llpsfb3
import config.make_templates as mt
import config.compute_PSF as CPSF

import pymultinest
import triangle
import matplotlib.pyplot as plt

import likelihood as llh

import logging

from config.config_file import config

from config.NPTF import NPTF

import time


class NPTF_1_ps(NPTF):
    def __init__(self, PS_dist_compressed, nexp, *args, **kwargs):

        self.nexp = nexp
        self.PS_dist_compressed = PS_dist_compressed

        NPTF.__init__(self, *args, **kwargs)
        if self.use_exposure:
            self.set_ps_exposure()
            if self.nbreak == 1:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_1break_exposure
            elif self.nbreak == 2:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_2break_exposure
            elif self.nbreak == 3:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_3break_exposure
        else:
            self.set_ps()
            if self.nbreak == 1:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_1break
            elif self.nbreak == 2:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_2break
            elif self.nbreak == 3:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_3break

    def set_ps(self):

        self.PS_dist_compressed = self.PS_dist_compressed[self.k_max_where]

    def set_ps_exposure(self):

        for i in range(self.nexp):
            self.PS_dist_compressed[i] = np.array(self.PS_dist_compressed[i])[np.array(self.k_max_where_exposure[i])]

    def make_ll_PDF_PSF_1break(self):
        """
        Computes the log of the likelihood
        :return:log (likelihood) of apparently just 1 Point source template
        """
        # TODO: Streamline the 1 PS, etc... to be in the same function
        # TODO: On ipython, there is a type matching error here, but it runs generally. Fix that!
        # print "CTB_map_compressed is ", self.CTB_map_compressed
        # print "the bkg is ", self.xbg_PSF_compressed
        ll = llpsfn.log_nu_k_ary_PSF_exact_1_PS(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                self.df_rho_div_f_ary, self.PS_dist_compressed,
                                                np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc)

        # print "the theta PS is ", self.theta_PS
        # print "The likelihood function is ", ll
        return ll

    def make_ll_PDF_PSF_1break_exposure(self):
        """
        Computes the log of the likelihood
        :param self.nexp: number of exposure regions. By default taken to be 1 which should give the same result as before
        :return:log (likelihood) of apparently just 1 Point source template
        """
        # TODO: Streamline the 1 PS, etc... to be in the same function
        # TODO: On ipython, there is a type matching error here, but it runs generally. Fix that!

        ll = 0.0

        for i in range(self.nexp):

            theta_PS = [self.theta_PS[0] * self.exposure_mean / self.exposure_means_list[i]] + list(
                self.theta_PS[1:3]) + \
                       [self.theta_PS[-1] * self.exposure_means_list[i] / self.exposure_mean]

            ll += llpsfn.log_nu_k_ary_PSF_exact_1_PS(self.xbg_PSF_compressed[i], np.array(theta_PS), self.f_ary,
                                                     self.df_rho_div_f_ary, self.PS_dist_compressed[i],
                                                     np.array(self.CTB_map_compressed[i], dtype='int32'), Sc=self.Sc)

        return ll


    def make_ll_PDF_PSF_2break(self):
        ll = llpsfb.log_nu_k_ary_PSF_exact_1_PS_2break(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                       self.df_rho_div_f_ary, self.PS_dist_compressed,
                                                       np.array(self.CTB_map_compressed, dtype='int32'))
        return ll

    def make_ll_PDF_PSF_2break_exposure(self):
        """
        Computes the log of the likelihood
        :return:log (likelihood) of apparently just 1 Point source template with 2 breaks
        """

        ll = 0.0
        for i in range(self.nexp):
            ####BS: fixed the line below to exposure correct both breaks!! make sure this works properly.
            theta_PS = [self.theta_PS[0] * self.exposure_mean / self.exposure_means_list[i]] + list(
                self.theta_PS[1:4]) + \
                       [self.theta_PS[-2] * self.exposure_means_list[i] / self.exposure_mean,self.theta_PS[-1] * self.exposure_means_list[i] / self.exposure_mean]

            ll += llpsfb.log_nu_k_ary_PSF_exact_1_PS_2break(self.xbg_PSF_compressed[i], np.array(theta_PS), self.f_ary,
                                                            self.df_rho_div_f_ary, self.PS_dist_compressed[i],
                                                            np.array(self.CTB_map_compressed[i], dtype='int32'))

        return ll

    def make_ll_PDF_PSF_3break(self):
        ll = llpsfb3.log_nu_k_ary_PSF_exact_1_PS_3break(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                       self.df_rho_div_f_ary, self.PS_dist_compressed,
                                                       np.array(self.CTB_map_compressed, dtype='int32'))
        return ll

    def make_ll_PDF_PSF_3break_exposure(self):
        """
        Computes the log of the likelihood
        :return:log (likelihood) of apparently just 1 Point source template with 2 breaks
        """

        ll = 0.0
        for i in range(self.nexp):
            ####BS: fixed the line below to exposure correct both breaks!! make sure this works properly.
            theta_PS = [self.theta_PS[0] * self.exposure_mean / self.exposure_means_list[i]] + list(
                self.theta_PS[1:5]) + \
                       [self.theta_PS[-3] * self.exposure_means_list[i] / self.exposure_mean, self.theta_PS[-2] * self.exposure_means_list[i] / self.exposure_mean,self.theta_PS[-1] * self.exposure_means_list[i] / self.exposure_mean]

            ll += llpsfb3.log_nu_k_ary_PSF_exact_1_PS_3break(self.xbg_PSF_compressed[i], np.array(theta_PS), self.f_ary,
                                                            self.df_rho_div_f_ary, self.PS_dist_compressed[i],
                                                            np.array(self.CTB_map_compressed[i], dtype='int32'))

        return ll

class NPTF_2_ps(NPTF):
    def __init__(self, PS_1_dist_compressed, PS_2_dist_compressed, *args, **kwargs):

        NPTF.__init__(self, *args, **kwargs)

        self.PS_1_dist_compressed = PS_1_dist_compressed
        self.PS_2_dist_compressed = PS_2_dist_compressed

        self.pre_load_1_PS = False

        self.set_ps()

    def set_ps(self):
        self.PS_1_dist_compressed = self.PS_1_dist_compressed[self.k_max_where]
        self.PS_2_dist_compressed = self.PS_2_dist_compressed[self.k_max_where]

    def init_1_PS(self, theta_1_PS_true):
        theta_1_PS_true = np.array(theta_1_PS_true)
        self.x_m_ary2, self.x_m_sum2 = llpsfn.return_xs(theta_1_PS_true, self.f_ary, self.df_rho_div_f_ary,
                                                        self.PS_2_dist_compressed,
                                                        np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc)

        self.pre_load_1_PS = True

    def make_ll_PDF_PSF(self):

        if self.pre_load_1_PS:
            ll = llpsfn.log_nu_k_ary_PSF_exact_2_PS(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                    self.df_rho_div_f_ary, self.PS_1_dist_compressed,
                                                    self.PS_2_dist_compressed,
                                                    np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc,
                                                    x_m_sum2_t=self.x_m_sum2, x_m_ary2_t=self.x_m_ary2)
        else:
            ll = llpsfn.log_nu_k_ary_PSF_exact_2_PS(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                    self.df_rho_div_f_ary, self.PS_1_dist_compressed,
                                                    self.PS_2_dist_compressed,
                                                    np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc)

        return ll


class NPTF_3_ps(NPTF):
    def __init__(self, PS_1_dist_compressed, PS_2_dist_compressed, PS_3_dist_compressed, *args, **kwargs):

        NPTF.__init__(self, *args, **kwargs)

        self.PS_1_dist_compressed = PS_1_dist_compressed
        self.PS_2_dist_compressed = PS_2_dist_compressed
        self.PS_3_dist_compressed = PS_3_dist_compressed

        self.pre_load_1_PS = False
        self.pre_load_2_PS = False

        self.set_ps()

    def set_ps(self):
        self.PS_1_dist_compressed = self.PS_1_dist_compressed[self.k_max_where]
        self.PS_2_dist_compressed = self.PS_2_dist_compressed[self.k_max_where]
        self.PS_3_dist_compressed = self.PS_3_dist_compressed[self.k_max_where]

    def init_1_PS(self, theta_1_PS_true):
        theta_1_PS_true = np.array(theta_1_PS_true)
        self.x_m_ary3, self.x_m_sum3 = llpsfn.return_xs(theta_1_PS_true, self.f_ary, self.df_rho_div_f_ary,
                                                        self.PS_3_dist_compressed,
                                                        np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc)

        self.pre_load_1_PS = True

    def init_2_PS(self, theta_1_PS_true, theta_2_PS_true):
        theta_1_PS_true = np.array(theta_1_PS_true)
        theta_2_PS_true = np.array(theta_2_PS_true)
        self.x_m_ary3, self.x_m_sum3 = llpsfn.return_xs(theta_1_PS_true, self.f_ary, self.df_rho_div_f_ary,
                                                        self.PS_3_dist_compressed,
                                                        np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc)
        self.x_m_ary2, self.x_m_sum2 = llpsfn.return_xs(theta_2_PS_true, self.f_ary, self.df_rho_div_f_ary,
                                                        self.PS_2_dist_compressed,
                                                        np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc)

        self.pre_load_2_PS = True

    def make_ll_PDF_PSF(self):

        if self.pre_load_2_PS:
            ll = llpsfn.log_nu_k_ary_PSF_exact_3_PS(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                    self.df_rho_div_f_ary, self.PS_1_dist_compressed,
                                                    self.PS_2_dist_compressed, self.PS_3_dist_compressed,
                                                    np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc,
                                                    x_m_sum2_t=self.x_m_sum2, x_m_ary2_t=self.x_m_ary2,
                                                    x_m_sum3_t=self.x_m_sum3, x_m_ary3_t=self.x_m_ary3)
        elif self.pre_load_1_PS:
            ll = llpsfn.log_nu_k_ary_PSF_exact_3_PS(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                    self.df_rho_div_f_ary, self.PS_1_dist_compressed,
                                                    self.PS_2_dist_compressed, self.PS_3_dist_compressed,
                                                    np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc,
                                                    x_m_sum3_t=self.x_m_sum3, x_m_ary3_t=self.x_m_ary3)
        else:
            # print 'Theta_ps: ', self.theta_PS
            # print 'len theta_PS: ', len(self.theta_PS)
            # print 'xbg: ', len(self.xbg_PSF_compressed),np.mean(self.xbg_PSF_compressed)
            # print 'PS1: ', len(self.PS_1_dist_compressed),np.mean(self.PS_1_dist_compressed)
            # print 'PS1: ', len(self.PS_1_dist_compressed),np.mean(self.PS_1_dist_compressed)
            # print 'PS2: ', len(self.PS_2_dist_compressed),np.mean(self.PS_2_dist_compressed)

            ll = llpsfn.log_nu_k_ary_PSF_exact_3_PS(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                    self.df_rho_div_f_ary, self.PS_1_dist_compressed,
                                                    self.PS_2_dist_compressed, self.PS_3_dist_compressed,
                                                    np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc)
        # print 'll = ', ll

        return ll
