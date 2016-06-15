# coding: utf-8

# In[1]:

import numpy as np
import scipy
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

import healpy as hp
import mpmath as mp
from astropy.io import fits

# import os

# current_dir = os.getcwd()
# change_path = ".."
# os.chdir(change_path)

# import pulsars.special as spc
import pulsars.masks as masks
import pulsars.CTB as CTB
import pulsars.psf as psf
import pulsars.diffuse_fermi as df
# import pulsars.likelihood_psf_pdep as llpsf
# os.chdir(current_dir)




import logging
from contextlib import contextmanager
import sys, os


# subroutines

# ======================================
# define PDF and PDF+PSF functions

# variables for creating PSF
do_make_PSF = 0
sigma_PSF_deg = 0.18 / 1.5  # check this!!!!
sigma_PSF = sigma_PSF_deg * np.pi / 180.

# psf_dir = 'psf/'



###subroutines


# This is from Mathematica fit
num = 0.321531
exp = 0.756881


##===================
def find_sim(spect, CTB_en_bins, En):
    if En < CTB_en_bins[0]:
        return spect[0]
    for i in range(np.size(CTB_en_bins) - 1):
        if En < CTB_en_bins[i + 1] and En > CTB_en_bins[i]:
            return spect[i]
    return spect[-1]


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

# class compute_PSF_kings:
#     def __init__(self,psf_dir,fcore,score,gcore,stail,gtail,SpE,frac=[1],nside=128,num_f_bins=10,save_tag='temp'):
#         self.sigma_PSF_deg = sigma_PSF_deg
#         self.psf_dir = psf_dir  
#         self.nside = nside
#         self.num_f_bins = num_f_bins
#         self.save_tag=save_tag

#         self.f_ary_file = self.psf_dir + 'f_ary-' + str(self.nside) + '-' + str(np.round(sigma_PSF_deg,3)) + '-' + str(self.num_f_bins) + '.dat'

#         self.make_or_load_psf()

#     def make_or_load_psf(self):
#         print 'The name of the PSF file is ', self.f_ary_file
#         if not os.path.exists(self.f_ary_file):
#             print 'we have to make PSF ...'
#             self.f_ary, self.df_rho_div_f_ary = psf.make_PSF(self.psf_dir,self.nside, self.sigma_PSF_deg)
#             print 'finished making PSF ...'
#         else:
#             print 'just loading PSF ...'
#             self.f_ary, self.df_rho_div_f_ary = psf.load_PSF(self.psf_dir, self.nside, self.sigma_PSF_deg)


# #Make PSF
# def main(CTB_en_bins,NSIDE,psf_dir,spect = 'False',sigma_PSF_deg = 'Null',num_f_bins=10,just_sigma = False,use_lowest_bin = False):

#     if spect == 'False':
#         spect = np.ones(len(CTB_en_bins))

#     if sigma_PSF_deg=='Null':
#         numBins = 10000
#         deltaE = (CTB_en_bins[-1] - CTB_en_bins[0])/numBins
#         CTB_en_bins_f = np.array([CTB_en_bins[0] + deltaE*i for i in range(numBins)])
#         if use_lowest_bin == True:
#             sigma_PSF_deg = num/CTB_en_bins_f[0]**exp #degrees
#         else:
#             spect_long = [find_sim(spect, CTB_en_bins,CTB_en_bins_f[i])for i in range(np.size(CTB_en_bins_f)-1)]
#             array2 = np.array([deltaE*spect_long[i]*num/CTB_en_bins_f[i]**exp for i in range(np.size(CTB_en_bins_f)-1)])
#             array3 = np.array([deltaE*spect_long[i] for i in range(np.size(CTB_en_bins_f)-1)])
#             sigma_PSF_deg = np.sum(array2)/np.sum(array3) #degrees

#     sigma_PSF = sigma_PSF_deg*np.pi/180.

#     if just_sigma == True:
#         return sigma_PSF_deg
#     else:
#         f_ary_file = psf_dir + 'f_ary-' + str(NSIDE) + '-' + str(np.round(sigma_PSF_deg,3)) + '-' + str(num_f_bins) + '.dat'

#         if not os.path.exists(f_ary_file):
#             print 'we have to make PSF ...'
#             f_ary, df_rho_div_f_ary = psf.make_PSF(psf_dir, NSIDE, sigma_PSF_deg)
#             print 'finished making PSF ...'
#         else:
#             print 'just loading PSF ...'
#             f_ary, df_rho_div_f_ary = psf.load_PSF(psf_dir, NSIDE, sigma_PSF_deg)

#         return f_ary, df_rho_div_f_ary
