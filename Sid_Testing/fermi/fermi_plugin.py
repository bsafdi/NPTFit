import sys, os

# current_dir = os.getcwd()
# change_path = ".."
# os.chdir(change_path)

# This file stores input
import healpy as hp
import numpy as np
from config.set_dirs import set_dirs as sd
import pulsars.masks as masks
#import pulsars.CTB as CTB
import config.make_templates as mt
import config.compute_PSF as CPSF

from config.make_NFW import make_NFW
from config.load_ps_mask import load_ps_mask as lpm
import config.smooth_with_psf as swp

import fermi.PSF_class as PSFC
import fermi.PSF_kings_class as PSFKC

import config.make_f_ary_kings as mfak

from astropy.io import fits

import logging

from pulsars.CTB_high_energy import CTB as CTBclass


class fermi_plugin:
    """ A plugin for Fermi data """

    def __init__(self, maps_dir, fermi_data_dir='False', work_dir='', data_name='p8', CTB_en_min=8, CTB_en_max=16,
                 nside=128, eventclass=5, eventtype=3, newstyle=0, highenergy=False, total_en_bins = []):
        # eventclass and eventtype are variables to determine which dataset you want to load
        # eventclass has the following options:
        #   1: Transient
        #   2: Source
        #   3: Clean
        #   4: Ultraclean
        #   5: UltracleanVeto
        # NB: only 2 and 5 are currently implemented - this is dealt with in pulsars/CTB.py
        # eventtype has the following options:
        #   0: All photons (Q4)
        #   1: Front-converting events
        #   2: Back-converting events
        #   3: Q1 - top quartile by PSF (also "bestpsf")
        #   4: Q2 - top two quartiles by PSF
        #   5: Q3 - top three quartiles by PSF
        # NB: for ultracleanveto, only option 3 (Q1/bestpsf) is currently implemented,
        # whilst for source 0, 3, 4 and 5 are available


        # Lina, May 2016: Adding the high energy data, 50 GeV- 2 TeV in 10 Logarithmic bins
        # These will be called bins 50-60 to keep far away from the previous ones

        self.work_dir = work_dir
        self.fermi_data_dir = fermi_data_dir

        self.nside = nside
        self.npix = hp.nside2npix(self.nside)

        self.eventclass = eventclass
        self.eventtype = eventtype
        self.newstyle = newstyle  # This is a keyword to activate the new style whilst we keep the old around for legacy reasons - Nick Rodd, 22 Sep 15

        # Parameters for the CTBCORE data
        self.maps_dir = maps_dir
        self.data_name = data_name
        self.set_data_tag(self.data_name)  # also can be p15_ultraclean_Q2_specbin/
        self.set_initial_dirs()

        self.CTB_en_min = CTB_en_min
        self.CTB_en_max = CTB_en_max
        self.nEbins = self.CTB_en_max - self.CTB_en_min

        self.additional_template_dict = {}
        self.templates = []
        self.comp_array = []
        self.template_dict = {}


        self.total_en_bins = 10 ** np.linspace(np.log10(0.3), np.log10(300), 31)
        #High energy data
        self.highenergy = highenergy
        if self.highenergy:
            self.total_en_bins = 10 ** np.linspace(np.log10(50), np.log10(2000), 11)

        if not (total_en_bins == []):
            self.total_en_bins = total_en_bins


        self.maxenergy = len(self.total_en_bins) - 2 #maximum energy bin

        # So far this is set up for the high energy run, but this file can be modified to accomodate any other energy. No changes are necessary for CTB_high_energy.py in pulsars
        self.CTB = CTBclass(high_energy = self.highenergy) #self.CTB_dir, 
        # evaluate functions
        self.load_CTBCORE_data()
        self.load_data_ps_model_and_maps()
        #self.make_ps_mask()

    def set_initial_dirs(self):

        if self.fermi_data_dir == 'False':
            self.CTB_dir = self.maps_dir + self.data_tag + '/'
            self.diff_dir = self.maps_dir + 'diffuse/'
        else:
            self.CTB_dir = self.fermi_data_dir
            self.diff_dir = self.fermi_data_dir  # Data and maps should be in the same directory
        self.ps_dir = self.maps_dir + '3FGL_masks_by_energy/'
        self.ps_flux_dir = self.maps_dir + '3FGL/'
        self.template_dir = self.maps_dir + 'additional_templates/'
        self.psf_data_dir = self.maps_dir + 'psf_data/'
        self.data_dir = self.work_dir + 'data/'
        if not os.path.exists(self.data_dir):
            os.mkdir(self.data_dir)

    def set_data_tag(self, data='old_CTBCORE'):  # flag: old
        self.is_p8 = False
        if data == 'old_CTBCORE':
            self.data_tag = 'p15_ultraclean_Q2_specbin'
        elif data == 'new_CTBCORE':
            self.data_tag = 'Apr3_Ultraclean_Q2'
        elif data == 'p8':
            self.is_p8 = True
            self.data_tag = 'PASS8_Jun15_UltracleanVeto_specbin'
        # elif data == 'p8_high_E': # Lina May 2016, Adding the high energy data
        #     self.is_p8 = True
        #     self.data_tag = 'PASS8_May16_Source_specbin'

    def load_CTBCORE_data(self):
        self.CTB_en_bins, self.CTB_count_maps, self.CTB_exposure_maps, self.CTB_psc_masks_temp = self.CTB.get_CTB(
            self.CTB_dir, self.nside, self.CTB_en_min, self.CTB_en_max, is_p8=self.is_p8, eventclass=self.eventclass, eventtype=self.eventtype, newstyle=self.newstyle)
        self.CTB_map_um = np.sum(self.CTB_count_maps, axis=0)
        self.CTB_count_psc_maps = [hp.ma(self.CTB_count_map) for self.CTB_count_map in self.CTB_count_maps]
        for i in range(len(self.CTB_count_psc_maps)):
            self.CTB_count_psc_maps[i].mask = self.CTB_psc_masks_temp[i]

        self.total_exposure_map = np.mean(self.CTB_exposure_maps, axis=0)
        self.NPIX_total = np.size(self.total_exposure_map)

    def load_data_ps_model_and_maps(self):
        if not self.is_p8:
            pass
        else:
            self.CTB_psc_masks = self.CTB.get_ps_masks(self.CTB_dir, self.nside, self.CTB_en_min, self.CTB_en_max,
                                                  eventclass=self.eventclass, eventtype=self.eventtype,
                                                  newstyle=self.newstyle)
            self.CTB_ps_models = self.CTB.get_ps_models(self.CTB_dir, self.nside, self.CTB_en_min, self.CTB_en_max,
                                                   eventclass=self.eventclass, eventtype=self.eventtype,
                                                   newstyle=self.newstyle)

            self.CTB_psc_masks_total = self.CTB.get_ps_masks(self.CTB_dir, self.nside, 0, self.maxenergy, eventclass=self.eventclass,
                                                        eventtype=self.eventtype, newstyle=self.newstyle)

    def load_template_from_file(self, filename,
                                alttempdir='False'):  # should already be an exposure-corrected counts map
        # By default load templates from template_dir:
        if alttempdir == 'False':
            fits_file_list = fits.open(self.template_dir + filename)[0].data
            return np.array([hp.ud_grade(fits_file, self.nside, power=-2) for fits_file in fits_file_list])[
                   self.CTB_en_min:self.CTB_en_max]
        else:
            # This option is currently only setup for a single file
            fits_file = fits.open(alttempdir + filename)[0].data
            return np.array([hp.ud_grade(fits_file, self.nside, power=-2)])

    def load_template_combination_from_file(self, filename, alttempdir='False', excluded=-1):
        # This is a function to load and combine a list of templates from one fits file
        if alttempdir == 'False':
            temp_list = fits.open(self.template_dir + filename)
        else:
            temp_list = fits.open(alttempdir + filename)
        forrange = range(len(temp_list))
        if excluded != -1:
            del forrange[excluded]
        combined_temp = np.zeros(self.npix)
        for i in forrange:
            combined_temp += temp_list[i].data
        return np.array([hp.ud_grade(combined_temp, self.nside, power=-2)])

    def reduce_to_energy_subset(self, en_minimum, en_maximum):
        en_maximum_p_1 = en_maximum + 1

        if self.is_p8:
            self.CTB_psc_masks = self.CTB_psc_masks[en_minimum:en_maximum]
            self.CTB_ps_models = self.CTB_ps_models[en_minimum:en_maximum]
            self.CTB_psc_masks_total = self.CTB_psc_masks_total[en_minimum:en_maximum]

        self.CTB_exposure_maps = self.CTB_exposure_maps[en_minimum:en_maximum]
        self.CTB_count_maps = self.CTB_count_maps[en_minimum:en_maximum]
        self.CTB_en_bins = self.CTB_en_bins[en_minimum:en_maximum_p_1]
        self.CTB_map_um = np.sum(self.CTB_count_maps, axis=0)
        self.CTB_count_psc_maps = self.CTB_count_psc_maps[en_minimum:en_maximum]
        self.total_exposure_map

        self.CTB_en_min = self.CTB_en_min + en_minimum;
        self.CTB_en_max = self.CTB_en_min + en_maximum - en_minimum;
        self.nEbins = self.CTB_en_max - self.CTB_en_min

        for key in self.template_dict.keys():
            self.template_dict[key] = self.template_dict[key][en_minimum:en_maximum]

    def add_diffuse(self, comp='p6', **kwargs):
        if 'template_dir' in kwargs.keys():
            template_dir = kwargs['template_dir']
        else:
            template_dir = self.diff_dir
        self.comp_array.append(comp)
        is_p8 = False
        if comp == 'p6' or comp == 'p7' or comp == 'p8':
            if self.data_name == 'new_CTBCORE':
                self.diff_file = 'diffmodelsApr3_Ultraclean_Q2-front.fits'
            elif self.data_name == 'old_CTBCORE':
                self.diff_file = 'diffmodelsp15_ultraclean_Q2-front.fits'
            elif self.data_name == 'p8':
                if comp == 'p6':
                    load_comp = 'p6v11'
                else:
                    load_comp = comp
                self.diff_file = 'diffuse_model_map_healpix_' + load_comp + '_PASS8_Jun15_UltracleanVeto_bestpsf.fits'
                is_p8 = True

            # self.templates.append(mt.return_diff_um(self.diff_dir + self.diff_file ,self.CTB_en_bins,self.CTB_en_min,self.CTB_en_max,NSIDE=self.nside,mode=comp,is_p8 = is_p8))
            self.template_dict[comp] = mt.return_diff_um(template_dir + self.diff_file, self.CTB_en_bins,
                                                         self.CTB_en_min, self.CTB_en_max, NSIDE=self.nside, mode=comp,
                                                         is_p8=is_p8)

    def add_diffuse_newstyle(self, comp='p6', eventclass=5, eventtype=3):
        # This code can be used to add in diffuse, pibrem or ics templates
        # First establish which directory - depends on eventclass. Only setup for ultracleanveto (5) and source (2) at the moment
        if eventclass == 2:
            subdir = self.CTB.dataset + 'Source/specbin/'
            diffmidname = self.CTB.dataset + 'Source'
        if eventclass == 5:
            subdir = self.CTB.dataset + 'UltracleanVeto/specbin/'
            diffmidname = self.CTB.dataset + 'UltracleanVeto'
        if eventclass == 4:
            subdir = self.CTB.dataset + 'Ultraclean/specbin/'
            diffmidname = self.CTB.dataset + 'Ultraclean'

        template_dir = self.diff_dir + subdir
        self.comp_array.append(comp)

        # Next determine which datatype - for ultracleanveto, only option 3 (Q4/bestpsf) is currently implemented, whilst for source 0, 3, 4 and 5 are available
        if eventtype == 0:
            datatype = ''
        if eventtype == 3:
            datatype = '_bestpsf'
        if eventtype == 4:
            datatype = '_psftop2'
        if eventtype == 5:
            datatype = '_psftop3'

        # Adjust for the fact p6 is called p6v11 in files
        if comp == 'p6':
            load_comp = 'p6v11'
        else:
            load_comp = comp
        # Define Filename
        self.diff_file = 'diffuse_model_map_healpix_' + load_comp + '_' + diffmidname + datatype + '.fits'
        # Add in file
        self.template_dict[comp] = mt.return_diff_p8_um(template_dir + self.diff_file, self.CTB_en_bins,
                                                        self.CTB_en_min, self.CTB_en_max, NSIDE=self.nside)

    def add_bubbles(self, comp='bubs', **kwargs):
        if 'template_dir' in kwargs.keys():
            template_dir = kwargs['template_dir']
        else:
            template_dir = self.template_dir
        self.comp_array.append(comp)
        self.template_dict[comp] = mt.return_bubbles(template_dir + 'bubbles_intensity_maps.npy', self.CTB_en_min,
                                                     self.CTB_en_max, self.CTB_exposure_maps, NSIDE=self.nside)

    def add_iso(self, comp='iso'):
        self.comp_array.append(comp)
        self.template_dict[comp] = mt.return_isotropic(self.CTB_en_min, self.CTB_en_max, self.CTB_exposure_maps,
                                                       NSIDE=self.nside)

    def add_nfw(self, gamma_nfw=1.25, GC_lng=0, nfw_maps_dir='generic', comp='nfw'):
        if nfw_maps_dir == 'generic':
            nfw_maps_dir = self.work_dir + 'data/maps/'
        if not os.path.exists(nfw_maps_dir):
            os.mkdir(nfw_maps_dir)

        self.comp_array.append(comp)

        self.NFW = make_NFW(gamma_nfw=gamma_nfw, GC_lng=GC_lng, nside=self.nside, data_maps_dir=nfw_maps_dir)
        self.nfw_profile_flux = self.NFW.NFW_profile[0]  # just a single map
        self.nfw_profile = self.flux_to_counts(np.array([self.nfw_profile_flux for i in range(self.nEbins)]))
        self.template_dict[comp] = self.nfw_profile

    def flux_to_counts(self, the_map):
        return the_map * self.CTB_exposure_maps * 4 * np.pi / self.npix

    def add_ps(self, ell, b, sigma_psf_deg, rescale=1, comp='ps', excluded=-1):
        # if not os.path.exists(psf_dir):
        #     os.mkdir(psf_dir)
        self.comp_array.append(comp)
        # if spect == 'False':
        #     spect = np.ones(self.nEbins)
        psf = sigma_psf_deg  # CPSF.main(self.CTB_en_bins, self.nside, psf_dir, sigma_pdf_deg, just_sigma = True)
        print 'The psf for this template is ', psf
        ps_inst = mt.ps_template(ell, b, self.nside, psf)
        ps_temps = [rescale * ps_inst.smooth_ps_map for i in range(self.nEbins)]
        self.template_dict[comp] = ps_temps

    def add_multiple_ps(self, ell, b, sigma_psf_deg, rescale, comp='ps_comb', excluded=-1):
        self.comp_array.append(comp)
        psf = sigma_psf_deg
        forrange = range(len(ell))
        if excluded != -1:
            del forrange[excluded]
        combined_temp = np.zeros(self.npix)
        for j in forrange:
            ps_inst = mt.ps_template(ell[j], b[j], self.nside, psf)
            combined_temp += rescale[j] * ps_inst.smooth_ps_map
        self.template_dict[comp] = [combined_temp for i in range(self.nEbins)]

    def add_ps_king_fast(self, ell, b, rescale=1, comp='ps'):
        # add in a king function psf, which is the functional form of the Fermi PSF
        # This is a faster version of ps_king, but requires joblib to run. Only relevant if
        # loading many different point sources repeatedly
        self.comp_array.append(comp)
        load_ps_king = mt.ps_template_king_fast(ell, b, self.nside, self.maps_dir, self.CTB_en_min, self.eventclass,
                                                self.eventtype)
        ps_temps = [rescale * load_ps_king for i in range(self.nEbins)]
        self.template_dict[comp] = ps_temps

    def add_multiple_ps_king_fast(self, ell, b, rescale, comp='ps_comb', excluded=-1):
        # This is a faster version of ps_king, but requires joblib to run. Only relevant if
        # loading many different point sources repeatedly
        self.comp_array.append(comp)
        forrange = range(len(ell))
        if excluded != -1:
            del forrange[excluded]
        combined_temp = np.zeros(self.npix)
        for j in forrange:
            load_ps_king = mt.ps_template_king_fast(ell, b, self.nside, self.maps_dir, self.CTB_en_min, self.eventclass,
                                                    self.eventtype)
            combined_temp += rescale[j] * load_ps_king
        self.template_dict[comp] = [combined_temp for i in range(self.nEbins)]

    def add_ps_king(self, ell, b, rescale=1, comp='ps'):
        # add in a king function psf, which is the functional form of the Fermi PSF
        self.comp_array.append(comp)
        psk_inst = mt.ps_template_king(ell, b, self.nside, self.maps_dir, self.CTB_en_min, self.eventclass,
                                       self.eventtype)
        ps_temps = [rescale * psk_inst.smooth_ps_map for i in range(self.nEbins)]
        self.template_dict[comp] = ps_temps

    def add_multiple_ps_king(self, ell, b, rescale, comp='ps_comb', excluded=-1):
        self.comp_array.append(comp)
        forrange = range(len(ell))
        if excluded != -1:
            del forrange[excluded]
        combined_temp = np.zeros(self.npix)
        for j in forrange:
            psk_inst = mt.ps_template_king(ell[j], b[j], self.nside, self.maps_dir, self.CTB_en_min, self.eventclass,
                                           self.eventtype)
            combined_temp += rescale[j] * psk_inst.smooth_ps_map
        self.template_dict[comp] = [combined_temp for i in range(self.nEbins)]

    def add_ps_model(self, comp='ps_model'):
        self.comp_array.append(comp)
        self.load_data_ps_model_and_maps()
        self.template_dict[comp] = self.CTB_ps_models

    def add_template_from_file(self, comp, template_name, alttempdir='False'):
        self.comp_array.append(comp)
        template = self.load_template_from_file(template_name, alttempdir=alttempdir)
        self.template_dict[comp] = template

    def add_template_combination_from_file(self, comp, template_name, alttempdir='False', excluded=-1):
        self.comp_array.append(comp)
        template = self.load_template_combination_from_file(template_name, alttempdir=alttempdir, excluded=excluded)
        self.template_dict[comp] = template

    def add_template_by_hand(self, comp, template):
        self.comp_array.append(comp)
        self.template_dict[comp] = template

    def use_template_normalization_file(self, file_path, key_suffix='False', use_keys='False', dont_use_keys='False'):
        self.template_norm_array = np.load(file_path)
        self.norm_array_en_bins = [ar[0] for ar in self.template_norm_array]

        self.CTB_en_bin_centers = 10 ** np.array(
            [(np.log10(self.CTB_en_bins[i]) + np.log10(self.CTB_en_bins[i + 1])) / 2. for i in
             range(len(self.CTB_en_bins) - 1)])
        lower = np.where(np.round(self.CTB_en_bin_centers, 3)[0] == np.round(np.array(self.norm_array_en_bins), 3))[0]
        upper = np.where(np.round(self.CTB_en_bin_centers, 3)[-1] == np.round(np.array(self.norm_array_en_bins), 3))[0]
        if len(lower) > 0 and len(upper > 0):
            lower_index = lower[0]
            upper_index = upper[0]
            if use_keys == 'False':
                keys = self.template_dict.keys()
            else:
                keys = list(use_keys)
            if dont_use_keys != 'False':
                for key in dont_use_keys:
                    if key in keys:
                        keys.remove(key)
            print 'Using normalization keys (reduced) ', keys

            for i in range(lower_index, upper_index + 1):
                if key_suffix != 'False':
                    for key in keys:
                        if key + key_suffix in self.template_norm_array[i][1].keys():
                            self.template_dict[key][i - lower_index] = self.template_dict[key][i - lower_index] * \
                                                                       self.template_norm_array[i][1][key + key_suffix]

                    for key in keys:
                        if key + key_suffix in self.template_norm_array[i][1].keys():
                            print 'was able to use normalization file for key ', key, ' at energy ', \
                            self.CTB_en_bin_centers[i - lower_index]
                else:

                    for key in keys:
                        if key in self.template_norm_array[i][1].keys():
                            self.template_dict[key][i - lower_index] = self.template_dict[key][i - lower_index] * \
                                                                       self.template_norm_array[i][1][key]

                    for key in keys:
                        if key in self.template_norm_array[i][1].keys():
                            print 'was able to use normalization file for key ', key, ' at energy ', \
                            self.CTB_en_bin_centers[i - lower_index]


        else:
            print 'can\'t use this normalization file. sorry! try again.'

    def make_ps_mask(self, mask_type='0.99', force_energy=False, energy=0.0, energy_bin=-1):
        if mask_type == 'False':
            pass
        else:
            if energy_bin > -1:
                total_en_bins = self.total_en_bins #10 ** np.linspace(np.log10(0.3), np.log10(300), 31)
                energy = total_en_bins[energy_bin]
            if mask_type == 'top300':
                if not force_energy:
                    self.ps_mask_array = self.CTB_psc_masks
                else:
                    print 'Taking PS mask from energy', energy
                    energys = np.vectorize(round)(self.total_en_bins, 3)#np.vectorize(round)(10 ** np.linspace(np.log10(0.3), np.log10(300), 31), 3)
                    bin_number = np.where(energys == round(energy, 3))[0][0]
                    self.ps_mask_array = [self.CTB_psc_masks_total[bin_number] for en in self.CTB_en_bins]
            elif mask_type == 'old_ps':
                the_ps_mask = np.loadtxt(self.maps_dir + 'additional_masks/PS_mask5SL-8-16.dat')
                self.ps_mask_array = [the_ps_mask for en in self.CTB_en_bins]
            else:
                if not force_energy:
                    self.ps_mask_array = [lpm(self.ps_dir, en, mask_type=mask_type, nside=self.nside).the_ps_mask for en in self.CTB_en_bins]
                else:
                    self.ps_mask_array = [lpm(self.ps_dir, energy, mask_type=mask_type, nside=self.nside).the_ps_mask for en in self.CTB_en_bins]

    def load_psf(self, data_name='p8', fits_file_path='False'):
        if fits_file_path != 'False':
            self.fits_file_name = fits_file_path
        elif data_name == 'p8':
            # Define the param and rescale indices for the various quartiles
            params_index_psf1 = 10
            rescale_index_psf1 = 11
            params_index_psf2 = 7
            rescale_index_psf2 = 8
            params_index_psf3 = 4
            rescale_index_psf3 = 5
            params_index_psf4 = 1
            rescale_index_psf4 = 2
            # Setup to load the correct PSF details depending on the dataset, and define appropriate theta_norm values, psf1 is bestpsf, psf4 is worst psf (but just a quartile each, so for Q1-3 need tocombine
            if self.eventclass == 2:
                psf_file_name = 'psf_P8R2_SOURCE_V6_PSF.fits'
                theta_norm_psf1 = [0.0000000, 9.7381019e-06, 0.0024811595, 0.022328802, 0.080147663, 0.17148392,
                                   0.30634315, 0.41720551]
                theta_norm_psf2 = [0.0000000, 0.00013001938, 0.010239333, 0.048691643, 0.10790632, 0.18585539,
                                   0.29140913, 0.35576811]
                theta_norm_psf3 = [0.0000000, 0.00074299273, 0.018672204, 0.062317201, 0.12894928, 0.20150553,
                                   0.28339386, 0.30441893]
                theta_norm_psf4 = [4.8923139e-07, 0.011167475, 0.092594658, 0.15382001, 0.16862869, 0.17309118,
                                   0.19837774, 0.20231968]
            elif self.eventclass == 5:
                psf_file_name = 'psf_P8R2_ULTRACLEANVETO_V6_PSF.fits'
                theta_norm_psf1 = [0.0000000, 9.5028121e-07, 0.00094418357, 0.015514370, 0.069725775, 0.16437751,
                                   0.30868705, 0.44075016]
                theta_norm_psf2 = [0.0000000, 1.6070284e-05, 0.0048551576, 0.035358049, 0.091767466, 0.17568974,
                                   0.29916159, 0.39315185]
                theta_norm_psf3 = [0.0000000, 0.00015569366, 0.010164870, 0.048955837, 0.11750811, 0.19840060,
                                   0.29488095, 0.32993394]
                theta_norm_psf4 = [0.0000000, 0.0036816313, 0.062240006, 0.14027030, 0.17077023, 0.18329804, 0.21722594,
                                   0.22251374]
            elif self.eventclass == 4:
                # NB: Currently these values are just copied from ultracleanveto - need to recalculate these for ultraclean - see /zfs/tslatyer/galactic/fermi/pro/incidence_angle_hist.pro
                psf_file_name = 'psf_P8R2_ULTRACLEAN_V6_PSF.fits'
                theta_norm_psf1 = [0.0000000, 9.5028121e-07, 0.00094418357, 0.015514370, 0.069725775, 0.16437751,
                                   0.30868705, 0.44075016]
                theta_norm_psf2 = [0.0000000, 1.6070284e-05, 0.0048551576, 0.035358049, 0.091767466, 0.17568974,
                                   0.29916159, 0.39315185]
                theta_norm_psf3 = [0.0000000, 0.00015569366, 0.010164870, 0.048955837, 0.11750811, 0.19840060,
                                   0.29488095, 0.32993394]
                theta_norm_psf4 = [0.0000000, 0.0036816313, 0.062240006, 0.14027030, 0.17077023, 0.18329804, 0.21722594,
                                   0.22251374]
            self.fits_file_name = self.psf_data_dir + psf_file_name
        if fits_file_path != 'False' or data_name == 'p8':
            self.f = fits.open(self.fits_file_name)
            # Now need to get the correct PSF for the appropriate quartile. 
            # If anything other than best psf, need to combine quartiles.
            # Quartiles aren't exactly equal in size, but approximate as so.
            self.PSFC_inst = PSFC.PSF(self.f, theta_norm=theta_norm_psf1, rescale_index=rescale_index_psf1,
                                      params_index=params_index_psf1)
            calc_sigma_PSF_deg = np.array(self.PSFC_inst.return_sigma_gaussian(self.CTB_en_bins))
            if ((self.eventtype == 4) or (self.eventtype == 5) or (self.eventtype == 0)):
                self.PSFC_inst = PSFC.PSF(self.f, theta_norm=theta_norm_psf2, rescale_index=rescale_index_psf2,
                                          params_index=params_index_psf2)
                calc_sigma_load = np.array(self.PSFC_inst.return_sigma_gaussian(self.CTB_en_bins))
                calc_sigma_PSF_deg = (calc_sigma_PSF_deg + calc_sigma_load) / 2.
                if ((self.eventtype == 5) or (self.eventtype == 0)):
                    self.PSFC_inst = PSFC.PSF(self.f, theta_norm=theta_norm_psf3, rescale_index=rescale_index_psf3,
                                              params_index=params_index_psf3)
                    calc_sigma_load = np.array(self.PSFC_inst.return_sigma_gaussian(self.CTB_en_bins))
                    calc_sigma_PSF_deg = (2. * calc_sigma_PSF_deg + calc_sigma_load) / 3.
                    if self.eventtype == 0:
                        self.PSFC_inst = PSFC.PSF(self.f, theta_norm=theta_norm_psf4, rescale_index=rescale_index_psf4,
                                                  params_index=params_index_psf4)
                        calc_sigma_load = np.array(self.PSFC_inst.return_sigma_gaussian(self.CTB_en_bins))
                        calc_sigma_PSF_deg = (3. * calc_sigma_PSF_deg + calc_sigma_load) / 4.
            # Now take mean of array and extract first element, which will be the PSF
            self.sigma_PSF_deg = calc_sigma_PSF_deg
            self.average_PSF()
        else:
            pass

    def load_king_param(self):
        """ A function to combine king function params at each energy value """
        CTB_en_bin_min = np.array([self.CTB_en_bins[i] for i in range(len(self.CTB_en_bins) - 1)])

        fcore, score, gcore, stail, gtail, SpE = self.load_king_param_en(CTB_en_bin_min[0])
        fcore_arr = np.array([fcore])
        score_arr = np.array([score])
        gcore_arr = np.array([gcore])
        stail_arr = np.array([stail])
        gtail_arr = np.array([gtail])
        SpE_arr = np.array([SpE])
        if len(CTB_en_bin_min > 0):
            for i in (np.arange(len(CTB_en_bin_min) - 1) + 1):
                fcore, score, gcore, stail, gtail, SpE = self.load_king_param_en(CTB_en_bin_min[i])
                fcore_arr = np.append(fcore_arr, np.array([fcore]), axis=0)
                score_arr = np.append(score_arr, np.array([score]), axis=0)
                gcore_arr = np.append(gcore_arr, np.array([gcore]), axis=0)
                stail_arr = np.append(stail_arr, np.array([stail]), axis=0)
                gtail_arr = np.append(gtail_arr, np.array([gtail]), axis=0)
                SpE_arr = np.append(SpE_arr, np.array([SpE]), axis=0)
        self.fcore_arr = fcore_arr
        self.score_arr = score_arr
        self.gcore_arr = gcore_arr
        self.stail_arr = stail_arr
        self.gtail_arr = gtail_arr
        self.SpE_arr = SpE_arr

    def load_king_param_en(self, energyval):
        """ A function to calculate the various king function parameters at one energy """

        # Define relevant parameters
        params_index_psf1 = 10
        rescale_index_psf1 = 11
        params_index_psf2 = 7
        rescale_index_psf2 = 8
        params_index_psf3 = 4
        rescale_index_psf3 = 5
        params_index_psf4 = 1
        rescale_index_psf4 = 2
        if self.eventclass == 2:
            psf_file_name = 'psf_P8R2_SOURCE_V6_PSF.fits'
            theta_norm_psf1 = [0.0000000, 9.7381019e-06, 0.0024811595, 0.022328802, 0.080147663, 0.17148392, 0.30634315,
                               0.41720551]
            theta_norm_psf2 = [0.0000000, 0.00013001938, 0.010239333, 0.048691643, 0.10790632, 0.18585539, 0.29140913,
                               0.35576811]
            theta_norm_psf3 = [0.0000000, 0.00074299273, 0.018672204, 0.062317201, 0.12894928, 0.20150553, 0.28339386,
                               0.30441893]
            theta_norm_psf4 = [4.8923139e-07, 0.011167475, 0.092594658, 0.15382001, 0.16862869, 0.17309118, 0.19837774,
                               0.20231968]
        elif self.eventclass == 5:
            psf_file_name = 'psf_P8R2_ULTRACLEANVETO_V6_PSF.fits'
            theta_norm_psf1 = [0.0000000, 9.5028121e-07, 0.00094418357, 0.015514370, 0.069725775, 0.16437751,
                               0.30868705, 0.44075016]
            theta_norm_psf2 = [0.0000000, 1.6070284e-05, 0.0048551576, 0.035358049, 0.091767466, 0.17568974, 0.29916159,
                               0.39315185]
            theta_norm_psf3 = [0.0000000, 0.00015569366, 0.010164870, 0.048955837, 0.11750811, 0.19840060, 0.29488095,
                               0.32993394]
            theta_norm_psf4 = [0.0000000, 0.0036816313, 0.062240006, 0.14027030, 0.17077023, 0.18329804, 0.21722594,
                               0.22251374]

        psf_file = fits.open(self.maps_dir + 'psf_data/' + psf_file_name)

        # Now need to load as many parameters as there are quartiles
        # First load 1st quartile as always needed, then load others if relevant
        kparam = PSFKC.PSF_king(psf_file, theta_norm=theta_norm_psf1, rescale_index=rescale_index_psf1,
                                params_index=params_index_psf1)
        fcore = np.array([kparam.return_king_params(energyval, 'fcore')])
        score = np.array([kparam.return_king_params(energyval, 'score')])
        gcore = np.array([kparam.return_king_params(energyval, 'gcore')])
        stail = np.array([kparam.return_king_params(energyval, 'stail')])
        gtail = np.array([kparam.return_king_params(energyval, 'gtail')])
        SpE = np.array([kparam.rescale_factor(energyval)])
        if ((self.eventtype == 4) or (self.eventtype == 5) or (self.eventtype == 0)):
            # Add in 2nd quartile if required
            kparam = PSFKC.PSF_king(psf_file, theta_norm=theta_norm_psf2, rescale_index=rescale_index_psf2,
                                    params_index=params_index_psf2)
            fcore = np.append(fcore, np.array([kparam.return_king_params(energyval, 'fcore')]))
            score = np.append(score, np.array([kparam.return_king_params(energyval, 'score')]))
            gcore = np.append(gcore, np.array([kparam.return_king_params(energyval, 'gcore')]))
            stail = np.append(stail, np.array([kparam.return_king_params(energyval, 'stail')]))
            gtail = np.append(gtail, np.array([kparam.return_king_params(energyval, 'gtail')]))
            SpE = np.append(SpE, np.array([kparam.rescale_factor(energyval)]))
            if ((self.eventtype == 5) or (self.eventtype == 0)):
                # Add in 3rd quartile if required
                kparam = PSFKC.PSF_king(psf_file, theta_norm=theta_norm_psf3, rescale_index=rescale_index_psf3,
                                        params_index=params_index_psf3)
                fcore = np.append(fcore, np.array([kparam.return_king_params(energyval, 'fcore')]))
                score = np.append(score, np.array([kparam.return_king_params(energyval, 'score')]))
                gcore = np.append(gcore, np.array([kparam.return_king_params(energyval, 'gcore')]))
                stail = np.append(stail, np.array([kparam.return_king_params(energyval, 'stail')]))
                gtail = np.append(gtail, np.array([kparam.return_king_params(energyval, 'gtail')]))
                SpE = np.append(SpE, np.array([kparam.rescale_factor(energyval)]))
                if self.eventtype == 0:
                    # Add in 4th quartile if required
                    kparam = PSFKC.PSF_king(psf_file, theta_norm=theta_norm_psf4, rescale_index=rescale_index_psf4,
                                            params_index=params_index_psf4)
                    fcore = np.append(fcore, np.array([kparam.return_king_params(energyval, 'fcore')]))
                    score = np.append(score, np.array([kparam.return_king_params(energyval, 'score')]))
                    gcore = np.append(gcore, np.array([kparam.return_king_params(energyval, 'gcore')]))
                    stail = np.append(stail, np.array([kparam.return_king_params(energyval, 'stail')]))
                    gtail = np.append(gtail, np.array([kparam.return_king_params(energyval, 'gtail')]))
                    SpE = np.append(SpE, np.array([kparam.rescale_factor(energyval)]))
        return fcore, score, gcore, stail, gtail, SpE

    def average_PSF(self):
        self.norm_dict = {
        comp: np.mean(self.template_dict[comp], axis=1) / np.sum(np.mean(self.template_dict[comp], axis=1)) for comp in
        self.template_dict.keys()}
        self.average_PSF_dict = {
        comp: np.sum([self.sigma_PSF_deg[i] * self.norm_dict[comp][i] for i in range(self.nEbins)]) for comp in
        self.template_dict.keys()}

    def load_psf_kings(self, eventclass=5, eventtype=3, **kwargs):
        if 'psf_energies' in kwargs.keys():
            self.psf_energies = kwargs['psf_energies']
        else:
            self.psf_energies = self.CTB_en_bins

        # self.fits_file_name = fits_file_path
        # Define the param and rescale indices for the various quartiles
        params_index_psf1 = 10
        rescale_index_psf1 = 11
        params_index_psf2 = 7
        rescale_index_psf2 = 8
        params_index_psf3 = 4
        rescale_index_psf3 = 5
        params_index_psf4 = 1
        rescale_index_psf4 = 2
        # Setup to load the correct PSF details depending on the dataset, and define appropriate theta_norm values, psf1 is bestpsf, psf4 is worst psf (but just a quartile each, so for Q1-3 need tocombine
        if eventclass == 2:
            psf_file_name = 'psf_P8R2_SOURCE_V6_PSF.fits'
            theta_norm_psf1 = [0.0000000, 9.7381019e-06, 0.0024811595, 0.022328802, 0.080147663, 0.17148392, 0.30634315,
                               0.41720551]
            theta_norm_psf2 = [0.0000000, 0.00013001938, 0.010239333, 0.048691643, 0.10790632, 0.18585539, 0.29140913,
                               0.35576811]
            theta_norm_psf3 = [0.0000000, 0.00074299273, 0.018672204, 0.062317201, 0.12894928, 0.20150553, 0.28339386,
                               0.30441893]
            theta_norm_psf4 = [4.8923139e-07, 0.011167475, 0.092594658, 0.15382001, 0.16862869, 0.17309118, 0.19837774,
                               0.20231968]
        elif eventclass == 5:
            psf_file_name = 'psf_P8R2_ULTRACLEANVETO_V6_PSF.fits'
            theta_norm_psf1 = [0.0000000, 9.5028121e-07, 0.00094418357, 0.015514370, 0.069725775, 0.16437751,
                               0.30868705, 0.44075016]
            theta_norm_psf2 = [0.0000000, 1.6070284e-05, 0.0048551576, 0.035358049, 0.091767466, 0.17568974, 0.29916159,
                               0.39315185]
            theta_norm_psf3 = [0.0000000, 0.00015569366, 0.010164870, 0.048955837, 0.11750811, 0.19840060, 0.29488095,
                               0.32993394]
            theta_norm_psf4 = [0.0000000, 0.0036816313, 0.062240006, 0.14027030, 0.17077023, 0.18329804, 0.21722594,
                               0.22251374]
        if 'fits_file_path' in kwargs.keys():
            self.fits_file_name = kwargs['fits_file_path']
        else:
            self.fits_file_name = self.psf_data_dir + psf_file_name

        self.f = fits.open(self.fits_file_name)
        # Now need to get the correct PSF for the appropriate quartile.
        # If anything other than best psf, need to combine quartiles.
        # Quartiles aren't exactly equal in size, but approximate as so.
        self.fcore_list = []
        self.score_list = []
        self.gcore_list = []
        self.stail_list = []
        self.gtail_list = []
        self.SpE_list = []

        self.append_to_kings_arrays(theta_norm_psf1, rescale_index_psf1, params_index_psf1)
        if ((eventtype == 4) or (eventtype == 5) or (eventtype == 0)):
            self.append_to_kings_arrays(theta_norm_psf2, rescale_index_psf2, params_index_psf2)
            if ((eventtype == 5) or (eventtype == 0)):
                self.append_to_kings_arrays(theta_norm_psf3, rescale_index_psf3, params_index_psf3)
                if eventtype == 0:
                    self.append_to_kings_arrays(theta_norm_psf4, rescale_index_psf4, params_index_psf4)
        self.average_PSF_kings()
        self.unfold_kings_lists()

    def append_to_kings_arrays(self, theta_norm_psf, rescale_index_psf, params_index_psf):
        PSFKC_inst = PSFKC.PSF_kings(self.f, theta_norm=theta_norm_psf, rescale_index=rescale_index_psf,
                                     params_index=params_index_psf)
        self.fcore_list += [PSFKC_inst.return_kings_params(self.psf_energies, 'fcore')]
        self.score_list += [PSFKC_inst.return_kings_params(self.psf_energies, 'score')]
        self.gcore_list += [PSFKC_inst.return_kings_params(self.psf_energies, 'gcore')]
        self.stail_list += [PSFKC_inst.return_kings_params(self.psf_energies, 'stail')]
        self.gtail_list += [PSFKC_inst.return_kings_params(self.psf_energies, 'gtail')]
        self.SpE_list += [PSFKC_inst.rescale_factor(self.psf_energies)]

    def average_PSF_kings(self):
        self.norm_dict = {
        comp: np.mean(self.template_dict[comp], axis=1) / np.sum(np.mean(self.template_dict[comp], axis=1)) for comp in
        self.template_dict.keys()}
        self.number_event_types = len(self.SpE_list)

        self.frac_list_unfolded_dict = {}
        self.frac_list_by_energy_dict = {}
        for comp in self.norm_dict.keys():
            self.frac_list_unfolded_dict[comp] = np.concatenate(
                np.array([self.norm_dict[comp] for i in range(self.number_event_types)]) / float(
                    self.number_event_types))
            self.frac_list_by_energy_dict[comp] = np.array([np.array(
                [self.norm_dict[comp][en] / float(self.number_event_types) for i in range(self.number_event_types)]) for
                                                            en in range(len(self.norm_dict[comp]))])

    def unfold_kings_lists(self):
        self.fcore_list_energy = np.transpose(self.fcore_list)
        self.score_list_energy = np.transpose(self.score_list)
        self.gcore_list_energy = np.transpose(self.gcore_list)
        self.stail_list_energy = np.transpose(self.stail_list)
        self.gtail_list_energy = np.transpose(self.gtail_list)
        self.SpE_list_energy = np.transpose(self.SpE_list)

        self.fcore_list = np.concatenate(self.fcore_list)
        self.score_list = np.concatenate(self.score_list)
        self.gcore_list = np.concatenate(self.gcore_list)
        self.stail_list = np.concatenate(self.stail_list)
        self.gtail_list = np.concatenate(self.gtail_list)
        self.SpE_list = np.concatenate(self.SpE_list)

    def make_f_ary_kings(self, ps_comp, energy_averaged=False, psf_save_tag='temp', num_f_bins=10, n_ps=10000,
                         n_pts_per_king=1000, **kwargs):
        self.psf_save_tag = psf_save_tag

        if 'psf_dir' in kwargs.keys():
            self.psf_dir = kwargs['psf_dir']
        else:
            self.psf_dir = self.work_dir + 'psf/'

        self.f_ary_list = []
        self.df_rho_div_f_ary_list = []
        if not energy_averaged:
            for en in range(self.CTB_en_bins[:-1]):
                mfak_inst = mfak.make_f_ary_kings(self.nside, self.psf_dir, self.fcore_list_energy[en],
                                                  self.score_list_energy[en], self.gcore_list_energy[en],
                                                  self.stail_list_energy[en], self.gtail_list_energy[en],
                                                  self.SpE_list_energy[en], self.frac_list_by_energy_dict[ps_comp][en],
                                                  save_tag=psf_save_tag, num_f_bins=num_f_bins, n_ps=n_ps,
                                                  n_pts_per_king=n_pts_per_king)

                self.f_ary_list += [mfak_inst.f_ary]
                self.df_rho_div_f_ary_list += [mfak_inst.df_rho_div_f_ary]

        else:
            mfak_inst = mfak.make_f_ary_kings(self.nside, self.psf_dir, self.fcore_list, self.score_list,
                                              self.gcore_list, self.stail_list, self.gtail_list, self.SpE_list,
                                              self.frac_list_unfolded_dict[ps_comp], save_tag=psf_save_tag,
                                              num_f_bins=num_f_bins, n_ps=n_ps, n_pts_per_king=n_pts_per_king)

            self.f_ary = mfak_inst.f_ary
            self.df_rho_div_f_ary = mfak_inst.df_rho_div_f_ary
            self.f_ary_list += [mfak_inst.f_ary]
            self.df_rho_div_f_ary_list += [mfak_inst.df_rho_div_f_ary]


class smooth_save_fits:
    def __init__(self, maps_dir, map_file_path, map_fits_number=0, fits_file_path='False', save_file_path='test.fits',
                 data_name='p8', nside=256, is_p8=True, eventclass=5, eventtype=3, highenergy = False, total_en_bins = []):
        self.nside = nside
        self.save_file_path = save_file_path
        self.maps_dir = maps_dir
        self.psf_data_dir = maps_dir + 'psf_data/'
        self.data_name = data_name
        self.map_file_path = map_file_path
        self.fits_file_path = fits_file_path
        self.eventclass = eventclass
        self.eventtype = eventtype


        self.total_en_bins = 10 ** np.linspace(np.log10(0.3), np.log10(300), 31)
        #High energy data
        self.highenergy = highenergy
        if self.highenergy:
            self.total_en_bins = 10 ** np.linspace(np.log10(50), np.log10(2000), 11)

        if not (total_en_bins == []):
            self.total_en_bins = total_en_bins

        self.maxenergy = len(self.total_en_bins) - 2 

        self.CTB_en_bins = self.total_en_bins# 10 ** np.linspace(np.log10(0.3), np.log10(300), 31)  # force energy bins over whole energy range

        self.load_psf(data_name=self.data_name, fits_file_path=self.fits_file_path, eventclass=self.eventclass,
                      eventtype=self.eventtype)

        self.maps_unsmoothed_fits = fits.open(self.map_file_path)[map_fits_number].data  # [map_fits_number]

        self.load_fits_map()
        self.load_psf(data_name=self.data_name, fits_file_path=self.fits_file_path, eventclass=self.eventclass,
                      eventtype=self.eventtype)
        self.set_initial_dirs(data=data_name)
        self.load_exposure()

        self.smoothed_maps = np.array([self.smooth_a_map(self.maps_unsmoothed[i], self.sigma_PSF_deg[i]) *
                                       self.CTB_exposure_maps[i] / np.mean(self.CTB_exposure_maps[i]) for i in
                                       range(len(self.CTB_en_bins) - 1)])

        self.save_fits_file()

    def load_exposure(self):
        self.CTB_exposure_maps = \
        self.CTB.get_CTB(self.CTB_dir, self.nside, 0, self.maxenergy, is_p8=self.is_p8, eventclass=self.eventclass,
                    eventtype=self.eventtype, newstyle=self.newstyle)[2]
        self.CTB_exposure_maps = np.array(list(self.CTB_exposure_maps) + [self.CTB_exposure_maps[-1]])

    def save_fits_file(self):
        hdu = fits.PrimaryHDU(self.smoothed_maps)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(self.save_file_path)

    def set_initial_dirs(self, data='old_CTBCORE'):
        self.is_p8 = False
        if data == 'old_CTBCORE':
            self.data_tag = 'p15_ultraclean_Q2_specbin'
        elif data == 'new_CTBCORE':
            self.data_tag = 'Apr3_Ultraclean_Q2'
        elif data == 'p8':
            self.is_p8 = True
            self.data_tag = 'PASS8_Jun15_UltracleanVeto_specbin'

        self.CTB_dir = self.maps_dir + self.data_tag + '/'

    def load_fits_map(self):
        if len(self.maps_unsmoothed_fits) == 30:
            self.maps_unsmoothed = self.maps_unsmoothed_fits
        elif len(np.shape(self.maps_unsmoothed_fits)) == 1:
            self.maps_unsmoothed = np.array([self.maps_unsmoothed_fits for i in range(30)])
        else:
            self.maps_unsmoothed = np.array([self.maps_unsmoothed_fits[0] for i in range(30)])

            # self.maps_smoothed = self.smooth_a_map(self.map_unsmoothed, sigma_PSF_deg)

    def smooth_a_map(self, the_map, sigma_PSF_deg):
        self.swp_inst = swp.smooth_gaussian_psf_quick(sigma_PSF_deg, the_map)
        # swp_inst.smooth_the_map()
        return self.swp_inst.the_smooth_map

    def load_psf(self, data_name='p8', fits_file_path='False', eventclass=5, eventtype=3):
        if fits_file_path != 'False':
            self.fits_file_name = fits_file_path
        elif data_name == 'p8':
            # Define the param and rescale indices for the various quartiles
            params_index_psf1 = 10
            rescale_index_psf1 = 11
            params_index_psf2 = 7
            rescale_index_psf2 = 8
            params_index_psf3 = 4
            rescale_index_psf3 = 5
            params_index_psf4 = 1
            rescale_index_psf4 = 2
            # Setup to load the correct PSF details depending on the dataset, and define appropriate theta_norm values, psf1 is bestpsf, psf4 is worst psf (but just a quartile each, so for Q1-3 need tocombine
            if eventclass == 2:
                psf_file_name = 'psf_P8R2_SOURCE_V6_PSF.fits'
                theta_norm_psf1 = [0.0000000, 9.7381019e-06, 0.0024811595, 0.022328802, 0.080147663, 0.17148392,
                                   0.30634315, 0.41720551]
                theta_norm_psf2 = [0.0000000, 0.00013001938, 0.010239333, 0.048691643, 0.10790632, 0.18585539,
                                   0.29140913, 0.35576811]
                theta_norm_psf3 = [0.0000000, 0.00074299273, 0.018672204, 0.062317201, 0.12894928, 0.20150553,
                                   0.28339386, 0.30441893]
                theta_norm_psf4 = [4.8923139e-07, 0.011167475, 0.092594658, 0.15382001, 0.16862869, 0.17309118,
                                   0.19837774, 0.20231968]
            elif eventclass == 5:
                psf_file_name = 'psf_P8R2_ULTRACLEANVETO_V6_PSF.fits'
                theta_norm_psf1 = [0.0000000, 9.5028121e-07, 0.00094418357, 0.015514370, 0.069725775, 0.16437751,
                                   0.30868705, 0.44075016]
                theta_norm_psf2 = [0.0000000, 1.6070284e-05, 0.0048551576, 0.035358049, 0.091767466, 0.17568974,
                                   0.29916159, 0.39315185]
                theta_norm_psf3 = [0.0000000, 0.00015569366, 0.010164870, 0.048955837, 0.11750811, 0.19840060,
                                   0.29488095, 0.32993394]
                theta_norm_psf4 = [0.0000000, 0.0036816313, 0.062240006, 0.14027030, 0.17077023, 0.18329804, 0.21722594,
                                   0.22251374]
            self.fits_file_name = self.psf_data_dir + psf_file_name
        if fits_file_path != 'False' or data_name == 'p8':
            self.f = fits.open(self.fits_file_name)
            # Now need to get the correct PSF for the appropriate quartile.
            # If anything other than best psf, need to combine quartiles.
            # Quartiles aren't exactly equal in size, but approximate as so.
            self.PSFC_inst = PSFC.PSF(self.f, theta_norm=theta_norm_psf1, rescale_index=rescale_index_psf1,
                                      params_index=params_index_psf1)
            calc_sigma_PSF_deg = np.array(self.PSFC_inst.return_sigma_gaussian(self.CTB_en_bins))
            if ((eventtype == 4) or (eventtype == 5) or (eventtype == 0)):
                self.PSFC_inst = PSFC.PSF(self.f, theta_norm=theta_norm_psf2, rescale_index=rescale_index_psf2,
                                          params_index=params_index_psf2)
                calc_sigma_load = np.array(self.PSFC_inst.return_sigma_gaussian(self.CTB_en_bins))
                calc_sigma_PSF_deg = (calc_sigma_PSF_deg + calc_sigma_load) / 2.
                if ((eventtype == 5) or (eventtype == 0)):
                    self.PSFC_inst = PSFC.PSF(self.f, theta_norm=theta_norm_psf3, rescale_index=rescale_index_psf3,
                                              params_index=params_index_psf3)
                    calc_sigma_load = np.array(self.PSFC_inst.return_sigma_gaussian(self.CTB_en_bins))
                    calc_sigma_PSF_deg = (2. * calc_sigma_PSF_deg + calc_sigma_load) / 3.
                    if eventtype == 0:
                        self.PSFC_inst = PSFC.PSF(self.f, theta_norm=theta_norm_psf4, rescale_index=rescale_index_psf4,
                                                  params_index=params_index_psf4)
                        calc_sigma_load = np.array(self.PSFC_inst.return_sigma_gaussian(self.CTB_en_bins))
                        calc_sigma_PSF_deg = (3. * calc_sigma_PSF_deg + calc_sigma_load) / 4.
            # Now take mean of array and extract first element, which will be the PSF
            self.sigma_PSF_deg = calc_sigma_PSF_deg
            self.average_PSF()
        else:
            pass
