import sys, os

# current_dir = os.getcwd()
# change_path = ".."
# os.chdir(change_path)

# This file stores input
import healpy as hp
import numpy as np
import numpy.ma as ma
from config.set_dirs import set_dirs as sd
import pulsars.masks as masks
import pulsars.CTB as CTB
import config.make_templates as mt
import config.compute_PSF as CPSF

from config.make_NFW import make_NFW
from config.load_ps_mask import load_ps_mask as lpm

from astropy.io import fits

import logging
import pickle


class config(sd):
    """The configuration file """

    # def __init__(self,tag='test',CTB_en_min = 8, CTB_en_max = 16, nside = 128,data_name = 'new_CTBCORE',make_plots_data_tag = True,**kwargs):#,**kwargs):
    def __init__(self, nside=128, nexp=1, use_exposure=False, *args, **kwargs):

        self.kwargs = kwargs

        # for saving and loading
        # self.tag = tag

        # nside
        self.nside = nside
        self.npix = hp.nside2npix(self.nside)

        # Parameters for the CTBCORE data
        # self.data_name = data_name

        # self.CTB_en_min = CTB_en_min;
        # self.CTB_en_max = CTB_en_max;
        # self.nEbins = self.CTB_en_max - self.CTB_en_min

        # setup the directories
        sd.__init__(self, *args, **kwargs)  #


        # #Now load the CTBCORE data
        # if 'CTB_count_maps' in self.kwargs.keys():
        #     self.load_external_data()
        # else:
        #     self.load_CTBCORE_data()
        # self.load_data_ps_model_and_maps()

        # templates will be added to self.templates
        self.templates = []
        self.comp_array = []
        self.templates_dict = {}
        self.templates_dict_nested = {}

        self.spect = 'False'
        self.norm_arrays = 'False'
        self.use_fixed_template = False

        self.additional_template_dict = {}

        self.fixed_ps_model_dict = {}
        self.nexp = nexp
        self.use_exposure = use_exposure
        if self.nexp > 1:
            self.use_exposure = True

    def make_plane_mask(self, band_mask_range):
        plane_mask_array = masks.mask_lat_band(90 + band_mask_range[0], 90 + band_mask_range[1], self.nside)
        return plane_mask_array

    def make_long_mask(self, lmin, lmax):
        long_mask_array = masks.mask_not_long_band(lmin, lmax, self.nside)
        return long_mask_array

    def make_lat_mask(self, bmin, bmax):
        lat_mask_array = masks.mask_not_lat_band(90 + bmin, 90 + bmax, self.nside)
        return lat_mask_array

    def make_ring_mask(self, inner, outer, lat, lng):
        ring_mask = np.logical_not(masks.mask_ring(inner, outer, lat, lng, self.nside))
        return ring_mask

    def make_mask_total(self, plane_mask = True, band_mask_range = [-30,30], lcut=False, lmin=-20, lmax=20,bcut=False, bmin=-20, bmax=20, mask_ring=False, inner=0, outer=30,lat= 90, lng = 0, ps_mask_array = 'False', custom_mask=False,custom_mask_arr=np.array([0,1])):
        """ Combine various masks into a single boolean array where True pixels are masked """
        # First add all variables to self as use some of them later
        self.plane_mask = plane_mask
        self.band_mask_range = band_mask_range
        self.lcut = lcut
        self.lmin = lmin
        self.lmax = lmax
        self.bcut = bcut
        self.bmin = bmin
        self.bmax = bmax
        self.mask_ring = mask_ring
        self.inner = inner
        self.outer = outer
        self.lat = lat
        self.lng = lng
        self.ps_mask_array = ps_mask_array
        self.custom_mask = custom_mask
        self.custom_mask_arr = custom_mask_arr

        # Now compute all the masks we might need
        self.plane_mask_ar = self.make_plane_mask(self.band_mask_range)
        self.long_mask_ar = self.make_long_mask(self.lmin, self.lmax)
        self.lat_mask_ar = self.make_lat_mask(self.bmin, self.bmax)
        self.ring_mask = self.make_ring_mask(self.inner, self.outer, self.lat, self.lng)
        # For custom masks, the custom_mask_arr must be an array of pixels in the nside
        # array you want to mask
        toadd_custom_mask = np.zeros(self.nside**2*12, dtype=bool)
        toadd_custom_mask[custom_mask_arr]=True

        # Next initialise a healpix array of false values as only true values are masked
        mask_array = np.zeros(self.nside**2*12, dtype=bool)

        if self.plane_mask:
            mask_array += self.plane_mask_ar

        if self.lcut:
            mask_array += self.long_mask_ar

        if self.bcut:
            mask_array += self.lat_mask_ar

        if self.mask_ring:
            mask_array += self.ring_mask

        if self.custom_mask:
            mask_array += toadd_custom_mask

        self.mask_geom_total = mask_array

        # Finally add in a point source mask if this was included
        if self.ps_mask_array == 'False':
            self.mask_total = self.mask_geom_total
            self.mask_total_array = np.array([self.mask_geom_total for en in self.CTB_en_bins])
        else:
            self.mask_total = np.sum(self.ps_mask_array, axis=0, dtype=bool) + self.mask_geom_total
            self.mask_total_array = np.array(self.ps_mask_array + self.mask_geom_total, dtype=bool)

        self.recompute_mean_exposure()

    def recompute_mean_exposure(self):
        """
        :return: computes the mean of the exposure over the ROI, as well as the total_exposure_map_masked
        """
        self.total_exposure_map_masked = hp.ma(self.total_exposure_map)
        self.total_exposure_map_masked.mask = self.mask_total
        self.exposure_mean = np.mean(self.total_exposure_map_masked.compressed())


            # Below is the original version of make_mask_total, replaced by the one above by Nick Rodd on 15 Sept 2015
            #    def make_mask_total(self,band_mask_range = [-30,30],mask_ring = False,lat= 90,lng = 0,inner = 0, outer= 30,ps_mask_array = 'False'):
            #        self.band_mask_range = band_mask_range
            #        self.band_mask = self.make_band_mask(self.band_mask_range)
            #
            #        self.mask_ring = mask_ring
            #        self.lat = lat
            #        self.lng = lng
            #        self.inner = inner
            #        self.outer = outer
            #
            #        if not self.mask_ring:
            #            self.mask_geom_total = self.band_mask #fix this!
            #        else:
            #            self.ring_mask = self.make_ring_mask(self.inner,self.outer,self.lat,self.lng)
            #            self.mask_geom_total = self.ring_mask + self.band_mask
            #
            #        # self.ps_mask_type = ps_mask_type
            #        self.ps_mask_array = ps_mask_array
            #        if self.ps_mask_array == 'False':
            #            self.mask_total = self.mask_geom_total
            #            self.mask_total_array = np.array([self.mask_geom_total for en in self.CTB_en_bins])
            #        else:
            #            self.mask_total = np.sum(self.ps_mask_array,axis=0,dtype=bool)+self.mask_geom_total
            #            self.mask_total_array = np.array(self.ps_mask_array + self.mask_geom_total,dtype=bool)

    def load_external_data(self, CTB_en_bins, CTB_count_maps, CTB_exposure_maps):
        self.CTB_en_bins = CTB_en_bins
        self.CTB_en_min = self.CTB_en_bins[0]
        self.CTB_en_max = self.CTB_en_bins[-1]
        self.nEbins = len(self.CTB_en_bins)
        self.CTB_exposure_maps = CTB_exposure_maps
        self.CTB_count_maps = CTB_count_maps

        self.total_exposure_map = np.mean(self.CTB_exposure_maps, axis=0)  # This is a mean over the energies
        self.NPIX_total = np.size(self.total_exposure_map)

        self.exposure_mean = np.mean(self.total_exposure_map)  # TODO: Incorporate energy dependence

    def rebin_external_data(self, nEbins,**kwargs):
        # TODO: Figure out if this is for energy dependent. If so, refactor with above function
        bin_len = (len(self.CTB_en_bins) - 1) / nEbins
        self.CTB_en_bins = np.array([self.CTB_en_bins[bin_len * i] for i in range(nEbins + 1)])
        self.CTB_en_min = self.CTB_en_bins[0]
        self.CTB_en_max = self.CTB_en_bins[-1]
        self.nEbins = len(self.CTB_en_bins)

        if 'spect_coeffs' in kwargs.keys():
            spect_coeffs = kwargs['spect_coeffs']
            '''
            spect_coeffs should be of the form [[a_1,a_2, ...],[b_1,b_2,...]],
            where \sum a_i = 1 and \sum b_i = 1, etc.
            '''
            self.CTB_exposure_maps = np.array(
                [np.sum(np.array([self.CTB_exposure_maps[bin_len * i:bin_len * (i + 1)][j]*np.array(spect_coeffs[i])[j] for j in range(len(spect_coeffs[i]))]), axis=0) for i in range(nEbins)])
        else:
            self.CTB_exposure_maps = np.array(
                [np.mean(self.CTB_exposure_maps[bin_len * i:bin_len * (i + 1)], axis=0) for i in range(nEbins)])
        self.CTB_count_maps = np.array(
            [np.sum(self.CTB_count_maps[bin_len * i:bin_len * (i + 1)], axis=0) for i in range(nEbins)])

        self.total_exposure_map = np.mean(self.CTB_exposure_maps, axis=0)  # This is a mean over the energies
        self.NPIX_total = np.size(self.total_exposure_map)

        for comp in self.templates_dict.keys():
            self.templates_dict[comp] = np.array(
                [np.sum(self.templates_dict[comp][bin_len * i:bin_len * (i + 1)], axis=0) for i in range(nEbins)])

        if self.use_fixed_template:
            for comp in self.fixed_template_dict.keys():
                self.fixed_template_dict[comp] = np.array(
                    [np.sum(self.fixed_template_dict[comp][bin_len * i:bin_len * (i + 1)], axis=0) for i in
                     range(nEbins)])

        self.mask_total_array = np.vectorize(bool)(
            np.array([np.sum(self.mask_total_array[bin_len * i:bin_len * (i + 1)], axis=0) for i in range(nEbins)]))

    def add_fixed_templates(self, template_dict):
        self.use_fixed_template = True
        self.fixed_template_dict = template_dict

    def add_new_template(self, input_template_dict, reset_templates=False):
        if reset_templates:
            del self.comp_array
            del self.templates
            self.comp_array = []
            self.templates = []
            self.templates_dict = {}

        # templates_correct_nside = np.array(  [hp.ud_grade(temp,self.nside, power=-2) for temp in templates]    )

        self.templates_dict.update(input_template_dict)

        [self.comp_array.append(comp) for comp in input_template_dict.keys()]
        [self.templates.append(temp) for temp in input_template_dict.values()]

    def load_template_from_file(self, filename):  # should already be an exposure-corrected counts map
        fits_file_list = fits.open(self.template_dir + filename)[0].data
        return np.array([hp.ud_grade(fits_file, self.nside, power=-2) for fits_file in fits_file_list])[
               self.CTB_en_min:self.CTB_en_max]

    def load_fixed_ps_model(self, filename, comp):
        with open(self.data_ps_dir + filename, 'rb') as f:
            the_list = pickle.load(f)
        self.fixed_ps_model_dict[comp] = {'medians': the_list[0], 'template': the_list[1]}

    def compress_templates(self):
        self.templates_masked = [[hp.ma(temp_2) for temp_2 in temp_1] for temp_1 in self.templates]
        for temp_masked in self.templates_masked:
            for temp_2 in temp_masked:
                temp_2.mask = self.mask_total
        self.templates_masked_compressed = [[temp_2.compressed() for temp_2 in temp_1] for temp_1 in
                                            self.templates_masked]

        the_dict = self.templates_dict
        keys = self.templates_dict.keys()

        #divide into exposure regions
        pixel_list, exposure_list = self.dividing_exposure(np.mean(self.CTB_exposure_maps, axis=0), self.nexp)
        self.pixel_list = pixel_list
        self.exposure_list = exposure_list

        self.templates_dict_nested = {
            key: {'templates': the_dict[key], 'templates_masked': self.return_masked(the_dict[key]),
                  'templates_masked_compressed': self.return_masked(the_dict[key], compressed=True),
                  'summed_templates': self.return_masked(the_dict[key], compressed=True, summed=True),
                  'summed_templates_not_compressed': np.sum(the_dict[key], axis=0),
                  'templates_masked_exposure': self.return_masked_exposure(the_dict[key]),
                  'flux_templates_masked_exposure': self.return_masked_exposure(the_dict[key],flux=True),
                  'templates_masked_compressed_exposure': self.return_masked_exposure(the_dict[key], compressed=True),
                  'flux_templates_masked_compressed_exposure': self.return_masked_exposure(the_dict[key],compressed=True,flux=True),
                  'summed_templates_exposure': self.return_masked_exposure(the_dict[key], compressed=True, summed=True),
                  'flux_summed_templates_exposure': self.return_masked_exposure(the_dict[key], compressed=True, summed=True,flux=True)}
            for key in keys}

        if self.use_fixed_template:
            the_dict_fixed = self.fixed_template_dict
            keys_fixed = self.fixed_template_dict.keys()
            new_template = np.zeros(np.shape(the_dict_fixed[keys_fixed[0]]))
            for key in keys_fixed:
                new_template += the_dict_fixed[key]
            self.fixed_template_dict_nested = {'templates': new_template,
                                               'templates_masked': self.return_masked(new_template),
                                               'templates_masked_compressed': self.return_masked(new_template, compressed=True),
                                               'summed_templates': self.return_masked(new_template, compressed=True, summed=True),
                                               'summed_templates_not_compressed': np.sum(new_template, axis=0),
                                               'template_masked_exposure': self.return_masked_exposure(new_template),
                                               'templates_masked_compressed_exposure': self.return_masked_exposure(new_template, compressed=True),
                                               'summed_templates_exposure': self.return_masked_exposure(new_template, compressed=True,summed=True),
                                               'flux_templates_masked_compressed_exposure': self.return_masked_exposure(new_template, compressed=True, flux=True),
                                               'flux_summed_templates_exposure': self.return_masked_exposure(new_template,compressed=True,summed=True, flux=True),
                                               'flux_templates_masked_exposure': self.return_masked_exposure(new_template, flux=True)}

    def return_masked(self, arrays, compressed=False, summed=False):
        """
        Masking arrays without the number of exposure regions.
        :param arrays: array to be masked (generally templates)
        :param compressed: whether or not to compress the array
        :param summed: whether or not to sum that array
        :return: returns a masked array of the same length, then may be compressed and/or summed
        """
        masked_arrays = [hp.ma(array) for array in arrays]
        for array in masked_arrays:
            array.mask = self.mask_total
        if not compressed:
            return masked_arrays
        else:
            masked_arrays_compressed = np.array([array.compressed() for array in masked_arrays])
        if summed:
            masked_arrays_compressed = np.sum(masked_arrays_compressed, axis=0)
        return masked_arrays_compressed

    def return_masked_exposure(self, arrays, compressed=False, summed=False, flux=False):
        """
        Takes an array and masks with the usual mask + exposure mask.
        :param arrays: Array to mask
        :param compressed: compresses the masked array
        :param summed:  sums the elements of the masked array
        :return: returns a list of self.nexp elements, each masked by total mask+ exposure, and compressed/summed
        """
        # TODO: Make this energy dependent. Change self.CTB_exposure_maps[0]

        # pixel_list, exposure_list = self.dividing_exposure(np.mean(self.CTB_exposure_maps, axis=0), self.nexp)
        # self.pixel_list = pixel_list
        # self.exposure_list = exposure_list
        masked_arrays_list = []
        #self.mask_array = []
        for i in range(self.nexp):
            if flux:
                masked_arrays = [hp.ma((array/flux_map)*np.mean(flux_map)) for array, flux_map in map(None, arrays, self.CTB_exposure_maps)]
            else:
                masked_arrays = [hp.ma(array) for array in arrays]
            for array in masked_arrays:
                #self.mask_array += [np.logical_or(self.mask_total, self.making_mask(pixel_list[i], self.CTB_exposure_maps[0]))]
                array.mask = np.logical_or(self.mask_total, self.making_mask(self.pixel_list[i], self.CTB_exposure_maps[0]))
            masked_arrays_list.append(masked_arrays)
        #self.masked_arrays_list = masked_arrays_list
        if not compressed:
            return masked_arrays_list
        else:
            masked_arrays_compressed = np.array([[arr.compressed() for arr in array] for array in masked_arrays_list])
        if summed:
            masked_arrays_compressed = np.sum(masked_arrays_compressed, axis=1)
        return masked_arrays_compressed

    # def compress_fixed_templates(self):
    #     self.fixed_templates_masked = [hp.ma(temp_1) for temp_1 in self.fixed_templates]
    #     for temp_masked in self.fixed_templates_masked:
    #         temp_masked.mask = self.mask_total
    #     self.fixed_templates_masked_compressed = [temp_2.compressed() for temp_2 in self.fixed_templates_masked]

    # def sum_fixed_templates_over_energy_range(self): #flag: old
    #     self.summed_fixed_templates = np.sum(self.fixed_templates_masked_compressed,axis=0)

    def compress_data(self):
        """
        Compresses the data, which means masks the data with the given overall maskS
        :return: self.CTB_masked_compressed_summed where it is summed over energy bins? Not sure...
        """
        # TODO: Streamline this to have the same function for templates and data
        # TODO: Make energy dependent by changing CTB_exposure_maps[0]

        self.CTB_masked = [hp.ma(the_data) for the_data in self.CTB_count_maps]
        for the_map in self.CTB_masked:
            the_map.mask = self.mask_total
        self.CTB_masked_compressed = [the_map.compressed() for the_map in self.CTB_masked]

        self.CTB_masked_compressed_summed = np.sum(self.CTB_masked_compressed, axis=0)

    def compress_data_exposure(self):
        """
        Compresses the data and break it by exposure, which means masks the data with the given overall mask and
        exposure map masks
        :return: self.CTB_masked_compressed_summed where it is summed over ?
        """
        # Reverting back here to find bug
        self.exposure_pixels_compressed()  # This compresses the exposure and gives out the means
        # self.CTB_masked = [hp.ma(the_data) for the_data in self.CTB_count_maps]
        self.CTB_masked_compressed_exposure = []
        for i in range(self.nexp):
            for the_map in self.CTB_masked:
                the_map.mask = np.logical_or(self.mask_total, self.making_mask(self.pixel_list[i], self.CTB_exposure_maps[0]))
            self.CTB_masked_compressed_exposure.append([the_map.compressed() for the_map in self.CTB_masked])

        self.CTB_masked_compressed_summed_exposure = np.sum(self.CTB_masked_compressed_exposure, axis=1)

    @staticmethod
    def convert_log_list(the_list, is_log):
        new_list = []
        for i in range(len(the_list)):
            if is_log[i]:
                new_list.append(10 ** the_list[i])
            else:
                new_list.append(the_list[i])
        return new_list

    def flatten(self, the_list):
        the_list = list(the_list)
        return [item for sublist in the_list for item in list(sublist)]

    def make_logging(self, log_dir):
        sd.make_logging(log_dir)

    def dividing_exposure(self, exposure_array, nexp):
        """
        This function divides the sky into nexp "almost" equal parts and returns an array of lists of pixels
        :param exposure_array: This should be the exposure array (before masking). It's called here CTB_exposure_map
        :param nexp: Number of divisions required
        :return: [ [list of pixels], [list of pixels]... ] with nexp list of pixels of almost equal exposure
        """
        pix_array = np.arange(len(exposure_array))
        new_array = np.array([[pix_array[i], exposure_array[i]] for i in range(len(exposure_array))])
        col = 1
        array_sorted = new_array[np.argsort(new_array[:, col])]
        array_split = np.array_split(array_sorted, nexp)
        just_pixels = [[array_split[i][j][0] for j in range(len(array_split[i]))] for i in range(len(array_split))]
        just_exposures = [[array_split[i][j][1] for j in range(len(array_split[i]))] for i in range(len(array_split))]
        return just_pixels, just_exposures

    @staticmethod
    def making_mask(mask_array, data_masked):
        """
        Starts from a list of pixels and masks everything but those pixels
        :param mask_array: This is the list of pixels
        :param data_masked: This is the dataset to be masked
        :return:
        """
        new_data_masked = np.zeros(len(data_masked))
        for i in range(len(mask_array)):
            new_data_masked[mask_array[i]] = 1
        return np.logical_not(new_data_masked)

    def masked_data(self, array_to_mask):
        """
        Masks an array by dividing the sky into nexp equal parts in exposure. An example of using this function
        is to plot the different exposure regions, or divide the data and templates
        :param array_to_mask: the array that we want to mask. This could be a data set
        :return: ith masked array (which could be data or templates)
        """
        # TODO: Works only for one energy bin for now since I am doing self.CTB_exposure_maps[0]. Fix that
        pixel_list, exposure_list = self.dividing_exposure(self.CTB_exposure_maps[0], self.nexp)
        masked_data_array = []
        for i in range(self.nexp):
            mask = pixel_list[i]  # i is for the splitting group
            data_mask = np.zeros(len(array_to_mask), dtype=bool)
            data_mask = self.making_mask(mask, data_mask)
            masked_array = ma.array(array_to_mask, mask=np.logical_not(data_mask))
            masked_data_array.append(masked_array)
        return masked_data_array

    def mask_compress_exposure(self):
        self.CTB_exposure_maps_masked_compressed = []
        for the_map in self.CTB_exposure_maps:
            temp = hp.ma(the_map)
            temp.mask = self.mask_total
            self.CTB_exposure_maps_masked_compressed += [temp.compressed()]
        self.CTB_exposure_maps_masked_compressed_mean = np.mean(self.CTB_exposure_maps_masked_compressed,axis=0)



    def exposure_pixels_compressed(self):
        """
        Takes the exposure pixel list and compresses it following the overall mask
        :return: self.pixel_list_compressed
        """

        self.mask_compress_exposure()
        self.pixel_list_compressed, self.exposure_list_compressed = self.dividing_exposure(self.CTB_exposure_maps_masked_compressed_mean, self.nexp)
        self.pixel_list, self.exposure_list = self.dividing_exposure(np.mean(self.CTB_exposure_maps, axis=0), self.nexp)
        #self.pixel_list_compressed = []
        #self.exposure_means_list = []

        self.exposure_means_list = [np.mean(self.exposure_list_compressed[i]) for i in range(len(self.exposure_list_compressed))]
        # for exp_group in range(len(self.pixel_list)):
        #     pixels_compressed = []
        #     exposure_list_compressed = []
        #     for i in range(exp_group):
        #         if not self.mask_total[i]:
        #             pixels_compressed.append(i)
        #             exposure_list_compressed.append(self.exposure_list[exp_group][i])
        #     self.pixel_list_compressed.append(pixels_compressed)
        #     self.exposure_means_list.append(np.mean(exposure_list_compressed))
