import numpy as np
import healpy as hp
from set_dirs import set_dirs as sd

class config(sd):
    """The configuration file """

    def __init__(self, nside=128, nexp=1, use_exposure=False, *args, **kwargs):

        self.kwargs = kwargs
        self.nside = nside

        # setup the directories - TODO: come back and fix this up!
        sd.__init__(self, *args, **kwargs)

        # Initialise template dictionary and array
        self.templates = []
        self.templates_dict = {}

        # TODO: remove any mention of this from other files
        self.use_fixed_template = False

        # TODO: Always use exposure regions, remove below once setup
        self.nexp = nexp
        self.use_exposure = use_exposure
        if self.nexp > 1:
            self.use_exposure = True

    # TODO: move below two functions to where they are needed, not this file
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


    def load_data(self, count_map, exposure_map):
        """ Function to input analysis data
            
            count_map: a healpix map of counts (integers)
            exposure map: a healpix map of the exposure
        """
        
        self.count_map = count_map
        self.exposure_map = exposure_map

    def load_mask(self, external_mask):
        """ Function to input analysis mask
        """
        
        # Append the mask to self
        self.mask_total = external_mask

        # Calculate mean exposure after masking - used in exposure corrections
        total_exposure_map_masked = hp.ma(self.exposure_map)
        total_exposure_map_masked.mask = self.mask_total
        self.exposure_mean = np.mean(total_exposure_map_masked.compressed())

    def add_template(self, template, label):
        """ Function to add a template to the template dictionary and array

            template: healpix map of the spatial template
            label: string used to identify the template in subsequent calls
        """
        
        self.templates_dict.update({label: template})
        self.templates.append(template)







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
        pixel_list, exposure_list = self.dividing_exposure(self.exposure_map, self.nexp)
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

        masked_arrays_list = []
        for i in range(self.nexp):
            if flux:
                masked_arrays = [hp.ma((array/flux_map)*np.mean(flux_map)) for array, flux_map in map(None, arrays, [self.exposure_map])]
            else:
                masked_arrays = [hp.ma(array) for array in arrays]
            for array in masked_arrays:
                array.mask = np.logical_or(self.mask_total, self.making_mask(self.pixel_list[i], self.exposure_map))
            masked_arrays_list.append(masked_arrays)
        if not compressed:
            return masked_arrays_list
        else:
            masked_arrays_compressed = np.array([[arr.compressed() for arr in array] for array in masked_arrays_list])
        if summed:
            masked_arrays_compressed = np.sum(masked_arrays_compressed, axis=1)
        return masked_arrays_compressed

    def compress_data(self):
        """
        Compresses the data, which means masks the data with the given overall maskS
        :return: self.CTB_masked_compressed_summed where it is summed over energy bins? Not sure...
        """
        # TODO: Streamline this to have the same function for templates and data
        # TODO: Make energy dependent by changing exposure_map[0]

        self.CTB_masked = [hp.ma(self.count_map)]
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
        temp = hp.ma(self.exposure_map)
        temp.mask = self.mask_total
        self.exposure_map_masked_compressed = [temp.compressed()]
        self.exposure_map_masked_compressed_mean = np.mean(self.exposure_map_masked_compressed,axis=0)
        self.pixel_list_compressed, self.exposure_list_compressed = self.dividing_exposure(self.exposure_map_masked_compressed_mean, self.nexp)
        self.pixel_list, self.exposure_list = self.dividing_exposure(self.exposure_map, self.nexp)

        self.exposure_means_list = [np.mean(self.exposure_list_compressed[i]) for i in range(len(self.exposure_list_compressed))]
        self.CTB_masked_compressed_exposure = []
        for i in range(self.nexp):
            for the_map in self.CTB_masked:
                the_map.mask = np.logical_or(self.mask_total, self.making_mask(self.pixel_list[i], self.exposure_map))
            self.CTB_masked_compressed_exposure.append([the_map.compressed() for the_map in self.CTB_masked])

        self.CTB_masked_compressed_summed_exposure = np.sum(self.CTB_masked_compressed_exposure, axis=1)

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

