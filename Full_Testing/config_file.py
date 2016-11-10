import numpy as np
import healpy as hp
from set_dirs import set_dirs as sd

class config(sd):
    """The configuration file """

    def __init__(self, nside=128, nexp=1, *args, **kwargs):

        self.kwargs = kwargs
        self.nside = nside
        self.nexp = nexp

        # setup the directories - TODO: come back and fix this up!
        # Remove run_tag
        sd.__init__(self, *args, **kwargs)

        # Initialise template dictionary and array
        self.templates = []
        self.templates_dict = {}

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
        
        self.mask_total = external_mask

    def add_template(self, template, label):
        """ Function to add a template to the template dictionary and array

            template: healpix map of the spatial template
            label: string used to identify the template in subsequent calls
        """
        
        self.templates_dict.update({label: template})
        self.templates.append(template)

    def return_masked_compressed(self, map_to_mask, expreg = False):
        """ Takes input healpix map, returns masked compressed version and if 
            expreg = True a set of such maps broken into exposure regions.
        """
        
        if not expreg:
            temp_masked_map = hp.ma(map_to_mask)
            temp_masked_map.mask = self.mask_total
            return temp_masked_map.compressed()
        else:
            temp_masked_map_list = []
            for i in range(self.nexp):
                temp_masked_map = hp.ma(map_to_mask)
                temp_masked_map.mask = np.logical_or(self.mask_total, 
                                       self.making_mask(self.pixel_list[i], 
                                       self.exposure_map))
                temp_masked_map_list.append(temp_masked_map.compressed())
            return temp_masked_map_list

    def compress_data_and_templates(self):
        """
        Compresses the data, which means masks the data with the given overall masks
        Does a version for poissonian and non-poissonian scans
        """

        # Compress data - this is used for a Poissonian scan
        temp_data = hp.ma(self.count_map)
        temp_data.mask = self.mask_total
        self.masked_compressed_data = temp_data.compressed()

        # Calculate mean exposure after masking - used in exposure corrections
        total_exposure_map_masked = hp.ma(self.exposure_map)
        total_exposure_map_masked.mask = self.mask_total
        self.exposure_mean = np.mean(total_exposure_map_masked.compressed())

        # Divide the map into exposure regions
        temp_exposure = hp.ma(self.exposure_map)
        temp_exposure.mask = self.mask_total
        exposure_map_masked_compressed = temp_exposure.compressed()
        
        self.pixel_list_compressed, self.exposure_list_compressed = \
        self.dividing_exposure(exposure_map_masked_compressed)
        
        self.pixel_list, self.exposure_list = \
        self.dividing_exposure(self.exposure_map)
        self.exposure_means_list = [np.mean(self.exposure_list_compressed[i]) \
                            for i in range(len(self.exposure_list_compressed))]
        
        self.masked_compressed_data_expreg = []
        for i in range(self.nexp):
            temp_data_expreg = hp.ma(self.count_map)
            temp_data_expreg.mask = np.logical_or(self.mask_total, \
                        self.making_mask(self.pixel_list[i], self.exposure_map))
            self.masked_compressed_data_expreg.append(temp_data_expreg.compressed())

        # Create a nested dictionary of different versions of the templates
        the_dict = self.templates_dict
        keys = self.templates_dict.keys()
        self.templates_dict_nested = {
            key: {'template': 
                  the_dict[key], 
                  'template_masked_compressed': 
                  self.return_masked_compressed(the_dict[key]),
                  'template_masked_compressed_expreg': 
                  self.return_masked_compressed(the_dict[key], expreg=True)}
            for key in keys}

    def dividing_exposure(self, exposure_array):
        """
        This function divides the sky into nexp "almost" equal parts and returns an array of lists of pixels
        :param exposure_array: This should be the exposure array (before masking). It's called here exposure_map
        :param nexp: Number of divisions required
        :return: [ [list of pixels], [list of pixels]... ] with nexp list of pixels of almost equal exposure
        """
        pix_array = np.arange(len(exposure_array))
        new_array = np.array([[pix_array[i], exposure_array[i]] for i in range(len(exposure_array))])
        col = 1
        array_sorted = new_array[np.argsort(new_array[:, col])]
        array_split = np.array_split(array_sorted, self.nexp)
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
