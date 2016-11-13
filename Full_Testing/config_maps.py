#------------------------------------------------------------------------------
# Configure Maps
#------------------------------------------------------------------------------
##
## Configure the maps - data, exposure and templates - required to perform the
## scan.
##
## These maps are input by the user as well as the mask for the analysis.
## The code contains routines to convert the maps to masked arrays, and if
## required divide it into exposure regions
##
#------------------------------------------------------------------------------

import numpy as np
import healpy as hp
from set_dirs import set_dirs as sd

class config_maps(sd):
    def __init__(self, nside=128, nexp=1, *args, **kwargs):
        #TODO if possible get rid of args/kwargs

        self.kwargs = kwargs
        self.nside = nside
        self.nexp = nexp

        # setup the directories - TODO: come back and fix this up! i.e. remove run tag
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

    def compress_data_and_templates(self):
        """ Compress data and templates
            Before calling must have loaded data, exposure, templates and mask
        """

        # Compress data - this is used for a Poissonian scan
        temp_data = hp.ma(self.count_map)
        temp_data.mask = self.mask_total
        self.masked_compressed_data = temp_data.compressed()

        # Divide the map into exposure regions
        self.divide_exposure()

        # Create an exposure region corrected version of the data
        self.masked_compressed_data_expreg = []
        for i in range(self.nexp):
            temp_data_expreg = hp.ma(self.count_map)
            temp_data_expreg.mask = self.expreg_mask[i] 
            self.masked_compressed_data_expreg.append(temp_data_expreg.\
                                                          compressed())

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

    def divide_exposure(self):
        """ Divide the ROI into nexp different regions of similar exposure
        """
        
        # Determine the pixels of the exposure regions
        pix_array = np.where(self.mask_total == False)[0]
        exp_array = np.array([[pix_array[i], self.exposure_map[pix_array[i]]] \
                                             for i in range(len(pix_array))])
        array_sorted = exp_array[np.argsort(exp_array[:, 1])]
        
        # Stop code from dividing by more pixels than there are in the ROI
        if (self.nexp > len(array_sorted)):
            print "nexp cannot be larger than the number of pixels in the ROI"
            print "Setting nexp = ROI size"
            self.nexp = len(array_sorted)

        # Convert from list of exreg pixels to masks
        array_split = np.array_split(array_sorted, self.nexp)
        expreg_array = [[array_split[i][j][0] \
                       for j in range(len(array_split[i]))] \
                       for i in range(len(array_split))]

        temp_expreg_mask = []
        for i in range(self.nexp):
            temp_mask = np.logical_not(np.zeros(hp.nside2npix(self.nside)))
            for j in range(len(expreg_array[i])):
                temp_mask[expreg_array[i][j]] = False
            temp_expreg_mask.append(temp_mask)

        self.expreg_mask = temp_expreg_mask

        # Store the total and region by region mean exposure
        expreg_values = [[array_split[i][j][1] \
                        for j in range(len(array_split[i]))] \
                        for i in range(len(array_split))]
        self.exposure_means_list = [np.mean(expreg_values[i]) \
                                   for i in range(self.nexp)]
        self.exposure_mean = np.mean(self.exposure_means_list)

    def return_masked_compressed(self, map_to_mask, expreg = False):
        """ Take input healpix map, return masked compressed version and if
            expreg = True a set of such maps broken into exposure regions
        """

        if not expreg:
            temp_masked_map = hp.ma(map_to_mask)
            temp_masked_map.mask = self.mask_total
            return temp_masked_map.compressed()
        else:
            temp_masked_map_list = []
            for i in range(self.nexp):
                temp_masked_map = hp.ma(map_to_mask)
                temp_masked_map.mask = self.expreg_mask[i]
                temp_masked_map_list.append(temp_masked_map.compressed())
            return temp_masked_map_list
