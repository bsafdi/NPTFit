###############################################################################
# config_maps.py
###############################################################################
#
# Configure the maps - data, exposure, and templates - required to perform the
# scan.
#
# These maps are input by the user as well as the mask for the analysis.
# The code contains routines to convert the maps to masked arrays, and if
# required divide it into exposure regions.
#
###############################################################################

from __future__ import print_function
from __future__ import absolute_import

import numpy as np
import numpy.ma as ma

from .set_dirs import SetDirs  # Module for creating directories


class ConfigMaps(SetDirs):
    def __init__(self, tag='Untagged', work_dir=None):
        """
            :param tag: Label associated with the scan
            :param work_dir: Location where all the output from the
            scan is stored
        """

        # Setup the directories needed by the code
        SetDirs.__init__(self, tag=tag, work_dir=work_dir)

        # Initialise data, exposure and mask
        self.count_map = []
        self.exposure_map = []
        self.mask_total = []

        # Initialise template dictionary and array
        self.templates = []
        self.templates_dict = {}

    def load_data(self, count_map, exposure_map):
        """ Function to input analysis data

            :param count_map: A map of counts (integers)
            :param exposure_map: A map of the exposure
        """

        self.count_map = count_map
        self.exposure_map = exposure_map

    def load_mask(self, external_mask):
        """ Function to input analysis mask

            :param external_mask: The analysis mask
        """

        self.mask_total = external_mask

    def add_template(self, template, label, units='counts'):
        """ Function to add a template to the template dictionary and array

            Note templates should be exposure corrected, so that they model
            the counts rather than flux, before being added

            :param template: Map of the spatial template
            :param label: String used to identify the template in
            subsequent calls
            :param units: Units of provided template. By default
            'counts': Template in photon counts
            'flux': Template in fluxes with units ph/cm^2/s

            .. note:: The exposure must be provided prior to adding a flux
            template. This is then multipleid by the exposure.
        """

        if units == 'flux':
            assert (len(self.exposure_map) != 0), \
                "Must provide exposure map before adding a flux template"
            assert (len(self.exposure_map) == len(template)), \
                "Template must be the same shape as the exposure map"
            template *= self.exposure_map
        self.templates_dict.update({label: template})
        self.templates.append(template)

    def compress_data_and_templates(self):
        """ Compress data, exposure and templates

            Before calling must have loaded data, exposure, templates and mask
        """

        # Check user has loaded data, exposure and templates
        assert((len(self.count_map) != 0) | (len(self.exposure_map) != 0)), \
            "Must load a count and exposure map before setting up the scan"
        assert(len(self.templates) != 0), \
            "Must load a template before setting up the scan"

        # Number of pixels is fixed to be the length of the count_map
        self.npix = len(self.count_map)

        # If no mask inserted, set to blank mask
        if len(self.mask_total) == 0:
            print("No mask set; defaulting to a blank mask")
            self.mask_total = np.zeros(self.npix, dtype=bool)

        # Check all inputs have the same length
        assert(len(self.exposure_map) == self.npix), \
            "Exposure map is a different shape to the data"
        assert(len(self.mask_total) == self.npix), \
            "Mask has a different shape to the data"
        for key in self.templates_dict.keys():
            assert(len(self.templates_dict[key]) == self.npix), \
                key + " has a different shape to the data"

        # Compress data - this is used for a Poissonian scan
        temp_data = ma.masked_array(data=self.count_map, mask=self.mask_total)
        # Ensure still an integer array
        self.masked_compressed_data = \
            np.array(temp_data.compressed(), dtype='int32')

        # Check the user has not accidentally masked the entire sky
        assert(len(self.masked_compressed_data != 0)), \
            "The entire sky has been masked - there is nothing to scan over"

        # Divide the map into exposure regions
        self.divide_exposure()

        # Create an exposure region corrected version of the data
        self.masked_compressed_data_expreg = []
        for i in range(self.nexp):
            temp_data_expreg = ma.masked_array(data=self.count_map,
                                               mask=self.expreg_mask[i])
            # Ensure still an integer array
            self.masked_compressed_data_expreg.append(np.array(
                 temp_data_expreg.compressed(), dtype='int32'))

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
        exp_array = np.array([[pix_array[i], self.exposure_map[pix_array[i]]]
                             for i in range(len(pix_array))])
        array_sorted = exp_array[np.argsort(exp_array[:, 1])]

        # Stop code from dividing by more pixels than there are in the ROI
        if self.nexp > len(array_sorted):
            print("nexp cannot be larger than the number of pixels in the ROI")
            print("Setting nexp = ROI size")
            self.nexp = len(array_sorted)

        # Convert from list of exreg pixels to masks (int as used to index)
        array_split = np.array_split(array_sorted, self.nexp)
        expreg_array = [np.array([array_split[i][j][0]
                        for j in range(len(array_split[i]))], dtype='int32')
                        for i in range(len(array_split))]

        temp_expreg_mask = []
        for i in range(self.nexp):
            temp_mask = np.logical_not(np.zeros(self.npix))
            for j in range(len(expreg_array[i])):
                temp_mask[expreg_array[i][j]] = False
            temp_expreg_mask.append(temp_mask)

        self.expreg_mask = temp_expreg_mask

        # Store the total and region by region mean exposure
        expreg_values = [[array_split[i][j][1]
                         for j in range(len(array_split[i]))]
                         for i in range(len(array_split))]
        self.exposure_means_list = [np.mean(expreg_values[i])
                                    for i in range(self.nexp)]
        self.exposure_mean = np.mean(self.exposure_means_list)

    def return_masked_compressed(self, map_to_mask, expreg=False):
        """ Take input map, return masked compressed version and if
            expreg = True a set of such maps broken into exposure regions

            :param map_to_mask: Map to compress
            :param expreg: Whether compressed by exposure region
            :return: A list of compressed maps, sublisted by exposure region
        """

        if not expreg:
            temp_masked_map = ma.masked_array(data=map_to_mask,
                                              mask=self.mask_total)
            return temp_masked_map.compressed()
        else:
            temp_masked_map_list = []
            for i in range(self.nexp):
                temp_masked_map = ma.masked_array(data=map_to_mask,
                                                  mask=self.expreg_mask[i])
                temp_masked_map_list.append(temp_masked_map.compressed())
            return temp_masked_map_list
