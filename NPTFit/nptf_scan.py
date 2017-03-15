###############################################################################
# nptf_scan.py
###############################################################################
#
# Setup and perform the (non-)Poissonian template fit with MultiNest.
# Alternatively configuring the scan defines the log prior cube and the
# likelihood function, which can be extracted for use separately.
#
# The file contains two main sections. The first of these sets up the
# likelihood function (configure_for_scan_internal) and then executes it
# (perform_scan). The second loads the output of the scan, and is used by
# Analysis
#
###############################################################################

from __future__ import print_function
from __future__ import absolute_import

import copy
import numpy as np
import pymultinest
from collections import OrderedDict

from . import pll  # The Poissonian likelihood function
from . import npll  # The non-Poissonian likelihood function
from .config_maps import ConfigMaps  # Setup maps and templates for the run


class NPTFScan(ConfigMaps):
    def __init__(self, tag='Untagged', work_dir=None):
        # Initialise ConfigMaps, creates base directories and allows user
        # to input maps, masks, and templates

        ConfigMaps.__init__(self, tag=tag, work_dir=work_dir)

        self.already_loaded = False

        # Initialise dictionary of models
        self.non_poiss_models = OrderedDict()
        self.poiss_models = OrderedDict()
        self.poiss_models_fixed = OrderedDict()

        # Initialise array of priors and whether they are log or linear flat
        self.priors = []
        self.is_poissonian = []
        self.poiss_list_is_log_prior = []
        self.non_poiss_list_is_log_prior = []

    def add_poiss_model(self, template_name, model_tag, prior_range=[],
                        log_prior=False, fixed=False, fixed_norm=1.0):
        """ Add a Poissonian model corresponding to a template.

            :param template_name: string corresponding to a template added via
                   b.add_template
            :param model_tag: label (LaTeX-ready string) for this model
            :param prior_range: [min_prior_value, max_prior_value]
            :param log_prior: boolean, = True for log spaced priors
            :param fixed: boolean, = True if the template is fixed, not floated
            :param fixed_norm: normalization of template, if fixed
        """

        # Fixed and non-fixed poiss models appended to different dictionaries
        if fixed:
            self.poiss_models_fixed[template_name] = {'fixed_norm': fixed_norm}
        else:
            assert (len(prior_range) != 0), "Fix template or insert a prior"
            self.poiss_models[template_name] = {'prior_range': prior_range,
                                                'log_prior': log_prior,
                                                'model_tag': model_tag}

    def add_non_poiss_model(self, template_name, model_tag, prior_range=[],
                            log_prior=False, dnds_model='specify_breaks',
                            fixed_params=None, units='counts'):
        """ Add a non-Poissonian model corresponding to a template.

            :param template_name: string corresponding to a template added via
                   b.add_template
            :param model_tag: label (LaTeX-ready string) for this model
            :param prior_range: for nb breaks, priors inserted as:
                   [A, n_1, n_2, ... n_{nb+1}, S1, S2, ... S_{nb}],
                   each element a [min, max] list.
                   e.g. [[-6,6],[2.05,30],[-2,1.95],[0.05,30]]
            :param log_prior: boolean, = True for log spaced priors
                   list corresponding to non-Poissonian parameters,
                   e.g. [True,False,False,False]
            :param dnds_model:
                'specify_breaks': Priors are locations of breaks
                'specify_relative_breaks': Priors are relative locations of
                                            breaks
            :param fixed_params: whether to fix non-Poissonian parameters, and if
                   so what value to set them to.
                   e.g. [[0,1.0],[2,4.5]] will set variable 0 to 1.0
                   and variable 2 to 4.5, whilst letting the others
                   float
            :param units: The units in which the prior or fixed value for the
                   breaks is specified. This should be either counts or flux
        """

        if log_prior is False:
            log_prior_list = [False for _ in range(len(model_tag))]
        else:
            log_prior_list = log_prior

        # If fixed parameters, use masks to account for which values to skip
        if fixed_params is None:
            n_params = len(model_tag)
        else:
            n_params = len(model_tag) - len(fixed_params)
            fixed_params = np.array(fixed_params)

        assert (units == 'counts' or units == 'flux'), \
            "units can only be counts or flux"

        # Add template to non_poiss_models dictionary
        self.non_poiss_models[template_name] = {'prior_range': prior_range,
                                                'log_prior': log_prior_list,
                                                'model_tag': model_tag,
                                                'n_params': n_params,
                                                'dnds_model': dnds_model,
                                                'fixed_params': fixed_params,
                                                'n_params_total': n_params,
                                                'units': units}

        # If fixed, account for the fact n_params_total != n_params
        if fixed_params is not None:
            self.non_poiss_models[template_name]['n_params_total'] += len(
                self.non_poiss_models[template_name]['fixed_params'])

    def configure_priors(self):
        """ Set up the cube of priors
        """

        # Configure the poissonian priors
        self.theta_min_poiss = [self.poiss_models[key]['prior_range'][0]
                                for key in self.poiss_model_keys]
        self.theta_max_poiss = [self.poiss_models[key]['prior_range'][1]
                                for key in self.poiss_model_keys]
        # Configure the non-poissonian priors
        self.theta_min_non_poiss = self.flatten(
            [np.array(val['prior_range'])[::, 0]
             for val in self.non_poiss_models.values()])
        self.theta_max_non_poiss = self.flatten(
            [np.array(val['prior_range'])[::, 1]
             for val in self.non_poiss_models.values()])

        # Setup the theta min and max arrays
        self.theta_min = self.theta_min_poiss + self.theta_min_non_poiss
        self.theta_max = self.theta_max_poiss + self.theta_max_non_poiss
        self.theta_interval = list(np.array(self.theta_max) -
                                   np.array(self.theta_min))

        # Define the prior cube
        self.lp = self.prior_cube

    def prior_cube(self, cube, ndim=1, nparams=1):
        """ Cube of priors; motivated by format required by MultiNest, but can
            be used for different MCMC evaluators
        """
        for i in range(ndim):
            cube[i] = cube[i] * self.theta_interval[i] + self.theta_min[i]
        return cube

    def configure_for_scan_internal(self):
        """ Configure data, templates and priors, for evaluation and run over
            log likelihood (with MultiNest or otherwise)
        """

        # Take data and templates, and compress, mask and if specified divide
        # into exposure regions. Function defined in ConfigMaps
        self.compress_data_and_templates()

        # Setup keys and number of parameters for the scan
        self.poiss_model_keys = self.poiss_models.keys()
        self.n_poiss = len(self.poiss_models.keys())
        self.n_non_poiss = np.sum([val['n_params'] for val in
                                   self.non_poiss_models.values()])
        self.n_non_poiss_models = len(self.non_poiss_models.keys())
        self.non_poiss_models_keys = self.non_poiss_models.keys()

        # At this stage we have the mean exposure, so if non-Poissonian
        # templates are defined in terms of flux, adjust this now
        for key in self.non_poiss_models_keys:
            is_flux = self.non_poiss_models[key]['units'] == 'flux'
            is_relative = self.non_poiss_models[key]['dnds_model'] == \
                'specify_relative_breaks'
            # If relative, just adjust the highest break
            if is_flux and is_relative:
                npt_params = self.non_poiss_models[key]['n_params_total']
                npt_breaks = int((npt_params - 2) / 2)
                break_locs = range(npt_breaks + 2, 2 * npt_breaks + 2)
                highest_break = npt_breaks + 2

                # Check if highest break is fixed and if so adjust
                fixed_breaks = 0
                highest_floated = True
                fp_list = self.non_poiss_models[key]['fixed_params']
                if fp_list is not None:
                    for fp in range(len(fp_list)):
                        if fp_list[fp][0] == highest_break:
                            self.non_poiss_models[key]['fixed_params'][fp][1] \
                                *= self.exposure_mean
                            highest_floated = False
                        if fp_list[fp][0] in break_locs:
                            fixed_breaks += 1

                # If floated then adjust
                if highest_floated:
                    floated_breaks = npt_breaks - fixed_breaks
                    floated_params = self.non_poiss_models[key]['n_params']
                    hloc = floated_params - floated_breaks
                    if self.non_poiss_models[key]['log_prior'][hloc] == True:
                        self.non_poiss_models[key]['prior_range'][hloc][0] += \
                            np.log10(self.exposure_mean)
                        self.non_poiss_models[key]['prior_range'][hloc][1] += \
                            np.log10(self.exposure_mean)
                    else:
                        self.non_poiss_models[key]['prior_range'][hloc][0] *= \
                            self.exposure_mean
                        self.non_poiss_models[key]['prior_range'][hloc][1] *= \
                            self.exposure_mean

            # If not relative, adjust all breaks
            if is_flux and not is_relative:
                npt_params = self.non_poiss_models[key]['n_params_total']
                npt_breaks = int((npt_params - 2) / 2)
                break_locs = range(npt_breaks + 2, 2 * npt_breaks + 2)

                # Check if there are any fixed breaks and adjust
                fixed_breaks = 0
                fp_list = self.non_poiss_models[key]['fixed_params']
                if fp_list is not None:
                    for fp in range(len(fp_list)):
                        if fp_list[fp][0] in break_locs:
                            self.non_poiss_models[key]['fixed_params'][fp][1] \
                                *= self.exposure_mean
                            fixed_breaks += 1

                # Adjust floated breaks
                floated_breaks = npt_breaks - fixed_breaks
                floated_params = self.non_poiss_models[key]['n_params']
                for fp in range(floated_params - floated_breaks, floated_params):
                    if self.non_poiss_models[key]['log_prior'][fp] == True:
                        self.non_poiss_models[key]['prior_range'][fp][0] += \
                            np.log10(self.exposure_mean)
                        self.non_poiss_models[key]['prior_range'][fp][1] += \
                            np.log10(self.exposure_mean)
                    else:
                        self.non_poiss_models[key]['prior_range'][fp][0] *= \
                            self.exposure_mean
                        self.non_poiss_models[key]['prior_range'][fp][1] *= \
                            self.exposure_mean

        # Define array of templates describing the spatial distribution of NPT
        # templates - each element is an array of the compressed template in
        # the various exposure regions
        self.NPT_dist_compressed_exp_ary = \
            [self.templates_dict_nested[list(self.non_poiss_models.keys())[i]]
             ['template_masked_compressed_expreg']
             for i in range(self.n_non_poiss_models)]

        # Unfold the full list of non-poissonian fixed parameters
        self.fixed_params_list = []
        n_params_d = 0
        for key in self.non_poiss_models_keys:
            m = self.non_poiss_models[key]
            if m['fixed_params'] is not None:
                fixed_params_d = copy.deepcopy(m['fixed_params'])
                fixed_params_d[::, 0] += n_params_d
                fixed_params_d_with_key = []
                for fix in fixed_params_d:
                    fixed_params_d_with_key += [[fix[0], [key, fix[1]]]]
                self.fixed_params_list += list(fixed_params_d_with_key)
            n_params_d += m['n_params_total']

        self.dnds_model_array = [self.non_poiss_models[key]['dnds_model']
                                 for key in self.non_poiss_models_keys]
        self.model_decompression_key = [[key,
                                         self.poiss_models[key]['log_prior']]
                                        for key in self.poiss_model_keys]

        for key in self.non_poiss_models.keys():
            for j in range(self.non_poiss_models[key]['n_params']):
                self.model_decompression_key += [
                    [key, self.non_poiss_models[key]['log_prior'][j]]]

        self.n_params = len(self.model_decompression_key)

        self.configure_priors()

        print('The number of parameters to be fit is', self.n_params)

    def make_pt_sum_theta(self, theta):
        """ Take the unfixed parameters being scanned by multinest and convert
            these back into objects needed to evaluate the likelihood
            (including adding back in fixed parameters)

            Specifically calculate the sum of the Poissonian templates (fixed
            and floated) and get the list of NPT parameters
        """

        # Setup array of template normalisations (poissonian)
        a_theta = []
        if self.n_poiss != 0:
            a_theta = [self.convert_log(self.model_decompression_key[i],
                                        theta[i]) for i in range(self.n_poiss)]

        # Using this form a compressed map of the sum of all poissonian
        # templates. If there are fixed Poissonian templates add them in too

        pt_sum_compressed_float = np.array(
            [np.sum(list(map(lambda i: a_theta[i][1] *
                             self.templates_dict_nested[a_theta[i][0]]
                             ['template_masked_compressed_expreg'][region],
                             range(len(a_theta)))), axis=0)
             for region in range(self.nexp)])

        pt_sum_compressed_fixed = \
            np.array([np.sum(list(map(lambda key:
                                      self.poiss_models_fixed[key][
                                          'fixed_norm'] *
                                      self.templates_dict_nested[key][
                                          'template_masked_compressed_expreg'][
                                          region],
                                      self.poiss_models_fixed.keys())), axis=0)
                      for region in range(self.nexp)])

        # Combine the two if necessary to build the full PT map
        if len(self.poiss_models_fixed) != 0:
            pt_sum_compressed = pt_sum_compressed_float + \
                                pt_sum_compressed_fixed
        else:
            pt_sum_compressed = pt_sum_compressed_float

        # If no Poissonian templates, pt_sum_compressed set to a zero array
        if (self.n_poiss + len(self.poiss_models_fixed)) == 0:
            pt_sum_compressed = np.array([np.zeros(
                len(self.masked_compressed_data_expreg[region]))
                                          for region in range(self.nexp)])

        # Setup array of NPT parameters - adding the fixed parameters back in
        # Note that through the use of convert_log, theta_ps here is an array
        # of [tag, value]
        theta_ps = []
        if self.n_non_poiss != 0:
            theta_ps = [self.convert_log(self.model_decompression_key[i],
                                         theta[i]) for i in
                        range(self.n_poiss, self.n_params)]
            # Add the fixed parameters back into theta_ps in the correct spot
            if len(self.fixed_params_list) != 0:
                index = 0
                for key in self.non_poiss_models_keys:
                    if self.non_poiss_models[key]['fixed_params'] is not None:
                        for fixed in self.non_poiss_models[key]['fixed_params']:
                            theta_ps.insert(int(fixed[0]) + index, [key,
                                                                    fixed[1]])
                    index += self.non_poiss_models[key]['n_params_total']

        return pt_sum_compressed, theta_ps

    def model_parameters_nbreak(self, theta_ps):
        """ Take the flat theta_ps array determined in make_pt_sum_theta and
            convert it into an array of separate NPT parameter arrays
        """

        # To interpret the organisation below, recall that an NPT template with
        # nb breaks has 2nb+2 parameters, organised as follows:
        # [A, n_{1}, ..., n_{nb+1}, S_{1}, ..., S_{nb}]
        # The relative positions of these in the array is:
        # [0, 1, ..., nb+1, nb+2, ..., 2nb+1]

        # Determine the number of breaks and total parameters for each NPT
        nbreak_ary = [int((val['n_params_total'] - 2) / 2.)
                      for val in self.non_poiss_models.values()]
        nparams_ary = [int(val['n_params_total'])
                       for val in self.non_poiss_models.values()]

        # Determine where each NPT model begins in the full theta_ps array
        offset_indices_ary = [0]
        for i in range(1, self.n_non_poiss_models):
            offset_indices_ary.append(offset_indices_ary[i - 1] +
                                      nparams_ary[i - 1])

        theta_ps_ary = []

        # Add in parameters for each NPT separately
        # Account for whether the user specified absolute or relative breaks
        for i in range(self.n_non_poiss_models):
            if self.dnds_model_array[i] == 'specify_breaks':
                # If absolute, just add the 2nb+2 parameters without alteration
                theta_ps_ary.append(theta_ps[offset_indices_ary[i]:
                                    offset_indices_ary[i] + nparams_ary[i]])

            elif self.dnds_model_array[i] == 'specify_relative_breaks':
                # If instead relative, weight all breaks by the last break
                # This is done by adding together three arrays:
                # 1. Unchanged parameters: [A, n_{1}, ..., n_{nb+1}], which
                # has length nb+2 (=nparams-nbreak)
                # 2. The adjusted breaks: [S_{1}, ..., S_{nb-1}], which has
                # length nb-1
                # 3. The unadjusted break: [S_{nb}], which has length 1
                theta_ps_ary.append(theta_ps[offset_indices_ary[i]:
                    offset_indices_ary[i] + nparams_ary[i]-nbreak_ary[i]] +
                    [theta_ps[offset_indices_ary[i] + nparams_ary[i] - 1] *
                    np.product(theta_ps[offset_indices_ary[i] + j:
                    offset_indices_ary[i] + nparams_ary[i] - 1]) for j in range(
                    offset_indices_ary[i] + nparams_ary[i] - nbreak_ary[i],
                    offset_indices_ary[i] + nparams_ary[i] - 1)] +
                    [theta_ps[offset_indices_ary[i]+nparams_ary[i] - 1]])
        return theta_ps_ary

    def make_ll(self):
        """ Pass all details to the likelihood evaluator, and define the total
            ll to be the sum of that in each of the exposure regions
        """

        ll = 0.0
        for i in range(self.nexp):
            # For each NPT template adjust the breaks to account for the
            # difference in exposure
            theta_ps_expreg = [[self.theta_ps[j][0]] +
                               list(self.theta_ps[j][1:self.nbreak_ary[j] + 2]) +
                               list(np.array(
                                    self.theta_ps[j][self.nbreak_ary[j] + 2:]) *
                                    (self.exposure_means_list[
                                        i] / self.exposure_mean))
                               for j in
                               range(len(self.NPT_dist_compressed_exp_ary))]

            # In evaluating the likelihood extract the exposure region i
            # version of each parameter
            ll += npll.log_like(self.PT_sum_compressed[i], theta_ps_expreg,
                                self.f_ary, self.df_rho_div_f_ary,
                                [NPT[i] for NPT in
                                 self.NPT_dist_compressed_exp_ary],
                                self.masked_compressed_data_expreg[i])
        return ll

    def log_like_nptf(self, theta, ndim=1, nparams=1):
        """ NPTF likelihood function
        """

        # Determine PT and NPT contribution and pass to the likelihood function
        pt_sum_compressed, theta_ps_marked = self.make_pt_sum_theta(theta)
        # theta_ps_marked is an array of [tag, value], extract values
        theta_ps = list(np.vectorize(float)(np.array(theta_ps_marked)[::, 1]))
        self.theta_ps = self.model_parameters_nbreak(theta_ps)
        self.nbreak_ary = [int((len(self.theta_ps[j]) - 2) / 2.)
                           for j in range(len(self.theta_ps))]
        self.PT_sum_compressed = pt_sum_compressed

        return self.make_ll()

    def log_like_ptf(self, theta, ndim=1, nparams=1):
        """ PTF likelihood function
        """

        # Determine PT and pass to the likelihood function
        pt_sum_compressed_arr, _ = self.make_pt_sum_theta(theta)
        # Take out the 0th element as the Poissonian likelihood does not loop
        # through exposure regions
        pt_sum_compressed = pt_sum_compressed_arr[0]

        return pll.log_like_poissonian(pt_sum_compressed,
                                       self.masked_compressed_data)

    def perform_scan(self, run_tag=None, nlive=100, pymultinest_options=None):
        """ When called creates the directories for and then performs the scan

            run_tag: label of file in which the chains output is stored
            nlive: number of live points to user during the scan. The default
            value of 100 is chosen to speed up testing, but for actual runs
            a larger value is recommended
            pymultinest_options: Custom user inputs passed to multinest; must
            be inserted in the form of a dictionary, as for the default below
        """

        self.run_tag = run_tag
        self.make_dirs_for_run(run_tag)

        if not pymultinest_options:
            # Set defaults
            pymultinest_options = {'importance_nested_sampling': False,
                                   'resume': False, 'verbose': True,
                                   'sampling_efficiency': 'model',
                                   'init_MPI': False, 'evidence_tolerance': 0.5,
                                   'const_efficiency_mode': False}
        else:
            pass  # Use passed parameters

        # Run MultiNest
        pymultinest.run(self.ll, self.lp, self.n_params,
                        outputfiles_basename=self.chains_dir_for_run,
                        n_live_points=nlive, **pymultinest_options)

    ############################################################
    # Functions below are associated with preparing the analysis
    ############################################################

    def load_scan(self, run_tag=None):
        """ Load the details of a scan, if not already loaded
            This sets up the output of the scan in a format useful for Analysis
        """

        if self.already_loaded:
            pass
        else:
            self.run_tag = run_tag
            self.make_dirs_for_run(run_tag)

            # Load the statistics and samples
            self.a = \
                pymultinest.Analyzer(n_params=self.n_params,
                                     outputfiles_basename=self.chains_dir_for_run)
            self.s = self.a.get_stats()
            self.chain_file = self.chains_dir_for_run + '/post_equal_weights.dat'
            self.samples = np.array(np.loadtxt(self.chain_file)[:, :-1])

            # Determine the median value of each parameter
            self.medians = [self.s['marginals'][i]['median']
                            for i in range(self.n_params)]

            # Determine the template normalisations
            self.calculate_norms()

            # Account for fixed parameters if necessary
            if len(self.fixed_params_list) > 0:
                self.fix_samples()
                self.fix_medians()
                self.n_params += len(self.fixed_params_list)
                self.fix_model_decompression()

            self.already_loaded = True

    def calculate_norms(self):
        """ Determine the normalisation associated with the models
        """

        # Determine the list of P and NP parameters and combine the list
        self.poiss_params = np.array(
            [self.poiss_models[key]['model_tag']
             for key in self.poiss_model_keys])

        self.non_poiss_params = self.flatten(np.array(
            [val['model_tag'] for val in self.non_poiss_models.values()]))

        self.params = np.array(list(self.poiss_params) +
                               list(self.non_poiss_params))

        # Determine which parameters had log priors and combine the list
        self.poiss_list_is_log_prior = np.array(
            [self.poiss_models[key]['log_prior']
             for key in self.poiss_model_keys])

        self.non_poiss_list_is_log_prior = np.array(
            self.flatten([val['log_prior']
                          for val in self.non_poiss_models.values()]))

        self.medians_not_log = \
            self.convert_log_list(self.medians,
                                  np.concatenate([self.poiss_list_is_log_prior,
                                                  self.non_poiss_list_is_log_prior]))

        # In terms of these extract the P and NP norms
        self.norms_poiss = {self.model_decompression_key[i][0]:
                            self.medians_not_log[i] for i in range(self.n_poiss)}

        self.norms_non_poiss = OrderedDict()
        index = 0  # Multiple parameters for each NPT, index counts these
        for i in range(len(self.non_poiss_models.keys())):
            self.norms_non_poiss[list(self.non_poiss_models.keys())[i]] = \
                [self.medians_not_log[self.n_poiss + index + j]
                 for j in range(
                    self.non_poiss_models[list(self.non_poiss_models.keys())[i]][
                        'n_params'])]
            index += \
                self.non_poiss_models[list(self.non_poiss_models.keys())[i]][
                    'n_params']

    def fix_samples(self):
        """ Account for fixed parameters in the samples
        """

        new_samples = []
        for i in range(len(self.samples)):
            sam = list(self.samples[i])
            index = self.n_poiss
            for key in self.non_poiss_models_keys:
                if self.non_poiss_models[key]['fixed_params'] is not None:
                    for fixed in self.non_poiss_models[key]['fixed_params']:
                        sam.insert(int(fixed[0]) + index, fixed[1])
                index += self.non_poiss_models[key]['n_params_total']
            new_samples += [np.array(sam)]
        self.samples = np.array(new_samples)

    def fix_medians(self):
        """ Account for fixed parameters in medians
        """

        index = self.n_poiss
        for key in self.non_poiss_models_keys:
            if self.non_poiss_models[key]['fixed_params'] is not None:
                for fixed in self.non_poiss_models[key]['fixed_params']:
                    self.medians.insert(int(fixed[0]) + index, fixed[1])
            index += self.non_poiss_models[key]['n_params_total']

    def fix_model_decompression(self):
        """ Account for fixed parameters in model_decompression_key
        """

        index = int(self.n_poiss)
        for key in self.non_poiss_models_keys:
            for fixed in self.fixed_params_list:
                if fixed[1][0] == key:
                    self.model_decompression_key.insert(index +
                                                        int(fixed[0]),
                                                        [key, False])
            index += self.non_poiss_models[key]['n_params_total']

    @staticmethod
    def convert_log_list(the_list, is_log):
        """ Take a list of parameters and accounts for those with log priors
        """

        new_list = []
        for i in range(len(the_list)):
            if is_log[i]:
                new_list.append(10 ** the_list[i])
            else:
                new_list.append(the_list[i])
        return new_list

    @staticmethod
    def convert_log(key, val):
        """ Convert values to physical values if a log prior is used
        """

        tag, log_prior = key
        if log_prior:
            return [tag, 10 ** val]
        else:
            return [tag, val]

    @staticmethod
    def flatten(the_list):
        """ Flatten out a nested array
        """

        the_list = list(the_list)
        return [item for sublist in the_list for item in list(sublist)]
