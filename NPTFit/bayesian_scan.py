###############################################################################
# bayesian_scan.py
###############################################################################
#
# Setup and perform the (non-)Poissonian template fit scan with MultiNest.
#
###############################################################################

import copy

import numpy as np
from mpi4py import MPI
import pymultinest

import PLL as log_like_poiss
from config_maps import config_maps


class bayesian_scan(config_maps):
    """ The file that sets up and executes the scan """

    def __init__(self, nlive=100, *args, **kwargs):
        config_maps.__init__(self, *args, **kwargs)

        self.priors = []
        self.is_log_prior = []
        self.is_poissonian = []

        self.poiss_list_template_number = []
        self.poiss_list_prior_range = []
        self.poiss_list_model_tag = []
        self.poiss_list_is_log_prior = []

        self.non_poiss_list_template_number = []
        self.non_poiss_list_prior_range = []
        self.non_poiss_list_model_tag = []
        self.non_poiss_list_is_log_prior = []

        self.nlive = nlive
        self.resume = False

        self.non_poiss_models = {}
        self.poiss_models = {}
        self.poiss_models_fixed = {}

        self.already_loaded = False

    def add_poiss_model(self, template_name, model_tag, prior_range, log_prior=False, fixed=False, fixed_norm=1.0):
        """ Add a Poissonian model corresponding to a template

            template_name: Name of poissonian template. Should be a string corresponding a template added using b.add_template
            model_tag: Label, for triangle plots (LaTeX-ready string) TODO: should we completely remove this?
            prior_range: [min_prior_value, max_prior_value]
            log_prior: Whether prior values are in log-space
            fixed: Whether template added is fixed or floating
            fixed_norm: Normalization of template, if fixed
        """
        if fixed:
            self.poiss_models_fixed[template_name] = {'fixed_norm': fixed_norm}
        else:
            self.poiss_models[template_name] = {'prior_range': prior_range, 'log_prior': log_prior, 'model_tag': model_tag, 'poiss': True}

    def add_non_poiss_model(self, template_name, model_tag, prior_range,
                            log_prior=False,
                            dnds_model='specify_relative_breaks',
                            fixed_params=None):
        """ Add a non-Poissonian model corresponding to a template

            template_name: Name of poissonian template. Should be a string corresponding a template added using b.add_template()
            model_tag: Label, for triangle plots (LaTeX-ready string) TODO: should we completely remove this?
            prior_range: Priors as [A, n_1, n_2, ... n_{nb+1}, S1, S2, ... S_{nb}], each element a [min, max] list. e.g. [[-6,6],[2.05,30],[-2,1.95],[0.05,30]]
            log_prior: Whether prior values are in log-space, list corresponding to non-Poissonian parameters e.g. [True,False,False,False]
            dnds_model:
                'specify_breaks': Priors are locations of breaks
                'specifiy_relative_breaks': Priors are relative locations of breaks
            fixed_params: Whether to fix non-Poissonian parameters.
                Liist like [[0,1.],[2,4.5]].  This will fixed var 0 to 1.0, var 1 to 4.5 while letting the other variables float
        """
        if log_prior is False:
            log_prior_list = [False for i in range(len(model_tag))]
        else:
            log_prior_list = log_prior
        if fixed_params is None:
            n_params = len(model_tag)
        else:
            n_params = len(model_tag) - len(fixed_params)

            total_params = np.arange(len(prior_range))
            fixed_params_wheres = np.array(fixed_params)[::, 0]
            mask = np.in1d(total_params, fixed_params_wheres)
            not_fixed_params = np.where(~mask)[0]
            fixed_params = np.array(fixed_params)

        self.non_poiss_models[template_name] = {'prior_range': prior_range,
                                                'log_prior': log_prior_list,
                                                'model_tag': model_tag,
                                                'n_params': n_params,
                                                'poiss': False,
                                                'dnds_model': dnds_model,
                                                'fixed_params': fixed_params}

        self.non_poiss_models[template_name]['n_params_total'] = \
            self.non_poiss_models[template_name]['n_params']
        if fixed_params is not None:
            self.non_poiss_models[template_name]['n_params_total'] += len(
                self.non_poiss_models[template_name]['fixed_params'])

    def configure_priors(self):
        """Sets up priors as sequired by MultiNest
        """

        self.configure_poissonian_priors()
        self.configure_non_poissonian_priors()
        self.theta_min = self.theta_min_poiss + self.theta_min_non_poiss
        self.theta_max = self.theta_max_poiss + self.theta_max_non_poiss
        self.theta_interval = list(np.array(self.theta_max) - np.array(self.theta_min))

        self.lp = self.log_prior

    def configure_poissonian_priors(self):
        """Sets up Poissonian priors
        """
        self.theta_min_poiss = [self.poiss_models[key]['prior_range'][0] for key in self.poiss_model_keys]
        self.theta_max_poiss = [self.poiss_models[key]['prior_range'][1] for key in self.poiss_model_keys]

    def configure_non_poissonian_priors(self):
        """Sets up non-Poissonian priors
        """
        self.theta_min_non_poiss = self.flatten(
            [np.array(val['prior_range'])[::, 0] for val in self.non_poiss_models.values()])
        self.theta_max_non_poiss = self.flatten(
            [np.array(val['prior_range'])[::, 1] for val in self.non_poiss_models.values()])

    def log_prior(self, cube, ndim, nparams):
        """Priors as required by MultiNest
        """
        for i in range(ndim):
            cube[i] = cube[i] * self.theta_interval[i] + self.theta_min[i]

    def configure_for_edep_scan(self):
        """Configure data/templates and priors for MultiNest run
        """
        self.compress_data_and_templates()
        self.poiss_model_keys = self.poiss_models.keys()
        self.n_poiss = len(self.poiss_models.keys())
        self.non_poiss_list_num_model_params = [d['n_params'] for d in self.non_poiss_models.values()]
        self.n_non_poiss = np.sum(self.non_poiss_list_num_model_params)
        self.n_non_poiss_models = len(self.non_poiss_models.keys())
        self.non_poiss_models_keys = self.non_poiss_models.keys()

        self.fixed_params_list = []
        n_params_d = 0
        for key in self.non_poiss_models_keys:
            d = self.non_poiss_models[key]
            if d['fixed_params'] is not None:
                fixed_params_d = copy.deepcopy(d['fixed_params'])
                fixed_params_d[::, 0] += n_params_d
                fixed_params_d_with_key = []
                for fix in fixed_params_d:
                    fixed_params_d_with_key += [[fix[0], [key, fix[1]]]]
                self.fixed_params_list += list(fixed_params_d_with_key)
            n_params_d += d['n_params_total']  # 'n_params'

        self.dnds_model_array = [self.non_poiss_models[key]['dnds_model'] for key in self.non_poiss_models_keys]
        self.model_decompression_key = [[key, self.poiss_models[key]['log_prior']] for key in
                                        self.poiss_model_keys]
        for key in self.non_poiss_models.keys():
            for j in range(self.non_poiss_models[key]['n_params']):
                self.model_decompression_key += [
                    [key, self.non_poiss_models[key]['log_prior'][j]]]

        self.n_params = len(self.model_decompression_key)
        self.configure_priors()

        print 'The number of parameters is ', self.n_params

    def convert_log(self, key, val):
        tag, log_prior = key
        if log_prior:
            return [tag, 10 ** val]
        else:
            return [tag, val]

    def make_xbg_edep_exposure(self, theta):
        """Return sum of Poissonian templates and theta_PS
        """
        theta_list = theta
        A_theta = []
        theta_PS = []

        if self.n_poiss > 0:
            A_theta = [self.convert_log(self.model_decompression_key[i], theta_list[i]) for i in range(self.n_poiss)]

        if self.n_non_poiss > 0:
            theta_PS = [self.convert_log(self.model_decompression_key[i], theta_list[i]) for i in range(self.n_poiss, self.n_params)]
            if len(self.fixed_params_list) > 0:
                index = 0
                for key in self.non_poiss_models_keys:
                    if self.non_poiss_models[key]['fixed_params'] is not None:
                        for fixed in self.non_poiss_models[key]['fixed_params']:
                            theta_PS.insert(int(fixed[0]) + index, [key, fixed[1]])
                    index += self.non_poiss_models[key]['n_params_total']  # n_params

        xbg_PSF_compressed_model = np.array([[np.sum(map(lambda i: A_theta[i][1] * self.templates_dict_nested[A_theta[i][0]]['template_masked_compressed_expreg'][region], range(len(A_theta))), axis=0) for region in range(self.nexp)]])

        xbg_PSF_compressed_fixed = np.array([[np.sum(map(lambda key: self.poiss_models_fixed[key]['fixed_norm'] * self.templates_dict_nested[key]['template_masked_compressed_expreg'][region], self.poiss_models_fixed.keys()), axis=0) for region in range(self.nexp)]])

        if len(self.poiss_models_fixed) != 0:
            xbg_PSF_compressed = xbg_PSF_compressed_model + xbg_PSF_compressed_fixed
        else:
            xbg_PSF_compressed = xbg_PSF_compressed_model

        return xbg_PSF_compressed, theta_PS

    def single_model_parameters_edep_nbreak(self, params, energy_bin=0):
        """Return parameters as required by MultiNest, for single non-Poissonian model
        """
        nbreak = self.nbreak
        nfix_index = (nbreak - 1) * 2 + 3
        return [params[0] / float(params[nfix_index + energy_bin])] + params[1:nfix_index - nbreak + 1] + [params[nfix_index + energy_bin] * np.product(params[i:nfix_index]) for i in range(nfix_index - nbreak + 1, nfix_index)] + [params[nfix_index + energy_bin]]

    def model_parameters_edep_nbreak(self, theta_PS, energy_bin=0):
        """Returns parameters as required by MultiNest, for arbitrary number of non-Poissonian models
        """
        nbreak_ary = [int((val['n_params_total'] - 2)/2.) for val in self.non_poiss_models.values()]
        nfix_index_ary = [(nbreak_ary[i] - 1) * 2 + 3 for i in range(len(nbreak_ary))]
        offsets = [0]
        for i in range(1, self.n_non_poiss_models):
            offsets.append(offsets[i-1] + val['n_params_total'] for val in self.non_poiss_models.values())

        if self.dnds_model_array[0] == 'specify_breaks':
            return [[theta_PS[offset] / float(theta_PS[offset + 2*nbreak + 1])] + theta_PS[1 + offset: offset + 2*nbreak + 2] for offset, nbreak in zip(offsets, nbreak_ary)]  # TODO: remove 0 when fixed multiple NPT!
        else:
            return [self.single_model_parameters_edep_nbreak(theta_PS[offset:nfix_index+offset] + theta_PS[nfix_index + energy_bin + offset: nfix_index + nbreak + energy_bin + offset]) for offset, nbreak, nfix_index in zip(offsets, nbreak_ary, nfix_index_ary)]

    def log_like_NPTF_edep(self, theta, ndim, nparams):
        """Reurn log likelihood as required by MultiNest
        """
        xbg_PSF_compressed, theta_PS_marked = self.make_xbg_edep_exposure(theta)
        theta_PS = list(np.vectorize(float)(np.array(theta_PS_marked)[::, 1]))

        self.NPTF_inst.theta_PS = self.model_parameters_edep_nbreak(theta_PS)
        self.NPTF_inst.xbg_PSF_compressed = xbg_PSF_compressed[0]  # Remove [0] when xbg_PSF_compressed is no longer an array

        ll = self.NPTF_inst.make_ll_PDF_PSF()

        return ll

    def log_like_poisson_edep(self, theta, ndim, nparams):
        xbg_PSF_compressed = self.make_xbg_edep(theta)[0]

        ll = log_like_poiss.log_likelihood_poissonian(xbg_PSF_compressed[0], np.array(self.masked_compressed_data, dtype='int32'))  # 0 index is energy -- remove when axed
        return ll

    def perform_scan(self, run_tag='test', edep=False, importance_nested_sampling=False, evidence_tolerance=0.5,
                     const_efficiency_mode=False, log_zero=None):
        self.run_tag = run_tag
        self.make_dirs_for_run(run_tag)

        if log_zero is None:
            pymultinest.run(self.ll, self.lp, self.n_params, importance_nested_sampling=importance_nested_sampling,
                            resume=self.resume, verbose=True, sampling_efficiency='model', n_live_points=self.nlive,
                            outputfiles_basename=self.chains_dir_for_run, init_MPI=False,
                            evidence_tolerance=evidence_tolerance, const_efficiency_mode=const_efficiency_mode)
        else:
            pymultinest.run(self.ll, self.lp, self.n_params, importance_nested_sampling=importance_nested_sampling,
                            resume=self.resume, verbose=True, sampling_efficiency='model', n_live_points=self.nlive,
                            outputfiles_basename=self.chains_dir_for_run, init_MPI=False,
                            evidence_tolerance=evidence_tolerance, const_efficiency_mode=const_efficiency_mode, log_zero=log_zero)

    def load_scan(self, run_tag='test', edep=False, *args, **kwargs):
        if self.already_loaded:
            pass
        else:
            self.run_tag = run_tag
            self.make_dirs_for_run(run_tag)

            self.a = pymultinest.Analyzer(n_params=self.n_params, outputfiles_basename=self.chains_dir_for_run)
            self.s = self.a.get_stats()
            self.chain_file = self.chains_dir_for_run + '/post_equal_weights.dat'
            self.samples = np.array(np.loadtxt(self.chain_file)[:, :-1])

            self.calculate_medians()

            print 'Len of fixed params list: ', len(self.fixed_params_list)
            self.calculate_norms(edep)

            if len(self.fixed_params_list) > 0:
                print 'fixing!'
                self.fix_fixed()
            self.already_loaded = True

    def fix_fixed(self):
        self.fix_samples()
        self.fix_medians()
        self.n_params += len(self.fixed_params_list)
        self.fix_model_decompression()

    def fix_model_decompression(self):
        index = int(self.n_poiss)
        for key in self.non_poiss_models_keys:
            for fixed in self.fixed_params_list:
                if fixed[1][0] == key:
                    self.model_decompression_key.insert(index + int(fixed[0]), [key, False])
            index += self.non_poiss_models[key]['n_params_total']  # n_params

    def fix_samples(self):
        new_samples = []
        for i in range(len(self.samples)):
            sam = list(self.samples[i])
            index = self.n_poiss
            for key in self.non_poiss_models_keys:
                if self.non_poiss_models[key]['fixed_params'] is not None:
                    for fixed in self.non_poiss_models[key]['fixed_params']:
                        sam.insert(int(fixed[0]) + index, fixed[1])
                index += self.non_poiss_models[key]['n_params_total']  # n_params
            new_samples += [np.array(sam)]
        self.samples = np.array(new_samples)

    def fix_medians(self):
        index = self.n_poiss
        for key in self.non_poiss_models_keys:
            if self.non_poiss_models[key]['fixed_params'] is not None:
                for fixed in self.non_poiss_models[key]['fixed_params']:
                    self.medians.insert(int(fixed[0]) + index, fixed[1])
            index += self.non_poiss_models[key]['n_params_total']  # n_params

    def calculate_medians(self):
        self.medians = [self.s['marginals'][i]['median'] for i in range(self.n_params)]
        print 'the medians are ', self.medians
        print 'the medians (from chains) are ', np.median(self.samples, axis=0)
        print 'the number of parameters is ', self.n_params

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

    def calculate_norms(self, edep):
        if edep:
            self.n_poiss_base = len(self.poiss_models.keys())
            self.poiss_comp_numbers = {self.poiss_model_keys[i]: i for i in range(self.n_poiss)}.copy()
            self.test_array = self.poiss_comp_numbers.copy()
            self.poiss_list_is_log_prior = np.array(
                [self.poiss_models[key]['log_prior'] for key in self.poiss_model_keys])
            self.poiss_params = np.array([self.poiss_models[key]['model_tag'] for key in self.poiss_model_keys])
        else:
            self.poiss_comp_numbers = {self.poiss_models.keys()[i]: i for i in range(self.n_poiss)}
            self.poiss_list_is_log_prior = np.array([val['log_prior'] for val in self.poiss_models.values()])
            self.poiss_params = np.array([val['model_tag'] for val in self.poiss_models.values()])
        self.non_poiss_params = self.flatten(np.array([val['model_tag'] for val in self.non_poiss_models.values()]))
        self.params = np.array(list(self.poiss_params) + list(self.non_poiss_params))
        self.non_poiss_list_is_log_prior = np.array(
            self.flatten([val['log_prior'] for val in self.non_poiss_models.values()]))
        self.medians_not_log = self.convert_log_list(self.medians, np.concatenate(
            [self.poiss_list_is_log_prior, self.non_poiss_list_is_log_prior]))
        self.norms_poiss = {self.model_decompression_key[i][0]: self.medians_not_log[i] for i in range(self.n_poiss)}

        self.norms_non_poiss = {}
        k = 0
        for i in range(len(self.non_poiss_models.keys())):
            self.norms_non_poiss[self.non_poiss_models.keys()[i]] = [self.medians_not_log[self.n_poiss + k + j] for j in
                                                                     range(self.non_poiss_models[
                                                                               self.non_poiss_models.keys()[i]][
                                                                               'n_params'])]
            k += self.non_poiss_models[self.non_poiss_models.keys()[i]]['n_params']
        print 'self.poiss_comp_numbers: ', self.poiss_comp_numbers
