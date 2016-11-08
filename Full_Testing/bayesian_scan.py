import sys, os

import healpy as hp
import numpy as np
import set_dirs as sd

import pymultinest
from mpi4py import MPI
import triangle
import matplotlib.pyplot as plt
import copy

import inspect
import types

import NPLL as log_like

import logging

from config_file import config

import time as time

import datetime


class bayesian_scan(config):
    """ The file that sets up and executes the scan """

    def __init__(self, nlive=100, k_max='False', Sc='False', *args, **kwargs):
        config.__init__(self, *args, **kwargs)

        self.k_max = k_max
        self.Sc = Sc

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



        # new stuff
        self.non_poiss_models = {}
        self.poiss_models = {}

        ###
        self.already_loaded = False


        ########
        ###temp
        self.index = 0

    def set_params(self, params):
        self.params = params
        self.n_params = len(params)

    def set_is_log_priors(self, is_logs):
        # prior should be of the form [[min, max]] or [[min, max],[min, max]], etc.
        for is_log in is_logs:
            self.is_log_prior.append(is_log)

    def adding_energy_bins(self, array1):
        new_array = []
        energy_array = [i * self.energy_bin_spacing for i in range(self.number_energy_bins)]
        for i in range(len(self.templates)):
            new_array.append(
                [np.sum(array1[i][energy_array[j]:energy_array[j] + self.energy_bin_spacing], axis=0) for j in
                 range(len(energy_array))])
        return new_array

    def make_poiss_models_edep(self):
        self.poiss_list_template_number_edep = []
        self.poiss_list_prior_range_edep = []
        self.poiss_list_model_tag_edep = []
        self.poiss_list_is_log_prior_edep = []
        self.templates_edep = []
        self.templates_edep_masked_compressed = []


        # Need to add compression and masking to templates_edep
        for i in range(len(self.templates)):
            for j in range(self.number_energy_bins):
                self.templates_edep.append(self.templates[i][j])

        self.templates_masked_compressed_organized = self.adding_energy_bins(self.templates_masked_compressed)

        for i in range(len(self.templates)):
            for j in range(self.number_energy_bins):
                self.templates_edep_masked_compressed.append(self.templates_masked_compressed_organized[i][j])

        for i in range(len(self.poiss_list_template_number)):
            for j in range(self.number_energy_bins):
                self.poiss_list_template_number_edep.append(self.poiss_list_template_number[i] + j)
                self.poiss_list_prior_range_edep.append(self.poiss_list_prior_range[i])
                self.poiss_list_model_tag_edep.append(self.poiss_list_model_tag[i] + '$^{(' + str(j) + ')}$')
                self.poiss_list_is_log_prior_edep.append(self.poiss_list_is_log_prior[i])

    def add_poiss_model(self, template_name, model_tag, prior_range, log_prior=False):
        self.poiss_models[template_name] = {'prior_range': prior_range, 'log_prior': log_prior, 'model_tag': model_tag,
                                            'poiss': True}

    def add_non_poiss_model(self, template_name, model_tag, prior_range, log_prior=False, dnds_model='specify_relative_breaks',
                            fixed_params=None):
        #dnds_model:
        #1. specify_breaks: Priors are locations of breaks (only works for 1 energy bin for now)
        #2. specifiy_relative_breaks: priors are relative locations of breaks (can handle multple energy bins)
        # Fixed params should be a list like [ [0,1],[2,4],[3,3.2]].  This will fixed var 0 to 1, var 1 to 4, and var 2 to 5, while letting the other variables float
        if log_prior == False:
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

            # prior_range = list( np.array(prior_range)[not_fixed_params] )
            # log_prior_list = list( np.array(log_prior_list)[not_fixed_params] )

        self.non_poiss_models[template_name] = {'prior_range': prior_range, 'log_prior': log_prior_list,
                                                'model_tag': model_tag, 'n_params': n_params, 'poiss': False,
                                                'dnds_model': dnds_model, 'fixed_params': fixed_params}

        self.non_poiss_models[template_name]['n_params_total'] = self.non_poiss_models[template_name]['n_params']
        if fixed_params is not None:
            self.non_poiss_models[template_name]['n_params_total'] += len(
                self.non_poiss_models[template_name]['fixed_params'])

    def set_is_poissonian(self, is_poiss):
        # prior should be of the form [[min, max]] or [[min, max],[min, max]], etc.
        for is_pois in is_poiss:
            self.is_poissonian.append(is_pois)

    def set_priors(self, prior):
        for pri in prior:
            self.priors.append(pri)

    def configure_priors(self, number_energy_bins=1):

        self.configure_poissonian_priors()
        self.configure_non_poissonian_priors()
        self.theta_min = self.theta_min_poiss + self.theta_min_non_poiss  # *number_energy_bins
        self.theta_max = self.theta_max_poiss + self.theta_max_non_poiss  # *number_energy_bins
        self.theta_interval = list(np.array(self.theta_max) - np.array(self.theta_min))

        self.lp = self.log_prior

    def configure_poissonian_priors(self):
        self.theta_min_poiss = [self.poiss_models[key]['prior_range'][0] for key in self.poiss_model_keys]
        self.theta_max_poiss = [self.poiss_models[key]['prior_range'][1] for key in self.poiss_model_keys]

    def configure_non_poissonian_priors(self):
        self.theta_min_non_poiss = self.flatten(
            [np.array(val['prior_range'])[::, 0] for val in self.non_poiss_models.values()])
        self.theta_max_non_poiss = self.flatten(
            [np.array(val['prior_range'])[::, 1] for val in self.non_poiss_models.values()])

    def log_prior(self, cube, ndim, nparams):
        for i in range(ndim):
            cube[i] = cube[i] * self.theta_interval[i] + self.theta_min[i]

    def log_prior_polychord(self, cube):
        for i in range(len(cube)):
            cube[i] = cube[i] * self.theta_interval[i] + self.theta_min[i]
        return cube

    def set_ll(self, mode='Poisson'):
        if mode == 'Poisson':
            self.ll = self.log_like_poisson

    def configure_for_edep_scan(self):

        self.compress_data()
        self.compress_data_exposure()
        self.compress_templates()
        self.number_energy_bins = 1

        for en in range(self.number_energy_bins):
            for key in self.templates_dict_nested.keys():
                self.templates_dict_nested[key + '-' + str(en)] = self.templates_dict_nested[key]
                self.templates_dict[key + '-' + str(en)] = self.templates_dict[key]

        self.edep_poiss_models = {}
        for en in range(self.number_energy_bins):
            for key in self.poiss_models.keys():
                self.edep_poiss_models[key + '-' + str(en)] = self.poiss_models[key].copy()
                print self.edep_poiss_models[key + '-' + str(en)]['model_tag']
                self.edep_poiss_models[key + '-' + str(en)]['model_tag'] += '$^{' + str(en) + '}$'
        self.poiss_models_bare_keys = self.poiss_models.keys()
        self.poiss_model_keys = []
        for en in range(self.number_energy_bins):
            for key in self.poiss_models_bare_keys:
                self.poiss_model_keys.append(key + '-' + str(en))
        self.poiss_models = self.edep_poiss_models.copy()


        # self.n_poiss = len(self.poiss_models.keys())*self.number_energy_bins
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
                                        self.poiss_model_keys]  # *self.number_energy_bins
        for key in self.non_poiss_models.keys():
            for j in range(self.non_poiss_models[key]['n_params']):
                print j, self.non_poiss_models[key]['log_prior']
                self.model_decompression_key += [
                    [key, self.non_poiss_models[key]['log_prior'][j]]]  # False means non-poiss template

        self.n_params = len(self.model_decompression_key)

        self.set_k_max_edep()
        self.configure_priors(number_energy_bins=self.number_energy_bins)

        print 'The number of parameters is ', self.n_params

    def configure_for_scan(self):

        self.compress_data()
        # self.sum_data_over_energy_range()
        self.compress_templates()
        # self.sum_templates_over_energy_range()


        self.n_poiss = len(self.poiss_models.keys())

        self.non_poiss_list_num_model_params = [d['n_params'] for d in self.non_poiss_models.values()]

        self.n_non_poiss = np.sum(self.non_poiss_list_num_model_params)
        self.n_non_poiss_models = len(self.non_poiss_models.keys())

        self.model_decompression_key = [[key, self.poiss_models[key]['log_prior']] for key in self.poiss_models.keys()]
        for key in self.non_poiss_models.keys():
            for j in range(self.non_poiss_models[key]['n_params']):
                self.model_decompression_key += [
                    [key, self.non_poiss_models[key]['log_prior'][j]]]  # False means non-poiss template

        self.poiss_model_keys = self.poiss_models.keys()

        self.n_params = len(self.model_decompression_key)

        self.set_k_max()

        self.priors = self.poiss_list_prior_range + self.non_poiss_list_prior_range
        self.configure_priors()

        print 'The number of parameters is ', self.n_params
        print 'The dynamical k_max is ', self.k_max

    def set_k_max(self):
        if self.k_max != False:
            if np.max(self.CTB_masked_compressed_summed) < self.k_max:
                self.k_max = np.max(self.CTB_masked_compressed_summed) + 1

    def set_k_max_edep(self):
        if self.k_max != False:
            self.k_max_array = []
            for en in range(self.number_energy_bins):
                if np.max(self.CTB_masked_compressed[en]) < self.k_max:
                    self.k_max_array.append(np.max(self.CTB_masked_compressed[en]) + 1)
                else:
                    self.k_max_array.append(self.k_max)
        else:
            self.k_max_array = [False for en in range(self.number_energy_bins)]

        print 'The k_max_array is ', self.k_max_array

    def convert_log(self, key, val):
        tag, log_prior = key
        if log_prior:
            return [tag, 10 ** val]
        else:
            return [tag, val]

    # @autojit
    def make_xbg(self, theta):
        theta_list = theta
        A_theta = []
        theta_PS = []
        if self.n_poiss > 0:
            A_theta = [self.convert_log(self.model_decompression_key[i], theta_list[i]) for i in range(self.n_poiss)]

        if self.n_non_poiss > 0:
            theta_PS = [self.convert_log(self.model_decompression_key[i], theta_list[i]) for i in
                        range(self.n_poiss, self.n_params)]
        else:
            theta_PS = 0

        if self.use_fixed_template:
            if self.n_poiss > 0:  # len(A_theta) > 0:
                xbg_PSF_compressed = self.fixed_template_dict_nested['summed_templates'] + np.sum(
                    map(lambda i: A_theta[i][1] * self.templates_dict_nested[A_theta[i][0]]['summed_templates'],
                        range(len(A_theta))), axis=0)
            else:
                xbg_PSF_compressed = self.fixed_template_dict_nested['summed_templates']
        else:
            xbg_PSF_compressed = np.sum(
                map(lambda i: A_theta[i][1] * self.templates_dict_nested[A_theta[i][0]]['summed_templates'],
                    range(len(A_theta))), axis=0)

            # np.sum(A_theta[i,1]*self.poiss_models[A_theta[0]]['summed_templates'] for i in range(len(A_theta)))

        return xbg_PSF_compressed, theta_PS

    # @autojit
    def sort_list(self, the_list, num):
        indices = range(0, len(the_list) + 1, num)
        new_list = [the_list[indices[i]:indices[i + 1]] for i in range(len(indices) - 1)]
        return new_list

    # @autojit
    def make_xbg_edep(self, theta):
        theta_list = theta
        A_theta = []
        theta_PS = []
        # print [ self.convert_log(self.model_decompression_key[i],theta_list[i]) for i in range(self.n_poiss) ]

        if self.n_poiss > 0:

            A_theta = self.sort_list(
                [self.convert_log(self.model_decompression_key[i], theta_list[i]) for i in range(self.n_poiss)],
                self.n_poiss / self.number_energy_bins)

        if self.n_non_poiss > 0:
            theta_PS = [self.convert_log(self.model_decompression_key[i], theta_list[i]) for i in
                        range(self.n_poiss, self.n_params)]
            if len(self.fixed_params_list) > 0:
                index = 0
                for key in self.non_poiss_models_keys:
                    if self.non_poiss_models[key]['fixed_params'] is not None:
                        for fixed in self.non_poiss_models[key]['fixed_params']:
                            theta_PS.insert(int(fixed[0]) + index, [key, fixed[1]])
                    index += self.non_poiss_models[key]['n_params_total']  # n_params

        if self.use_fixed_template:
            if self.n_poiss > 0:  # len(A_theta[0]) > 0:
                xbg_PSF_compressed = [self.fixed_template_dict_nested['templates_masked_compressed'][en] + np.sum(map(
                    lambda i: A_theta[en][i][1] *
                              self.templates_dict_nested[A_theta[en][i][0]]['templates_masked_compressed'][en],
                    range(len(A_theta[en]))), axis=0) for en in range(self.number_energy_bins)]
            else:
                xbg_PSF_compressed = [self.fixed_template_dict_nested['templates_masked_compressed'][en] for en in
                                      range(self.number_energy_bins)]
        else:
            xbg_PSF_compressed = [np.sum(map(lambda i: A_theta[en][i][1] *
                                                       self.templates_dict_nested[A_theta[en][i][0]][
                                                           'templates_masked_compressed'][en], range(len(A_theta[en]))),
                                         axis=0) for en in range(self.number_energy_bins)]

        return xbg_PSF_compressed, theta_PS

    def make_xbg_edep_exposure(self, theta):
        theta_list = theta
        A_theta = []
        theta_PS = []

        if self.n_poiss > 0:
            A_theta = self.sort_list(
                [self.convert_log(self.model_decompression_key[i], theta_list[i]) for i in range(self.n_poiss)],
                self.n_poiss / self.number_energy_bins)

        if self.n_non_poiss > 0:
            theta_PS = [self.convert_log(self.model_decompression_key[i], theta_list[i]) for i in
                        range(self.n_poiss, self.n_params)]
            if len(self.fixed_params_list) > 0:
                index = 0
                for key in self.non_poiss_models_keys:
                    if self.non_poiss_models[key]['fixed_params'] is not None:
                        for fixed in self.non_poiss_models[key]['fixed_params']:
                            theta_PS.insert(int(fixed[0]) + index, [key, fixed[1]])
                    index += self.non_poiss_models[key]['n_params_total']  # n_params

        if self.use_fixed_template:
            if self.n_poiss > 0:  # len(A_theta[0]) > 0:
                xbg_PSF_compressed = [self.fixed_template_dict_nested['templates_masked_compressed_exposure'][en] + np.sum(map(
                    lambda i: A_theta[en][i][1] *
                              self.templates_dict_nested[A_theta[en][i][0]]['templates_masked_compressed_exposure'][en],
                    range(len(A_theta[en]))), axis=0) for en in range(self.number_energy_bins)]
            else:
                xbg_PSF_compressed = [self.fixed_template_dict_nested['templates_masked_compressed_exposure'][en] for en in
                                      range(self.number_energy_bins)]
        else:
            xbg_PSF_compressed = [[np.sum(map(lambda i: A_theta[en][i][1] *
                                                       self.templates_dict_nested[A_theta[en][i][0]][
                                                           'templates_masked_compressed_exposure'][region][en], range(len(A_theta[en]))), axis=0)
                                        for region in range(self.nexp)] for en in range(self.number_energy_bins)]

        return xbg_PSF_compressed, theta_PS

    def log_like_poisson(self, theta, ndim, nparams):

        xbg_PSF_compressed = self.make_xbg(theta)[0]

        # mt.make_xbg_PSF([A_diff, A_bubs, A_iso],x_vec)

        ll = log_like.log_poisson(xbg_PSF_compressed, np.array(self.CTB_masked_compressed_summed, dtype='int32'))

        return ll

    def log_like_poisson_edep(self, theta, ndim, nparams):
        xbg_PSF_compressed = self.make_xbg_edep(theta)[0]

        ll = np.sum([log_like.log_likelihood_poissonian(xbg_PSF_compressed[i],
                                                   np.array(self.CTB_masked_compressed[i], dtype='int32')) for i in
                     range(self.number_energy_bins)])

        return ll

    def log_like_NPTF(self, theta, ndim, nparams):
        xbg_PSF_compressed, theta_PS = self.make_xbg(theta)

        self.NPTF_inst.theta_PS = list(np.vectorize(float)(np.array(theta_PS)[::, 1]))
        self.NPTF_inst.xbg_PSF_compressed = xbg_PSF_compressed
        self.NPTF_inst.set_xbg()

        ll = self.NPTF_inst.make_ll_PDF_PSF()

        del xbg_PSF_compressed, self.NPTF_inst.xbg_PSF_compressed, self.NPTF_inst.xbg_k_bins
        return ll

    def single_model_parameters_edep(self, params, energy_bin=0):
        return [params[0] / params[3 + energy_bin], params[1], params[2], params[3 + energy_bin]]

    def single_model_parameters_edep_nbreak(self, params, energy_bin=0):
        nbreak = self.nbreak
        nfix_index = (nbreak - 1) * 2 + 3
        return [params[0] / float(params[nfix_index + energy_bin])] + params[1:nfix_index - nbreak + 1] + [
            params[nfix_index + energy_bin] * np.product(params[i:nfix_index]) for i in range(nfix_index - nbreak + 1, nfix_index)] + [
                   params[nfix_index + energy_bin]]
        # return [params[0] / float(params[nfix_index + energy_bin])] + params[1:nfix_index - nbreak + 1] + [
        #     params[nfix_index + energy_bin] * params[i] for i in range(nfix_index - nbreak + 1, nfix_index)] + [
        #            params[nfix_index + energy_bin]]

    # def single_model_parameters_specify_nbreak(self, params):
    #     nbreak = self.nbreak
    #     nfix_index = (nbreak - 1) * 2 + 3
    #     return [params[0] / float(params[-1])] + params[1:]

    def model_parameters_edep(self, theta_PS, energy_bin):
        if self.n_non_poiss_models == 1:
            return self.single_model_parameters_edep(theta_PS[0:3] + [theta_PS[3 + energy_bin]])
        elif self.n_non_poiss_models == 2:
            off_set = self.non_poiss_models[self.model_decompression_key[self.n_poiss][0]]['n_params_total']
            return self.single_model_parameters_edep(
                theta_PS[0:3] + [theta_PS[3 + energy_bin]]) + self.single_model_parameters_edep(
                theta_PS[0 + off_set:3 + off_set] + [theta_PS[3 + energy_bin + off_set]])
        else:
            off_set = self.non_poiss_models[self.model_decompression_key[self.n_poiss][0]]['n_params_total']
            off_set_2 = off_set + self.non_poiss_models[self.model_decompression_key[self.n_poiss + off_set][0]][
                'n_params_total']
            return self.single_model_parameters_edep(
                theta_PS[0:3] + [theta_PS[3 + energy_bin]]) + self.single_model_parameters_edep(
                theta_PS[0 + off_set:3 + off_set] + [
                    theta_PS[3 + energy_bin + off_set]]) + self.single_model_parameters_edep(
                theta_PS[0 + off_set_2:3 + off_set_2] + [theta_PS[3 + energy_bin + off_set_2]])

    def model_parameters_edep_nbreak(self, theta_PS, energy_bin):
        nbreak = self.nbreak
        nfix_index = (nbreak - 1) * 2 + 3
        if self.n_non_poiss_models == 1:
            if self.dnds_model_array[0] == 'specify_breaks': #only works for energy_bin=0 for now
                return [theta_PS[0] / float(theta_PS[-1])] + theta_PS[1:]
            else:
                return self.single_model_parameters_edep_nbreak(
                    theta_PS[0:nfix_index] + theta_PS[nfix_index + energy_bin: nfix_index + nbreak + energy_bin])
        elif self.n_non_poiss_models == 2:
            off_set = self.non_poiss_models[self.model_decompression_key[self.n_poiss][0]]['n_params_total']
            return self.single_model_parameters_edep_nbreak(theta_PS[0:nfix_index] + theta_PS[
                                                                                     nfix_index + energy_bin: nfix_index + nbreak + energy_bin]) + self.single_model_parameters_edep_nbreak(
                theta_PS[0 + off_set:nfix_index + off_set] + theta_PS[
                                                             nfix_index + energy_bin + off_set: nfix_index + nbreak + energy_bin + off_set])
        else:
            off_set = self.non_poiss_models[self.model_decompression_key[self.n_poiss][0]]['n_params_total']
            off_set_2 = off_set + self.non_poiss_models[self.model_decompression_key[self.n_poiss + off_set][0]][
                'n_params_total']
            return self.single_model_parameters_edep_nbreak(theta_PS[0:nfix_index] + theta_PS[
                                                                                     nfix_index + energy_bin: nfix_index + nbreak + energy_bin]) + self.single_model_parameters_edep_nbreak(
                theta_PS[0 + off_set:nfix_index + off_set] + theta_PS[
                                                             nfix_index + energy_bin + off_set: nfix_index + nbreak + energy_bin + off_set]) + self.single_model_parameters_edep_nbreak(
                theta_PS[0 + off_set_2:nfix_index + off_set_2] + theta_PS[
                                                                 nfix_index + energy_bin + off_set_2: nfix_index + nbreak + energy_bin + off_set_2])

            # @jit

    def prepare_log_like_NPTF_edep(self, theta):
        xbg_PSF_compressed, theta_PS_marked = self.make_xbg_edep(theta)
        theta_PS = list(np.vectorize(float)(np.array(theta_PS_marked)[::, 1]))
        return xbg_PSF_compressed, theta_PS

    def log_like_NPTF_edep(self, theta, ndim, nparams):
        if self.use_exposure:
            xbg_PSF_compressed, theta_PS_marked = self.make_xbg_edep_exposure(theta)
        else:
            xbg_PSF_compressed, theta_PS_marked = self.make_xbg_edep(theta)
        theta_PS = list(np.vectorize(float)(np.array(theta_PS_marked)[::, 1]))

        ll = 0.0

        for energy_bin in range(self.number_energy_bins):
            self.NPTF_inst_array[energy_bin].theta_PS = self.model_parameters_edep_nbreak(theta_PS, energy_bin)
            self.NPTF_inst_array[energy_bin].xbg_PSF_compressed = xbg_PSF_compressed[energy_bin]
            self.NPTF_inst_array[energy_bin].set_xbg()

            ll += self.NPTF_inst_array[energy_bin].make_ll_PDF_PSF()

        # print 'On index = ', self.index
        # print 'theta_PS: ', theta_PS
        # print 'theta: ', np.asarray(theta)
        # print 'min: ', np.min(xbg_PSF_compressed)
        # print 'max: ',np.max(xbg_PSF_compressed)
        # print 'mean: ', np.mean(xbg_PSF_compressed)
        # print 'll: ',ll

        self.index += 1

        return ll

    def perform_scan(self, run_tag='test', edep=False, importance_nested_sampling=False, evidence_tolerance=0.5,
                     const_efficiency_mode=False,log_zero=None):
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
                        evidence_tolerance=evidence_tolerance, const_efficiency_mode=const_efficiency_mode,log_zero=log_zero)

    def lp_pc(self, params):
        return self.lp(np.array(params), self.n_params, self.n_params)

    def ll_pc(self, params):
        temp = self.ll(np.array(params), self.n_params, self.n_params)
        return temp

    def perform_scan_polychord(self, run_tag='test', edep=False, nlive='False', n_chords=1):
        self.make_dirs_for_run(run_tag, polychord=True)
        if nlive != 'False':
            self.nlive = nlive
        import pypolychord
        print 'Performing scan with polychord.  The number of parameters is ', self.n_params
        pypolychord.run(self.ll_pc, self.log_prior_polychord, self.n_params, n_live=self.nlive, n_chords=n_chords,
                        output_basename=self.chains_dir_for_run)

        # self.load_scan(run_tag = run_tag, edep = edep)

    def perform_scan_minuit(self, run_tag='test', edep=False, print_level=1, **kwargs):
        self.run_tag = run_tag
        self.prepare_for_minuit(print_level=print_level, **kwargs)
        self.m.migrad()

    def save_minuit(self):
        self.matrix = self.m.matrix()
        self.fitarg = self.m.fitarg.copy()
        self.fitarg['ll_value'] = self.m.fval
        # print self.fitarg
        import json
        with open(self.plots_dir_for_run + 'fitarg.txt', 'wb') as f:
            json.dump(self.fitarg, f)

        with open(self.plots_dir_for_run + 'll.txt', 'wb') as f:
            json.dump(self.fitarg['ll_value'], f)
        # np.savetxt(self.plots_dir_for_run + 'fitarg.txt',self.fitarg,delimiter=" ", fmt="%s")

        with open(self.plots_dir_for_run + 'covariance.txt', 'wb') as f: json.dump(self.matrix, f)

        self.calculate_samples_minuit()
        np.savetxt(self.chains_dir_for_run + '.txt', self.samples)

    def load_minuit(self, run_tag='test', print_level=1, **kwargs):
        self.run_tag = run_tag
        self.prepare_for_minuit(print_level=print_level, **kwargs)
        print 'fitarg file for loading minuit: ', self.plots_dir_for_run + 'fitarg.txt'
        import json
        with open(self.plots_dir_for_run + 'fitarg.txt', 'r') as f:
            self.fitarg = json.load(f)

        print 'fitarg is ', self.fitarg

        with open(self.plots_dir_for_run + 'll.txt', 'r') as f:
            self.fitarg['ll_value'] = json.load(f)

        with open(self.plots_dir_for_run + 'covariance.txt', 'r') as f:
            self.matrix = json.load(f)

        self.m.fitarg = self.fitarg.copy()

        self.samples = np.loadtxt(self.chains_dir_for_run + '.txt')

    def plot_ll_profile(self, param='All'):
        if param == 'All':
            for param in self.params:
                self.m.draw_mnprofile(param)
                plt.savefig(self.plots_dir_for_run + param + '.pdf')
                plt.close()
        else:
            self.m.draw_mnprofile(param)
            plt.savefig(self.plots_dir_for_run + param + '.pdf')
            plt.close()

    def prepare_for_minuit(self, print_level=1, **kwargs):
        self.prepare_minuit_params()
        self.make_dirs_for_run(self.run_tag, polychord=True)
        from iminuit import Minuit, describe, Struct, util
        import config.iminuit_functions as imf
        self.ll_minuit = imf.call_ll(self.n_params, self.ll_pc, self.params)
        kwdargs = {'limit_' + self.params[i]: (self.theta_min[i], self.theta_max[i]) for i in range(self.n_params)}
        other_kwargs = {'print_level': print_level}
        z = kwdargs.copy()
        z.update(other_kwargs)
        z.update(kwargs)
        self.m = Minuit(self.ll_minuit, **z)

    def prepare_minuit_params(self):
        self.poiss_params = np.array([self.poiss_models[key]['model_tag'] for key in self.poiss_model_keys])
        if len(self.non_poiss_models.values()) > 0:
            self.non_poiss_params = self.flatten(np.array([val['model_tag'] for val in self.non_poiss_models.values()]))
            self.params = np.array(list(self.poiss_params) + list(self.non_poiss_params))
        else:
            self.params = np.array(self.poiss_params)

    def load_scan(self, run_tag='test', edep=False, polychord=False, minuit=False, *args, **kwargs):
        if self.already_loaded:
            pass
        else:
            self.run_tag = run_tag
            self.make_dirs_for_run(run_tag)

            print 'In load_scan: polychord, minuit:', polychord, minuit

            if polychord:
                pass    
            elif minuit:
                self.load_minuit(run_tag=run_tag, *args, **kwargs)
            else:
                self.a = pymultinest.Analyzer(n_params=self.n_params, outputfiles_basename=self.chains_dir_for_run)
                self.s = self.a.get_stats()
                self.chain_file = self.chains_dir_for_run + '/post_equal_weights.dat'
                self.samples = np.array(np.loadtxt(self.chain_file)[:, :-1])
            # self.samples
            if not minuit:
                self.calculate_medians()
            else:
                self.calculate_best_fit_minuit()

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
        # self.save_dictionary()

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

    def calculate_best_fit_minuit(self):
        # self.medians = [self.m.values[param] for param in self.params]
        self.medians = [self.m.fitarg[param] for param in self.params]

        # Lina: Calculate errors for minuit

    def calculate_samples_minuit(self, n_eval=1000):
        self.calculate_best_fit_minuit()
        # covariance = np.array(np.loadtxt(self.plots_dir_for_run + 'covariance.txt'))
        self.samples = np.random.multivariate_normal(self.medians, self.matrix, n_eval)

    def calculate_norms(self, edep):
        if edep:
            # self.calculate_norms_edep()
            self.n_poiss_base = len(self.poiss_models.keys())
            self.poiss_comp_numbers = {self.poiss_model_keys[i]: i for i in range(self.n_poiss)}.copy()
            self.test_array = self.poiss_comp_numbers.copy()
            # print 'comp #s: ', self.poiss_comp_numbers
            self.poiss_list_is_log_prior = np.array(
                [self.poiss_models[key]['log_prior'] for key in self.poiss_model_keys])  # *self.number_energy_bins )
            self.poiss_params = np.array([self.poiss_models[key]['model_tag'] for key in self.poiss_model_keys])
        else:
            self.poiss_comp_numbers = {self.poiss_models.keys()[i]: i for i in range(self.n_poiss)}
            self.poiss_list_is_log_prior = np.array([val['log_prior'] for val in self.poiss_models.values()])
            self.poiss_params = np.array([val['model_tag'] for val in self.poiss_models.values()])
        self.non_poiss_params = self.flatten(np.array([val['model_tag'] for val in self.non_poiss_models.values()]))
        self.params = np.array(list(self.poiss_params) + list(self.non_poiss_params))
        # self.poiss_list_is_log_prior = np.array( [val['log_prior'] for val in self.poiss_models.values() ] )
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

    def save_dictionary(self):
        self.the_dict = {k: v for k, v in self.__dict__.items() if not (k.startswith('__') and k.endswith('__'))}
        for key, item in self.the_dict.items():
            if inspect.ismethod(item) or isinstance(item, types.InstanceType) or isinstance(item, types.FunctionType):
                del self.the_dict[key]

        self.keys = np.array(list(self.the_dict.keys()))
        self.values = np.array(list(self.the_dict.values()))
        np.savez_compressed(self.dict_dir + self.run_tag + '-keys', self.keys)
        np.savez_compressed(self.dict_dir + self.run_tag + '-values', self.values)
