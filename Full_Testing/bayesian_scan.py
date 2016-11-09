import numpy as np
import pymultinest
from mpi4py import MPI
import copy
import NPLL as log_like
from config_file import config


class bayesian_scan(config):
    """ The file that sets up and executes the scan """

    def __init__(self, nlive=100, *args, **kwargs):
        config.__init__(self, *args, **kwargs)

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


        self.non_poiss_models[template_name] = {'prior_range': prior_range, 'log_prior': log_prior_list,
                                                'model_tag': model_tag, 'n_params': n_params, 'poiss': False,
                                                'dnds_model': dnds_model, 'fixed_params': fixed_params}

        self.non_poiss_models[template_name]['n_params_total'] = self.non_poiss_models[template_name]['n_params']
        if fixed_params is not None:
            self.non_poiss_models[template_name]['n_params_total'] += len(
                self.non_poiss_models[template_name]['fixed_params'])

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

        self.configure_priors(number_energy_bins=self.number_energy_bins)

        print 'The number of parameters is ', self.n_params

    def convert_log(self, key, val):
        tag, log_prior = key
        if log_prior:
            return [tag, 10 ** val]
        else:
            return [tag, val]

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

        xbg_PSF_compressed = np.sum(
            map(lambda i: A_theta[i][1] * self.templates_dict_nested[A_theta[i][0]]['summed_templates'],
                range(len(A_theta))), axis=0)


        return xbg_PSF_compressed, theta_PS

    def sort_list(self, the_list, num):
        indices = range(0, len(the_list) + 1, num)
        new_list = [the_list[indices[i]:indices[i + 1]] for i in range(len(indices) - 1)]
        return new_list

    def make_xbg_edep(self, theta):
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

        xbg_PSF_compressed = [[np.sum(map(lambda i: A_theta[en][i][1] *
                                                    self.templates_dict_nested[A_theta[en][i][0]][
                                                        'templates_masked_compressed_exposure'][region][en], range(len(A_theta[en]))), axis=0)
                                    for region in range(self.nexp)] for en in range(self.number_energy_bins)]

        return xbg_PSF_compressed, theta_PS

    def log_like_poisson_edep(self, theta, ndim, nparams):
        xbg_PSF_compressed = self.make_xbg_edep(theta)[0]

        ll = np.sum([log_like.log_likelihood_poissonian(xbg_PSF_compressed[i],
                                                   np.array(self.CTB_masked_compressed[i], dtype='int32')) for i in
                     range(self.number_energy_bins)])

        return ll

    def single_model_parameters_edep(self, params, energy_bin=0):
        return [params[0] / params[3 + energy_bin], params[1], params[2], params[3 + energy_bin]]

    def single_model_parameters_edep_nbreak(self, params, energy_bin=0):
        nbreak = self.nbreak
        nfix_index = (nbreak - 1) * 2 + 3
        return [params[0] / float(params[nfix_index + energy_bin])] + params[1:nfix_index - nbreak + 1] + [
            params[nfix_index + energy_bin] * np.product(params[i:nfix_index]) for i in range(nfix_index - nbreak + 1, nfix_index)] + [
                   params[nfix_index + energy_bin]]

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


    def log_like_NPTF_edep(self, theta, ndim, nparams):
        xbg_PSF_compressed, theta_PS_marked = self.make_xbg_edep_exposure(theta)
        theta_PS = list(np.vectorize(float)(np.array(theta_PS_marked)[::, 1]))

        ll = 0.0

        for energy_bin in range(self.number_energy_bins):
            self.NPTF_inst_array[energy_bin].theta_PS = self.model_parameters_edep_nbreak(theta_PS, energy_bin)
            self.NPTF_inst_array[energy_bin].xbg_PSF_compressed = xbg_PSF_compressed[energy_bin]

            ll += self.NPTF_inst_array[energy_bin].make_ll_PDF_PSF()

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
                [self.poiss_models[key]['log_prior'] for key in self.poiss_model_keys])  # *self.number_energy_bins )
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
