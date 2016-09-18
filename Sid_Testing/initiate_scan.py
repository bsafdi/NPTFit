import healpy as hp
import numpy as np

import config.NPTF_models as NPTF_models
from config.config_file import config
from config.bayesian_scan import bayesian_scan

class bayesian_scan_NPTF(bayesian_scan):
    """ The file that sets up and executes the scan """

    def __init__(self, *args, **kwargs):

        bayesian_scan.__init__(self, *args, **kwargs)

    def initiate_ps_from_f_ary(self, f_ary, df_rho_div_f_ary, nbreak=1):

        if self.use_exposure:
            self.initiate_ps_from_f_ary_exposure(sigma_psf_deg, nbreak, **kwargs)
        else:
            self.nbreak = nbreak
            self.configure_for_edep_scan()

            self.NPT_compressed_ary = [self.templates_dict_nested[self.non_poiss_models.keys()[i]]['summed_templates'] for i in range(self.number_ps_models)]

            self.NPTF_inst_array = [
                    NPTF_models.NPTF_ps(self.NPT_compressed_ary, self.nexp, self.psf_dir, self.CTB_masked_compressed[en],
                                          k_max=self.k_max_array[en], Sc=self.Sc, nside=self.nside, nbreak=self.nbreak,
                                          f_ary=f_ary[en], df_rho_div_f_ary=df_rho_div_f_ary[en]) for en in range(self.number_energy_bins)]

            self.ll = self.log_like_NPTF_edep

    # def initiate_ps_from_f_ary_exposure(self, sigma_psf_deg, nbreak=1, **kwargs):
    #     self.nbreak = nbreak
    #     self.exposure_pixels_compressed()
    #     self.configure_for_edep_scan()
    #     self.use_exposure = True

    #     if 'f_ary' and 'df_rho_div_f_ary' in kwargs.keys():
    #         f_ary = kwargs['f_ary']
    #         df_rho_div_f_ary = kwargs['df_rho_div_f_ary']
    #         kwarg_array = [{'f_ary': f_ary[en], 'df_rho_div_f_ary': df_rho_div_f_ary[en]} for en in
    #                        range(self.number_energy_bins)]
    #     else:
    #         kwarg_array = [{'sigma_psf_deg': sigma_psf_deg[en]} for en in range(self.number_energy_bins)]

    #     if 'PS_1_dist' in kwargs.keys():
    #         PS_1_dist = kwargs['PS_1_dist']  # summed, but not masked or compressed
    #         self.PS_1_dist_compressed_exposure = self.return_masked_exposure([PS_1_dist], compressed=True)[0]
    #         # TODO: Check the statement right before
    #     elif 'PS_1_dist_compressed' in kwargs.keys():
    #         self.PS_1_dist_compressed_exposure = kwargs['PS_1_dist_compressed']
    #     else:
    #         # The hackiest worst way to flip the energy and the number of regions
    #         self.PS_1_dist_compressed_exposure = self.templates_dict_nested[self.non_poiss_models.keys()[0]]['flux_summed_templates_exposure']


    #     self.NPTF_inst_array = [
    #             NPTF_models.NPTF_1_ps(self.PS_1_dist_compressed_exposure, self.nexp, self.psf_dir,
    #                                   np.array(self.CTB_masked_compressed_exposure)[::,en],
    #                                   k_max=self.k_max_array[en], Sc=self.Sc, nside=self.nside, nbreak=self.nbreak,
    #                                   use_exposure=True, exposure=self.exposure_means_list,
    #                                   exposure_mean=self.exposure_mean, **kwarg_array[en])
    #             for en in range(self.number_energy_bins)]


    #     # self.NPTF_inst.set_ps()

    #     self.ll = self.log_like_NPTF_edep

    def initiate_poissonian_edep(self):
        self.configure_for_edep_scan()
        self.ll = self.log_like_poisson_edep
