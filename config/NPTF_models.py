import numpy as np
import pulsars.likelihood_psf_pdep_new as llpsfn
import pulsars.likelihood_psf_pdep_2break as llpsfb
import pulsars.likelihood_psf_pdep_3break as llpsfb3
from config.NPTF import NPTF


class NPTF_1_ps(NPTF):
    def __init__(self, PS_dist_compressed, nexp, *args, **kwargs):

        self.nexp = nexp
        self.PS_dist_compressed = PS_dist_compressed

        NPTF.__init__(self, *args, **kwargs)
        if self.use_exposure:
            self.set_ps_exposure()
            if self.nbreak == 1:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_1break_exposure
            elif self.nbreak == 2:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_2break_exposure
            elif self.nbreak == 3:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_3break_exposure
        else:
            self.set_ps()
            if self.nbreak == 1:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_1break
            elif self.nbreak == 2:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_2break
            elif self.nbreak == 3:
                self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_3break

    def set_ps(self):

        self.PS_dist_compressed = self.PS_dist_compressed[self.k_max_where]

    def set_ps_exposure(self):

        for i in range(self.nexp):
            self.PS_dist_compressed[i] = np.array(self.PS_dist_compressed[i])[np.array(self.k_max_where_exposure[i])]

    def make_ll_PDF_PSF_1break(self):
        """
        Computes the log of the likelihood
        :return:log (likelihood) of apparently just 1 Point source template
        """
        # TODO: On ipython, there is a type matching error here, but it runs generally. Fix that!
        ll = llpsfn.log_nu_k_ary_PSF_exact_1_PS(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                self.df_rho_div_f_ary, self.PS_dist_compressed,
                                                np.array(self.CTB_map_compressed, dtype='int32'), Sc=self.Sc)

        return ll

    def make_ll_PDF_PSF_1break_exposure(self):
        """
        Computes the log of the likelihood
        :param self.nexp: number of exposure regions. By default taken to be 1 which should give the same result as before
        :return:log (likelihood) of apparently just 1 Point source template
        """
        # TODO: On ipython, there is a type matching error here, but it runs generally. Fix that!

        ll = 0.0

        for i in range(self.nexp):

            theta_PS = [self.theta_PS[0] * self.exposure_mean / self.exposure_means_list[i]] + list(
                self.theta_PS[1:3]) + \
                       [self.theta_PS[-1] * self.exposure_means_list[i] / self.exposure_mean]

            ll += llpsfn.log_nu_k_ary_PSF_exact_1_PS(self.xbg_PSF_compressed[i], np.array(theta_PS), self.f_ary,
                                                     self.df_rho_div_f_ary, self.PS_dist_compressed[i],
                                                     np.array(self.CTB_map_compressed[i], dtype='int32'), Sc=self.Sc)

        return ll


    def make_ll_PDF_PSF_2break(self):
        ll = llpsfb.log_nu_k_ary_PSF_exact_1_PS_2break(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                       self.df_rho_div_f_ary, self.PS_dist_compressed,
                                                       np.array(self.CTB_map_compressed, dtype='int32'))
        return ll

    def make_ll_PDF_PSF_2break_exposure(self):
        """
        Computes the log of the likelihood
        :return:log (likelihood) of apparently just 1 Point source template with 2 breaks
        """

        ll = 0.0
        for i in range(self.nexp):
            ####BS: fixed the line below to exposure correct both breaks!! make sure this works properly.
            theta_PS = [self.theta_PS[0] * self.exposure_mean / self.exposure_means_list[i]] + list(
                self.theta_PS[1:4]) + \
                       [self.theta_PS[-2] * self.exposure_means_list[i] / self.exposure_mean,self.theta_PS[-1] * self.exposure_means_list[i] / self.exposure_mean]

            ll += llpsfb.log_nu_k_ary_PSF_exact_1_PS_2break(self.xbg_PSF_compressed[i], np.array(theta_PS), self.f_ary,
                                                            self.df_rho_div_f_ary, self.PS_dist_compressed[i],
                                                            np.array(self.CTB_map_compressed[i], dtype='int32'))

        return ll

    def make_ll_PDF_PSF_3break(self):
        ll = llpsfb3.log_nu_k_ary_PSF_exact_1_PS_3break(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                       self.df_rho_div_f_ary, self.PS_dist_compressed,
                                                       np.array(self.CTB_map_compressed, dtype='int32'))
        return ll

    def make_ll_PDF_PSF_3break_exposure(self):
        """
        Computes the log of the likelihood
        :return:log (likelihood) of apparently just 1 Point source template with 2 breaks
        """

        ll = 0.0
        for i in range(self.nexp):
            ####BS: fixed the line below to exposure correct both breaks!! make sure this works properlyt.
            theta_PS = [self.theta_PS[0] * self.exposure_mean / self.exposure_means_list[i]] + list(
                self.theta_PS[1:5]) + \
                       [self.theta_PS[-3] * self.exposure_means_list[i] / self.exposure_mean, self.theta_PS[-2] * self.exposure_means_list[i] / self.exposure_mean,self.theta_PS[-1] * self.exposure_means_list[i] / self.exposure_mean]

            ll += llpsfb3.log_nu_k_ary_PSF_exact_1_PS_3break(self.xbg_PSF_compressed[i], np.array(theta_PS), self.f_ary,
                                                            self.df_rho_div_f_ary, self.PS_dist_compressed[i],
                                                            np.array(self.CTB_map_compressed[i], dtype='int32'))

        return ll
