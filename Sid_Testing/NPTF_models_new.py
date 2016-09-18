import numpy as np
import pulsars.likelihood_psf_pdep_new as llpsfn
import pulsars.likelihood_psf_pdep_2break as llpsfb
import pulsars.likelihood_psf_pdep_3break as llpsfb3
from config.NPTF import NPTF


class NPTF_ps(NPTF):
    def __init__(self, PS_dist_compressed, nexp, *args, **kwargs):

        self.nexp = nexp
        self.PS_dist_compressed = PS_dist_compressed

        NPTF.__init__(self, *args, **kwargs)
        # if self.use_exposure:
        #     self.set_ps_exposure()
        #     if self.nbreak == 1:
        #         self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_1break_exposure
        #     elif self.nbreak == 2:
        #         self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_2break_exposure
        #     elif self.nbreak == 3:
        #         self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_3break_exposure
        # else:
        self.set_ps()
        self.make_ll_PDF_PSF = self.make_ll_PDF_PSF_1break

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
        ll = call general exposure code

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

            # modify for arbitrary theta PS
            theta_PS = [self.theta_PS[0] * self.exposure_mean / self.exposure_means_list[i]] + list(
                self.theta_PS[1:3]) + \
                       [self.theta_PS[-1] * self.exposure_means_list[i] / self.exposure_mean]

            ll += call general exposure code

        return ll
