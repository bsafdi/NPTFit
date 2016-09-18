import numpy as np
from NPLL import log_like


class NPTF_ps(NPTF):
    def __init__(self, PS_dist_compressed):
        self.PS_dist_compressed = PS_dist_compressed

    def make_ll_PDF_PSF(self):
        ll = log_like(self.xbg_PSF_compressed, np.array(self.theta_PS), self.f_ary,
                                                       self.df_rho_div_f_ary, self.PS_dist_compressed,
                                                       np.array(self.CTB_map_compressed, dtype='int32'))
        return ll
