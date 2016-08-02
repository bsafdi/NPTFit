import numpy as np
import config.compute_PSF as CPSF

class NPTF:
    def __init__(self, psf_dir, CTB_map_compressed, k_max='False', Sc='False', nside=128, norm_arrays='False',
                 num_f_bins=10, n_ps=10000, n_pts_per_psf=1000, nbreak=1, use_exposure=False, **kwargs):  
        self.nbreak = nbreak
        self.psf_dir = psf_dir
        self.nside = nside
        self.CTB_map_compressed = CTB_map_compressed
        self.num_f_bins = num_f_bins
        self.n_ps = n_ps
        self.n_pts_per_psf = n_pts_per_psf
        self.use_exposure = use_exposure
        if use_exposure:
            self.set_k_max_exposure(k_max)
        else:
            self.set_k_max(k_max)
        self.set_Sc(Sc=Sc)
        self.npixROI = len(self.CTB_map_compressed)
        # TODO: I don't think this is right since the dimensions of (CTB_map_compressed is (1,6000))

        self.xbg_PSF_compressed = []
        self.theta_PS = []  # This is the wrong type!!! This should end up being an array!!! TODO: Fix these types!!!

        self.pre_load_iso = False
        self.exposure_means_list = [1.0]
        self.exposure_mean = 1.0
        self.sort_kwargs_psf(kwargs)
        self.sort_exposure(kwargs)

    def sort_kwargs_psf(self, kwargs):
        if 'sigma_psf_deg' in kwargs.keys():
            self.sigma_PSF_deg = kwargs['sigma_psf_deg']
            self.make_f_ary_gaussian()
        elif 'f_ary' and 'df_rho_div_f_ary' in kwargs.keys():
            self.f_ary = kwargs['f_ary']
            self.df_rho_div_f_ary = kwargs['df_rho_div_f_ary']
        else:
            print 'Incorrect PSF info!  Stoping the scan.  But, the kwargs are ', kwargs
            sys.exit()

    def sort_exposure(self, kwargs):
        if 'exposure' in kwargs.keys():
            self.exposure_means_list = kwargs['exposure']
        if 'exposure_mean' in kwargs.keys():
            self.exposure_mean = kwargs['exposure_mean']

    def set_k_max(self,k_max):
        if k_max == 'False':
            self.k_max = max(self.CTB_map_compressed)
        else:
            self.k_max = k_max

        self.k_max_where = np.where(self.CTB_map_compressed < self.k_max)
        self.CTB_map_compressed = self.CTB_map_compressed[self.k_max_where]

    def set_k_max_exposure(self, k_max):

        if k_max == 'False':
            # TODO: Fix this back
            print "I got to k_max == False in set_kmax_exposure in NPTF.py and it has not been fixed yet"
            self.k_max = max([max(self.CTB_map_compressed[i]) for i in range(len(self.CTB_map_compressed))])
        else:
            self.k_max = k_max

        self.k_max_where_exposure = []

        for i in range(len(self.CTB_map_compressed)):
            self.k_max_where_exposure += [np.where(self.CTB_map_compressed[i] < self.k_max)[0]]
            print np.shape(self.k_max_where_exposure)
            self.CTB_map_compressed[i] = self.CTB_map_compressed[i][self.k_max_where_exposure[i]]

    def set_Sc(self, Sc='False'):
        if Sc == 'False':
            self.Sc = float(100000.0)  # fudge for now
        else:
            self.Sc = float(Sc)

    def set_xbg_exposure(self):

        for i in range(len(self.xbg_PSF_compressed)):
            self.xbg_PSF_compressed[i] = self.xbg_PSF_compressed[i][self.k_max_where_exposure[i]]

    def set_xbg(self):
        if self.use_exposure:
            self.set_xbg_exposure()
        else:
            self.xbg_PSF_compressed = self.xbg_PSF_compressed[self.k_max_where]

    def make_f_ary_gaussian(self):
        self.CPSF_inst = CPSF.compute_PSF(self.sigma_PSF_deg, self.psf_dir, nside=self.nside,
                                          num_f_bins=10)  # num_f_bins hardcoded for now

        self.f_ary, self.df_rho_div_f_ary = self.CPSF_inst.f_ary, self.CPSF_inst.df_rho_div_f_ary
