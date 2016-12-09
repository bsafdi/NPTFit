# PSF-making tests

import sys
sys.path.append("../NPTFit")

import numpy as np

from NPTFit import psf_correction as pc # Module for determining the PSF correction

def test_psf():

    pc_inst = pc.PSFCorrection(psf_sigma_deg=0.1812,n_psf=5000,n_pts_per_psf=100,num_f_bins=20,f_trunc=0.1,nside=64)
    f_ary_gauss = pc_inst.f_ary
    df_rho_div_f_ary_gauss = pc_inst.df_rho_div_f_ary

    # Define parameters that specify the Fermi-LAT PSF at 2 GeV
    fcore = 0.748988248179
    score = 0.428653790656
    gcore = 7.82363229341
    stail = 0.715962650769
    gtail = 3.61883748683
    spe = 0.00456544262478

    # Define the full PSF in terms of two King functions
    def king_fn(x, sigma, gamma):
        return 1./(2.*np.pi*sigma**2.)*(1.-1./gamma)*(1.+(x**2./(2.*gamma*sigma**2.)))**(-gamma)

    def Fermi_PSF(r):
        return fcore*king_fn(r/spe,score,gcore) + (1-fcore)*king_fn(r/spe,stail,gtail)

    # Modify the relevant parameters in pc_inst and then make or load the PSF
    pc_inst = pc.PSFCorrection(delay_compute=True)
    pc_inst.psf_r_func = lambda r: Fermi_PSF(r)
    pc_inst.sample_psf_max = 10.*spe*(score+stail)/2.
    pc_inst.psf_samples = 10000
    pc_inst.psf_tag = 'Fermi_PSF_2GeV'
    pc_inst.make_or_load_psf_corr()

    # Extract f_ary and df_rho_div_f_ary as usual
    f_ary_9 = pc_inst.f_ary
    df_rho_div_f_ary_9 = pc_inst.df_rho_div_f_ary
