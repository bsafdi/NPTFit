# Poissonian and non-Poissonian scan tests

import sys
sys.path.append("../NPTFit")

import numpy as np

from NPTFit import nptfit # module for performing scan
from NPTFit import create_mask as cm # module for creating the mask
from NPTFit import psf_correction as pc # module for determining the PSF correction


def test_scan_non_poiss():
    n = nptfit.NPTF(tag='Test_NPoiss')

    fermi_data = np.load('fermi_data/fermidata_counts.npy').astype(np.int32)
    fermi_exposure = np.load('fermi_data/fermidata_exposure.npy')
    n.load_data(fermi_data, fermi_exposure)

    analysis_mask = cm.make_mask_total(mask_ring = True, inner = 0, outer = 5, ring_b = 90, ring_l = 0)
    n.load_mask(analysis_mask)

    dif = np.load('fermi_data/template_dif.npy')
    iso = np.load('fermi_data/template_iso.npy')

    n.add_template(dif, 'dif')
    n.add_template(iso, 'iso')

    n.add_poiss_model('iso','$A_\mathrm{iso}$', False, fixed=True, fixed_norm=1.47)
    n.add_non_poiss_model('iso',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$S_b$'],
                          [[-6,1],[2.05,30],[-2,1.95]],
                          [True,False,False],
                          fixed_params = [[3,22.]])

    pc_inst = pc.PSFCorrection(psf_sigma_deg=0.1812)
    f_ary, df_rho_div_f_ary = pc_inst.f_ary, pc_inst.df_rho_div_f_ary

    n.configure_for_scan(f_ary=f_ary, df_rho_div_f_ary=df_rho_div_f_ary, nexp=1)

    n.perform_scan(nlive=50)

    n.load_scan()


def test_scan_poiss():
    n = nptfit.NPTF(tag='Test_Poiss')

    fermi_data = np.load('fermi_data/fermidata_counts.npy').astype(np.int32)
    fermi_exposure = np.load('fermi_data/fermidata_exposure.npy')
    n.load_data(fermi_data, fermi_exposure)

    analysis_mask = cm.make_mask_total(mask_ring = True, inner = 0, outer = 5, ring_b = 90, ring_l = 0)
    n.load_mask(analysis_mask)

    dif = np.load('fermi_data/template_dif.npy')
    iso = np.load('fermi_data/template_iso.npy')

    n.add_template(dif, 'dif')
    n.add_template(iso, 'iso')

    n.add_poiss_model('dif', '$A_\mathrm{dif}$', [0, 30], False)
    n.add_poiss_model('iso', '$A_\mathrm{iso}$', [0, 5], False)

    n.configure_for_scan()

    n.perform_scan(nlive=50)

    n.load_scan()
