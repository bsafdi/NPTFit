# Poissonian and non-Poissonian scan tests

import sys
sys.path.append("../NPTFit")

import numpy as np

from NPTFit import nptfit # module for performing scan
from NPTFit import psf_correction as pc # module for determining the PSF correction

def test_scan_non_poiss():
    n = nptfit.NPTF(tag='Test_NPoiss')

    fermi_data = np.array([2, 1, 1, 1, 4, 10]).astype(np.int32)
    fermi_exposure = np.array([1., 1., 1., 2., 2., 2.])
    n.load_data(fermi_data, fermi_exposure)

    analysis_mask = np.array([False, False, False, False, False, True])

    dif = np.array([1., 2., 3., 4., 5., 6.])
    iso = np.array([1., 1., 1., 1., 1., 1.])

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
    n.load_scan()

    n = nptfit.NPTF(tag='Test_NPoiss')
    n.load_data(fermi_data, fermi_exposure)
    n.add_template(iso, 'iso')
    n.add_non_poiss_model('iso',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$S_{b1}$'],
                          [[0,1],[2.05,2.5],[-2.,-1.5],[1.,4.]])
    n.configure_for_scan()
    n.perform_scan(nlive=50)

    n = nptfit.NPTF(tag='Test_NPoiss')
    n.load_data(fermi_data, fermi_exposure)
    n.add_template(iso, 'iso')
    n.add_non_poiss_model('iso',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$S_{b1}$'],
                          [[0,1]],
                          units='flux',fixed_params=[[1,2.],[2,-2.],[3,100.]])
    n.configure_for_scan(nexp=len(dif)+1)
    n.perform_scan(nlive=50)

    n = nptfit.NPTF(tag='Test_NPoiss')
    n.load_data(fermi_data, fermi_exposure)
    n.add_template(iso, 'iso')
    n.add_non_poiss_model('iso',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$n_3','$S_{b1}$','$S_{b2}$'],
                          [[0,1],[0,10],[0,1.]],
                          units='flux',fixed_params=[[1,2.],[2,-2.],[3,1.0]],
                          dnds_model='specify_relative_breaks')
    n.configure_for_scan()
    n.perform_scan(nlive=50)

    n = nptfit.NPTF(tag='Test_NPoiss')
    n.load_data(fermi_data, fermi_exposure)
    n.add_template(iso, 'iso')
    n.add_non_poiss_model('iso',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$n_3','$S_{b1}$','$S_{b2}$'],
                          [[0,200]],
                          units='flux',fixed_params=[[0,1.],[1,2.],[2,1.],[3,-2],[5,250]])
    n.configure_for_scan()
    n.perform_scan(nlive=50)

    n = nptfit.NPTF(tag='Test_NPoiss')
    n.load_data(fermi_data, fermi_exposure)
    n.add_template(iso, 'iso')
    n.add_non_poiss_model('iso',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$n_3','$S_{b1}$','$S_{b2}$'],
                          [[200,400]],
                          units='flux',dnds_model='specify_relative_breaks',
                          fixed_params=[[0,1.],[1,2.],[2,1.],[3,-2],[4,150]])
    n.configure_for_scan()
    n.perform_scan(nlive=50)

    n = nptfit.NPTF(tag='Test_NPoiss')
    n.load_data(fermi_data, fermi_exposure)
    n.add_template(iso, 'iso')
    n.add_non_poiss_model('iso',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$n_3','$S_{b1}$','$S_{b2}$'],
                          [[0.,1.]],[True]
                          units='flux',fixed_params=[[0,1.],[1,2.],[2,1.],[3,-2],[5,250]])
    n.configure_for_scan()
    n.perform_scan(nlive=50)

    n = nptfit.NPTF(tag='Test_NPoiss')
    n.load_data(fermi_data, fermi_exposure)
    n.add_template(iso, 'iso')
    n.add_non_poiss_model('iso',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$n_3','$S_{b1}$','$S_{b2}$'],
                          [[0.1,1]],[True]
                          units='flux',dnds_model='specify_relative_breaks',
                          fixed_params=[[0,1.],[1,2.],[2,1.],[3,-2],[4,150]])
    n.configure_for_scan()
    n.perform_scan(nlive=50)

    n = nptfit.NPTF(tag='Test_NPoiss')
    n.load_data(fermi_data, fermi_exposure)
    n.add_template(iso, 'iso')
    n.add_template(dif, 'dif')
    n.add_non_poiss_model('iso',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$S_{b1}$'],
                          [[0,1]],
                          units='flux',fixed_params=[[1,2.],[2,-2.],[3,100.]])
    n.add_non_poiss_model('dif',
                          ['$A^\mathrm{ps}_\mathrm{iso}$','$n_1$','$n_2$','$S_{b1}$'],
                          [[0,1]],
                          units='flux',fixed_params=[[1,2.],[2,-2.],[3,100.]])
    n.configure_for_scan()
    n.perform_scan(nlive=50)

def test_scan_poiss():
    n = nptfit.NPTF(tag='Test_Poiss')

    fermi_data = np.array([2, 1, 1, 1, 4, 10]).astype(np.int32)
    fermi_exposure = np.array([1., 1., 1., 2., 2., 2.])
    n.load_data(fermi_data, fermi_exposure)

    analysis_mask = np.array([False, False, False, False, False, True])
    n.load_mask(analysis_mask)

    dif = np.array([1., 2., 3., 4., 5., 6.])
    iso = np.array([1., 1., 1., 1., 1., 1.])

    n.add_template(dif, 'dif')
    n.add_template(iso, 'iso')
    n.add_template(iso, 'iso_f', units='flux')
    n.add_template(iso, 'iso_PS', units='PS')

    n.add_poiss_model('dif', '$A_\mathrm{dif}$', [0, 30], False)
    n.add_poiss_model('iso', '$A_\mathrm{iso}$', [0, 5], False)

    n.configure_for_scan()

    n.perform_scan(nlive=50)

    n.load_scan()

    n = nptfit.NPTF(tag='Test_Poiss')

    n.load_data(fermi_data, fermi_exposure)

    n.add_template(dif, 'dif')

    n.add_poiss_model('dif', '$A_\mathrm{dif}$', [0, 30], False)

    n.configure_for_scan()

    n.perform_scan(nlive=50)

    n.load_scan()
