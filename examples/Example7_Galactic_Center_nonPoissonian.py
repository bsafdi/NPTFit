# This file is called by Example7_Galactic_Center_Batch.batch
# The scan performs a run over the inner galaxy

# NB: this example makes use of the Fermi Data, which needs to already be installed. See Example 1 for details.

import numpy as np

from NPTFit import nptfit # module for performing scan
from NPTFit import create_mask as cm # module for creating the mask
from NPTFit import psf_correction as pc # module for determining the PSF correction

n = nptfit.NPTF(tag='GCE_Example')

fermi_data = np.load('fermi_data/fermidata_counts.npy')
fermi_exposure = np.load('fermi_data/fermidata_exposure.npy')
n.load_data(fermi_data, fermi_exposure)

pscmask=np.array(np.load('fermi_data/fermi_pscmask_cons.npy'), dtype=bool)
analysis_mask = cm.make_mask_total(band_mask = True, band_mask_range = 2,
                                   mask_ring = True, inner = 0, outer = 30,
                                   custom_mask = pscmask)
n.load_mask(analysis_mask)

dif = np.load('fermi_data/template_dif.npy')
iso = np.load('fermi_data/template_iso.npy')
bub = np.load('fermi_data/template_bub.npy')
gce = np.load('fermi_data/template_gce.npy')
dsk = np.load('fermi_data/template_dsk.npy')

n.add_template(dif, 'dif')
n.add_template(iso, 'iso')
n.add_template(bub, 'bub')
n.add_template(gce, 'gce')
n.add_template(dsk, 'dsk')

n.add_poiss_model('dif', '$A_\mathrm{dif}$', fixed=True, fixed_norm=14.67)
n.add_poiss_model('iso', '$A_\mathrm{iso}$', [0,2], False)
n.add_poiss_model('gce', '$A_\mathrm{gce}$', [0,2], False)
n.add_poiss_model('bub', '$A_\mathrm{bub}$', [0,2], False)

n.add_non_poiss_model('gce',
                      ['$A_\mathrm{gce}^\mathrm{ps}$','$n_1^\mathrm{gce}$','$n_2^\mathrm{gce}$','$S_b^{(1), \mathrm{gce}}$'],
                      [[-6,1],[2.05,30],[-2,1.95],[0.05,40]],
                      [True,False,False,False])
n.add_non_poiss_model('dsk',
                      ['$A_\mathrm{dsk}^\mathrm{ps}$','$n_1^\mathrm{dsk}$','$n_2^\mathrm{dsk}$','$S_b^{(1), \mathrm{dsk}}$'],
                      [[-6,1],[2.05,30],[-2,1.95],[0.05,40]],
                      [True,False,False,False])

pc_inst = pc.PSFCorrection(psf_sigma_deg=0.1812)
f_ary, df_rho_div_f_ary = pc_inst.f_ary, pc_inst.df_rho_div_f_ary

n.configure_for_scan(f_ary, df_rho_div_f_ary, nexp=1)

n.perform_scan(nlive=100)
