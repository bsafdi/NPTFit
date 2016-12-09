# Mask-making tests

import sys
sys.path.append("../NPTFit")

import numpy as np

from NPTFit import create_mask as cm # Module for creating masks

def test_masks():

    example1 = cm.make_mask_total()
    example2 = cm.make_mask_total(band_mask = True, band_mask_range = 30)

    example3a = cm.make_mask_total(l_mask = False, l_deg_min = -30, l_deg_max = 30, 
                                   b_mask = True, b_deg_min = -30, b_deg_max = 30)
    example3b = cm.make_mask_total(l_mask = True, l_deg_min = -30, l_deg_max = 30, 
                                   b_mask = False, b_deg_min = -30, b_deg_max = 30)
    example3c = cm.make_mask_total(l_mask = True, l_deg_min = -30, l_deg_max = 30, 
                                  b_mask = True, b_deg_min = -30, b_deg_max = 30)

    example4a = cm.make_mask_total(mask_ring = True, inner = 0, outer = 30, ring_b = 0, ring_l = 0)
    example4b = cm.make_mask_total(mask_ring = True, inner = 30, outer = 180, ring_b = 0, ring_l = 0)
    example4c = cm.make_mask_total(mask_ring = True, inner = 30, outer = 90, ring_b = 0, ring_l = 0)
    example4d = cm.make_mask_total(mask_ring = True, inner = 0, outer = 30, ring_b = 45, ring_l = 45)

    pscmask=np.array(np.load('fermi_data/fermidata_pscmask.npy'), dtype=bool)
    example6 = cm.make_mask_total(band_mask = True, band_mask_range = 2,
                                  mask_ring = True, inner = 0, outer = 30,
                                  custom_mask = pscmask)
