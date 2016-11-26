###############################################################################
# create_mask.py
###############################################################################
#
# Creates a Boolean mask where pixels labelled as true are masked and those
# labelled false are unmasked.
#
# Note throughout we adjust from b to theta = 90-b, as this is what healpy
# uses.
#
# Note also all inputs are in degrees.
#
# NB: this module fundamentally assumes that the analysis is being performed
# on a spherical region pixelised using HEALPix. If this is not the case, the
# mask must be computed differently.
#
###############################################################################

import numpy as np
import healpy as hp


def make_plane_mask(band_mask_range, nside):
    """ Masks within |b| < band_mask_range
    """
    mask_none = np.arange(hp.nside2npix(nside))
    return (np.radians(90-band_mask_range) < hp.pix2ang(nside, mask_none)[0]) * \
           (hp.pix2ang(nside, mask_none)[0] < np.radians(90+band_mask_range))


def make_long_mask(l_deg_min, l_deg_max, nside):
    """ Masks outside l_deg_min < l < l_deg_max
    """
    mask_none = np.arange(hp.nside2npix(nside))
    return (np.radians(l_deg_max) < hp.pix2ang(nside, mask_none)[1]) * \
           (hp.pix2ang(nside, mask_none)[1] < np.radians(360 + l_deg_min))


def make_lat_mask(b_deg_min, b_deg_max, nside):
    """ Masks outside b_deg_min < b < b_deg_max
    """
    mask_none = np.arange(hp.nside2npix(nside))
    return np.logical_not(
           (np.radians(90-b_deg_max) < hp.pix2ang(nside, mask_none)[0]) *
           (hp.pix2ang(nside, mask_none)[0] < np.radians(90-b_deg_min)))


def make_ring_mask(inner, outer, ring_b, ring_l, nside):
    """ Masks outside inner < r < outer, of a ring centred at (ring_b,ring_l)
    """
    mask_none = np.arange(hp.nside2npix(nside))
    return np.logical_not(
           (np.cos(np.radians(inner)) >=
            np.dot(hp.ang2vec(np.radians(90-ring_b),
                   np.radians(ring_l)), hp.pix2vec(nside, mask_none))) *
           (np.dot(hp.ang2vec(np.radians(90-ring_b),
            np.radians(ring_l)), hp.pix2vec(nside, mask_none)) >=
            np.cos(np.radians(outer))))


def make_mask_total(nside=128,
                    band_mask=False, band_mask_range=30,
                    l_mask=False, l_deg_min=-30, l_deg_max=30,
                    b_mask=False, b_deg_min=-30, b_deg_max=30,
                    mask_ring=False, inner=0, outer=30,
                    ring_b=0, ring_l=0,
                    custom_mask=None):
    """ Combines band, l, b, ring, and custom masks into a single mask
    """

    # Initialise an array where no pixels are masked
    mask_array = np.zeros(nside**2*12, dtype=bool)

    # Add masks depending on input
    if band_mask:
        mask_array += make_plane_mask(band_mask_range, nside)

    if l_mask:
        mask_array += make_long_mask(l_deg_min, l_deg_max, nside)

    if b_mask:
        mask_array += make_lat_mask(b_deg_min, b_deg_max, nside)

    if mask_ring:
        mask_array += make_ring_mask(inner, outer, ring_b, ring_l, nside)

    if custom_mask is not None:
        mask_array += custom_mask

    return mask_array
