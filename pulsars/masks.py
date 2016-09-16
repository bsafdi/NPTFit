# This code takes input desired masking regions and returns a boolean array of npix values, where pixels labelled as true are masked and those labelled as false are not 

import numpy as np
import healpy as hp

def mask_lat_band(lat_deg_min, lat_deg_max, nside):
    #make mask of band lat_deg_min < b < lat_deg_max
    mask_none = np.arange(hp.nside2npix(nside))
    return (np.radians(lat_deg_min) < hp.pix2ang(nside, mask_none)[0]) * \
            (hp.pix2ang(nside, mask_none)[0] < np.radians(lat_deg_max))

### Comment, Nick Rodd, MIT, 15 Sep 2015
# Note given the way l is defined by default (from 0 to 360, with the galactic centre at 0), it is much easier to define the not_long_band code first
# Might think a trick like l = ((l +180) mod 360)-180 would work to redefine l to be in [-180,180] with the galactic centre at 0, but it seems this gives slightly different answers to idl, so I've used the below method which agrees exactly 

def mask_not_long_band(long_deg_min, long_deg_max, nside):
    mask_none = np.arange(hp.nside2npix(nside))
    return (np.radians(long_deg_max) < hp.pix2ang(nside, mask_none)[1]) * \
           (hp.pix2ang(nside, mask_none)[1] < np.radians(360 + long_deg_min))

def mask_long_band(long_deg_min, long_deg_max, nside):
    #make mask of region outside band long_deg_min < l < long_deg_max
    return np.logical_not(mask_not_long_band(long_deg_min, long_deg_max, nside))

def mask_not_lat_band(lat_deg_min, lat_deg_max, nside):
    #make mask of region outside band lat_deg_min < b < lat_deg_max
    return np.logical_not(mask_lat_band(lat_deg_min, lat_deg_max, nside))

def mask_ring(ring_deg_min, ring_deg_max, center_theta_deg, center_phi_deg, nside):
    #make mask of region outside ring_deg_min < theta < ring_deg_max
    mask_none = np.arange(hp.nside2npix(nside))
    return (np.cos(np.radians(ring_deg_min)) >= np.dot(hp.ang2vec(np.radians(center_theta_deg), np.radians(center_phi_deg)), hp.pix2vec(nside, mask_none))) * \
            (np.dot(hp.ang2vec(np.radians(center_theta_deg), np.radians(center_phi_deg)), hp.pix2vec(nside, mask_none)) >= np.cos(np.radians(ring_deg_max)))
