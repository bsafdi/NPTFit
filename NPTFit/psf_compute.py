###############################################################################
# psf_compute.py
###############################################################################
#
# The statistics of non-poissonian templates is modified by the non-zero point
# spread functions associated with real instruments. Here we calculate that
# correction for an arbitrary user specified PSF.
#
###############################################################################

from __future__ import absolute_import

import numpy as np
import healpy as hp
from . import pdf_sampler


def psf_corr(nside, num_f_bins, n_psf, n_pts_per_psf, f_trunc, psf_r_func,
             sample_psf_max, psf_samples):
    # Setup pdf of the psf
    radial_pdf = lambda r: r * psf_r_func(r)
    rvals = np.linspace(0, sample_psf_max, psf_samples)
    pofr = radial_pdf(rvals)
    dist = pdf_sampler.PDFSampler(rvals, pofr)

    # Create an array of n_psf points to put down psfs
    # Establish an array of n_ps unit vectors
    # By sampling vals from a normal, end up with uniform normed vectors
    xyz = np.random.normal(size=(n_psf, 3))
    xyz_unit = np.divide(xyz, np.linalg.norm(xyz, axis=1)[:, None])

    # Convert to array of theta and phi values
    # theta = arccos(z/r), and here r=1. Similar expression for phi
    theta_c = np.arccos(xyz_unit[:, 2])
    phi_c = np.arctan2(xyz_unit[:, 1], xyz_unit[:, 0])

    # Now put a point source down at each of these locations
    outlist = []
    for ps_i in range(n_psf):
        # For each point source put down n_pts_per_psf counts
        # Determine where they are placed on the map as determine by the psf
        dr = dist(n_pts_per_psf)
        dangle = np.random.uniform(0, 2 * np.pi, n_pts_per_psf)
        dtheta = dr * np.sin(dangle)
        dphi = dr * np.cos(dangle) / (np.sin(theta_c[ps_i] + dtheta / 2))

        # Now combine with position of point source to get the exact location
        theta_base = theta_c[ps_i] + dtheta
        phi_base = phi_c[ps_i] + dphi

        # Want 0 <= theta < pi; 0 <= phi < 2pi
        # Carefully map to ensure this is true
        theta_remap_north = np.where(theta_base > np.pi)[0]
        theta_base[theta_remap_north] = 2 * np.pi - theta_base[theta_remap_north]
        theta_remap_south = np.where(theta_base < 0)[0]
        theta_base[theta_remap_south] = -theta_base[theta_remap_south]

        phi_base[theta_remap_north] += np.pi
        phi_base[theta_remap_south] += np.pi
        phi_base = np.mod(phi_base, 2 * np.pi)

        # As the PSF extends to infinity, if draw a value a long way from the
        # centre can occasionally still have a theta value outside the default
        # range above. For any sensible PSF (much smaller than the size of the
        # sky) this happens rarely. As such we just cut these values out.
        good_val = np.where((theta_base <= np.pi) & (theta_base >= 0))[0]
        theta = theta_base[good_val]
        phi = phi_base[good_val]

        # Convert these values back to a healpix pixel
        pixel = hp.ang2pix(nside, theta, phi)

        # From this information determine the flux fraction per pixel
        mn = np.min(pixel)
        mx = np.max(pixel) + 1
        pixel_hist = np.histogram(pixel, bins=mx - mn, range=(mn, mx), normed=1)[
            0]
        outlist.append(pixel_hist)

    f_values = np.concatenate(outlist)
    # f_values is now the full list of flux fractions from all psfs
    # Ignore values which fall below the cutoff f_trunc
    f_values_trunc = f_values[f_values >= f_trunc]

    # Rebin into the user defined number of bins
    rho_ary, f_bin_edges = np.histogram(f_values_trunc, bins=num_f_bins,
                                        range=(0., 1.))

    # Convert to output format
    df = f_bin_edges[1] - f_bin_edges[0]
    f_ary = (f_bin_edges[:-1] + f_bin_edges[1:]) / 2.
    rho_ary = rho_ary / (df * n_psf)
    rho_ary /= np.sum(df * f_ary * rho_ary)
    df_rho_div_f_ary = df * rho_ary / f_ary

    return f_ary, df_rho_div_f_ary
