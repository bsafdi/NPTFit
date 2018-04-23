###############################################################################
# psf_compute_cart.py
###############################################################################
#
# The statistics of non-poissonian templates is modified by the non-zero point
# spread functions associated with real instruments. Here we calculate that
# correction for an arbitrary user specified PSF.
#
# The computation is performed on a periodic cartesian grid, rather than a
# healpix map.
#
###############################################################################

from __future__ import absolute_import

import numpy as np
import healpy as hp
from . import pdf_sampler

def psf_corr(gridsize, pixarea, num_f_bins, n_psf, n_pts_per_psf, f_trunc, 
             psf_r_func, sample_psf_max, psf_samples):

    # Calculation is performed on a grid of size gridsize x gridsize

    # Setup pdf of the psf
    radial_pdf = lambda r: r * psf_r_func(r)
    rvals = np.linspace(0, sample_psf_max, psf_samples)
    pofr = radial_pdf(rvals)
    dist = pdf_sampler.PDFSampler(rvals, pofr)

    # Work on a periodic gridsize x gridsize Cartesian grid, of pixel size 
    # pixarea
    pixwidth = np.sqrt(pixarea)
    gridwidth = gridsize*pixwidth

    # Create an array of n_psf points to put down psfs
    x_c = np.random.uniform(0, gridwidth, n_psf)
    y_c = np.random.uniform(0, gridwidth, n_psf)

    # Put a point source down at each of these locations
    outlist = []
    for ps_i in range(n_psf):
        # For each point source put down n_pts_per_psf counts
        # Determine where they are placed on the map as determine by the psf
        dr = dist(n_pts_per_psf)
        dangle = np.random.uniform(0, 2 * np.pi, n_pts_per_psf)
        dx = dr * np.sin(dangle)
        dy = dr * np.cos(dangle)

        # Combine with position of point source to get the exact location
        x_arr = x_c[ps_i] + dx
        y_arr = y_c[ps_i] + dy

        # Convert these arrays into positions on the grid
        # Commented out code was old, where we gave the grid periodic boundary
        # conditions. Now work to account for the fact flux can get lost off the edge
        #x_loc = np.mod(np.floor(x_arr/pixwidth).astype(int),gridsize) 
        #y_loc = np.mod(np.floor(y_arr/pixwidth).astype(int),gridsize)
        x_loc = np.floor(x_arr/pixwidth).astype(int)
        y_loc = np.floor(y_arr/pixwidth).astype(int)

        # Keep only pixels within the grid
        keep = np.where((x_loc > -1) & (x_loc < gridsize) & 
                        (y_loc > -1) & (y_loc < gridsize))[0]

        x_keep = x_loc[keep]
        y_keep = y_loc[keep]

        # Move from a 2D array to a 1D array
        pixel = gridsize*y_keep + x_keep

        # From this information determine the flux fraction per pixel
        mn = np.min(pixel)
        mx = np.max(pixel) + 1
        pixel_hist = np.histogram(pixel, bins=mx - mn, range=(mn, mx))[0]

        # Normalize manually by the number of points generated, not the
        # number that survived
        pixel_hist = pixel_hist.astype(float)/float(n_pts_per_psf)
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
