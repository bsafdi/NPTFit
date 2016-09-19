#------------------------------------------------------------------------------
# NPTF PSF correction calculation
#------------------------------------------------------------------------------
##
## The statistics of non-poissonian templates is modified by the non-zero
## point spread functions associated with real instruments. Here we calculate
## that correction for an arbitrary user specified PSF.
##
#------------------------------------------------------------------------------


import numpy as np
import healpy as hp
import pdf_sampler

def psf_corr(num_f_bins, n_psf, n_pts_per_psf, f_trunc, psf_r_func, 
             sample_psf_max, psf_sample):

    # Setup pdf of the psf
    radial_pdf = lambda r: r*psf_r_func(r)
    rvals = np.linspace(0, sample_psf_max, psf_sample)
    pofr = radial_pdf(rvals)
    dist = pdf_sampler.pdf_sampler(rvals, pofr)    

    # First establish an array of n_ps unit vectors
    # By sampling vals from a normal, end up with uniform normed vectors
    xyz = np.random.normal(size=(n_psf, 3))
    xyz_unit = np.divide(xyz, np.linalg.norm(xyz, axis=1)[:, None])

    # Convert to array of theta and phi values
    # theta = arccos(z/r), and here r=1. Similar expression for phi
    theta_c = np.arccos(xyz_unit[:,2])
    phi_c = np.arctan2(xyz_unit[:,1], xyz_unit[:,0])

    # Loop through each point source
    for ps_i in range(n_psf):
        dp = dist(n_pts_per_psf)
        dangle = np.random.uniform(0,2*np.pi,n_pts_per_psf)
        dtheta = dp*np.sin(dangle)
        dphi = dp*np.cos(dangle)/(np.sin(theta_c[ps_i]+dtheta/2))

        # Now combine with position of point source
        theta_base = theta_c[ps_i] + dtheta
        phi_base = phi_c[ps_i] + dphi

        # Need 
        
        # Do the theta mapping carefully
        theta_change_where_up = np.where(theta_base > np.pi)[0]
        theta_base[theta_change_where_up] = 2*np.pi-theta_base[theta_change_where_up]
        theta_change_where_down = np.where(theta_base < 0)[0]
        theta_base[theta_change_where_down] = -theta_base[theta_change_where_down]

        phi_base[theta_change_where_down] += np.pi
        phi_base[theta_change_where_up] += np.pi

        phi_base = np.mod(phi_base, 2*np.pi)

        # Now as the kings function extends to infty, chance the above still won't
        # get everything within the right theta range. As this happens very rarely
        # better to just cut such pixels out rather than slow code down to deal with
        good_val = np.where((theta_base <= np.pi) & (theta_base >= 0))[0]
        theta = theta_base[good_val]
        phi = phi_base[good_val]

        # Now map to a healpix pixel
        pixel = hp.ang2pix(nside, theta, phi)

        # Now want to determine the flux fraction per pixel to get rho(f)
        mn = np.min(pixel)
        mx = np.max(pixel) + 1
        pixel_hist = np.histogram(pixel, bins=mx-mn, range=(mn, mx), normed=1)[0]
        outlist.append(pixel_hist)

    f_values = np.concatenate(outlist)

    f_values_trunc = f_values[f_values >= f_trunc]
    rho_ary, f_bin_edges = np.histogram(f_values_trunc, bins=num_f_bins, range=(0.,1.))
    df = f_bin_edges[1] - f_bin_edges[0]
    f_ary = (f_bin_edges[:-1] + f_bin_edges[1:])/2.
    rho_ary = rho_ary / (df * n_psf)
    f_ary_edge = f_bin_edges[:-1]
    rho_ary = rho_ary / np.sum(df*f_ary*rho_ary)
    df_rho_div_f_ary = df*rho_ary / f_ary
    return f_ary, df_rho_div_f_ary
