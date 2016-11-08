import numpy as np
import healpy as hp

import logging
logger = logging.getLogger(__name__)

num_f_bins = 10
n_gauss = 50000
n_pts_per_gauss = 1000
f_trunc = 0.01 #10. / n_pts_per_gauss

#cdef double make_PSF_int(str psf_dir, int nside, double sigma_PSF_deg):


def make_PSF(psf_dir, nside, sigma_PSF_deg):
    sigma_PSF = sigma_PSF_deg*np.pi/180.
    f_values = np.array([])
    f_values_temp = np.array([])
    
    xyz = np.random.normal(size=(n_gauss, 3))
    xyz_unit = np.divide(xyz, np.linalg.norm(xyz, axis=1)[:, None])
    theta_c = np.arccos(xyz_unit[:,2])
    phi_c = np.arctan2(xyz_unit[:,1], xyz_unit[:,0])
        
    for gauss in range(n_gauss):
        dtheta = np.random.normal(scale=sigma_PSF, size=n_pts_per_gauss)
        dphi = np.random.normal(scale=sigma_PSF, size=n_pts_per_gauss)
        theta = np.mod(theta_c[gauss] + dtheta, np.pi)
        phi = np.mod(phi_c[gauss] + dphi, 2*np.pi)
        pixel = hp.ang2pix(nside, theta, phi)

        mn = np.min(pixel)
        mx = np.max(pixel) + 1
        pixel_hist = np.histogram(pixel, bins=mx-mn, range=(mn, mx), normed=1)[0]
        f_values_temp = np.append(f_values_temp, pixel_hist)
        if (gauss + 1) % 1000 == 0:
            f_values = np.append(f_values, np.array(f_values_temp).ravel())
            f_values_temp = np.array([])
            logging.info('%s', 'Gaussian ' + str(gauss+1) + ' done.')

    f_values_trunc = f_values[f_values >= f_trunc]
    rho_ary, f_bin_edges = np.histogram(f_values_trunc, bins=num_f_bins, range=(0.,1.))
    df = f_bin_edges[1] - f_bin_edges[0]
    f_ary = (f_bin_edges[:-1] + f_bin_edges[1:])/2.
    rho_ary = rho_ary / (df * n_gauss)
    f_ary_edge = f_bin_edges[:-1]        
    np.savetxt(psf_dir + 'f_ary-' + str(nside) + '-' + str(np.round(sigma_PSF_deg,3)) + '-' + str(num_f_bins) + '.dat',
               f_ary, fmt='%1.5f')
    np.savetxt(psf_dir + 'rho_ary-' + str(nside) + '-' + str(np.round(sigma_PSF_deg,3)) + '-' + str(num_f_bins) + '.dat',
               rho_ary, fmt='%1.5f')
    rho_ary = rho_ary / np.sum(df*f_ary*rho_ary)
    df_rho_div_f_ary = df*rho_ary / f_ary
    return f_ary, df_rho_div_f_ary
               
def load_PSF(psf_dir, nside, sigma_PSF_deg):
    f_ary = np.loadtxt(psf_dir + 'f_ary-' + str(nside) + '-' + str(np.round(sigma_PSF_deg,3)) + '-' + str(num_f_bins) + '.dat')
    df = f_ary[1] - f_ary[0]
    f_ary_edge = f_ary - df/2 + 0.0001
    rho_ary = np.loadtxt(psf_dir + 'rho_ary-' + str(nside) + '-' + str(np.round(sigma_PSF_deg,3)) + '-' + str(num_f_bins) + '.dat')
    
    rho_ary = rho_ary / np.sum(df*f_ary*rho_ary)
    df_rho_div_f_ary = df*rho_ary / f_ary
    return f_ary, df_rho_div_f_ary
