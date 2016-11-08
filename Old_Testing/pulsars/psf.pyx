import numpy as np
import healpy as hp
cimport cython
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

import logging
logger = logging.getLogger(__name__)

cdef int num_f_bins = 10
cdef int n_gauss = 50000
cdef int n_pts_per_gauss = 1000
cdef double f_trunc = 0.01 #10. / n_pts_per_gauss

cdef extern from "math.h":
    double log(double x) nogil
    double exp(double x) nogil
    double pow(double x, double y) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef np.ndarray make_PSF_int(str psf_dir, int nside, double sigma_PSF_deg,double sigma_PSF, np.ndarray theta_c, np.ndarray phi_c):
    cdef int gauss
    cdef np.ndarray dtheta = np.zeros(n_pts_per_gauss,dtype=DTYPE) 
    cdef np.ndarray dphi = np.zeros(n_pts_per_gauss,dtype=DTYPE)
    cdef np.ndarray theta_temp = np.zeros(n_pts_per_gauss,dtype=DTYPE)
    cdef np.ndarray theta_temp_formap = np.zeros(n_pts_per_gauss,dtype=DTYPE)
    cdef np.ndarray phi_temp = np.zeros(n_pts_per_gauss,dtype=DTYPE)
    cdef np.ndarray theta = np.zeros(n_pts_per_gauss,dtype=DTYPE) 
    cdef np.ndarray phi = np.zeros(n_pts_per_gauss,dtype=DTYPE) 
    cdef np.ndarray pixel_hist #= np.zeros(mx-mn,dtype=DTYPE) 
    cdef double mn,mx
    cdef int[::1] pixel
    cdef list outlist=[] #=np.array((n_gauss,mx-mn),dtype=DTYPE)
    for gauss in range(n_gauss):
        dtheta = np.random.normal(scale=sigma_PSF, size=n_pts_per_gauss)
        dphi = np.random.normal(scale=sigma_PSF, size=n_pts_per_gauss)/(np.sin(theta_c[gauss]+dtheta/2))
        theta_temp_formap = theta_c[gauss] + dtheta
        theta_temp = theta_c[gauss] + dtheta
        phi_temp = phi_c[gauss] + dphi

        theta_temp[np.where(theta_temp_formap > np.pi)[0]] = 2*np.pi-theta_temp[np.where(theta_temp_formap > np.pi)[0]]
        theta_temp[np.where(theta_temp_formap < 0)[0]] = -theta_temp[np.where(theta_temp_formap < 0)[0]]

        phi_temp[np.where(theta_temp_formap < 0)[0]] += np.pi
        phi_temp[np.where(theta_temp_formap > np.pi)[0]] += np.pi

        phi_temp = np.mod(phi_temp, 2*np.pi)

        theta = theta_temp[np.where((theta_temp <= np.pi) & (theta_temp >= 0))[0]]
        phi = phi_temp[np.where((theta_temp <= np.pi) & (theta_temp >= 0))[0]]

        pixel = np.array(hp.ang2pix(nside, theta, phi),dtype='int32')

        mn = np.min(pixel)
        mx = np.max(pixel) + 1
        pixel_hist = np.histogram(pixel, bins=mx-mn, range=(mn, mx), normed=1)[0]
        outlist.append(pixel_hist)
        # f_values_temp = np.append(f_values_temp, pixel_hist)
        # if (gauss + 1) % 1000 == 0:
        #     f_values = np.append(f_values, np.array(f_values_temp).ravel())
        #     f_values_temp = np.array([])
    cdef np.ndarray f_values = np.concatenate(outlist) 
    return f_values


def make_PSF(psf_dir, nside, sigma_PSF_deg):
    sigma_PSF = sigma_PSF_deg*np.pi/180.
    f_values = np.array([])
    f_values_temp = np.array([])
    
    xyz = np.random.normal(size=(n_gauss, 3))
    xyz_unit = np.divide(xyz, np.linalg.norm(xyz, axis=1)[:, None])
    theta_c = np.arccos(xyz_unit[:,2])
    phi_c = np.arctan2(xyz_unit[:,1], xyz_unit[:,0])
        
    f_values=make_PSF_int(psf_dir,nside,sigma_PSF_deg,sigma_PSF,theta_c,phi_c)

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
