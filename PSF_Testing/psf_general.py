#
# NAME:
#  psf_king.py
#
# PURPOSE:
#  Define a class that rapidly calculates rho(f), a function that accounts for
#  the effect of the PSF in an NPTF run.
#  For details of rho(f) and its use see Sec 2.2 of 1104.0010
#  For details of the Fermi PSF, see: http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_LAT_IRFs/IRF_PSF.html 
#
# INPUT:
#  Parameters to define a king function (fcore,score,gcore,stail,gtail,SpE)
#  These parameters can be an array if sampling from multiple king functions is required
#  If so can also frac array which determines what fraction of the sampling points come
#  from each king function
#
# OUTPUT:
#  rho(f) and f values
#
# REVISION HISTORY:
#  2015-Sep-30 Created by modifying psf.py by Nick Rodd, MIT

import numpy as np
import healpy as hp
import gauss_pdf as gp

class psf_general:
    """ A class to rapidly evaluate and save, or load, the Fermi PSF in a form needed for NPTF scans"""
    def __init__(self,sigma_PSF_deg,psf_dir='/zfs/nrodd/NPTFit/PSF_Testing/psf_dir/',nside=128,num_f_bins=10,n_psf=50000,n_pts_per_psf=1000,f_trunc=0.01):

        self.sigma_PSF_deg = sigma_PSF_deg
        self.psf_dir = psf_dir
        self.nside = nside
        self.num_f_bins = num_f_bins
        self.n_psf = n_psf
        self.n_pts_per_psf = n_pts_per_psf
        self.f_trunc = f_trunc

        # Load instance of PSF
        self.gpi = gp.gauss_pdf(sigma_PSF_deg)

    def make_PSF(self):
        f_values = np.array([])
        f_values_temp = np.array([])

        # First establish an array of n_ps unit vectors
        # By sampling vals from a normal, end up with uniform normed vectors
        xyz = np.random.normal(size=(self.n_psf, 3))
        xyz_unit = np.divide(xyz, np.linalg.norm(xyz, axis=1)[:, None])

        # Convert to array of theta and phi values
        # Recall theta = arccos(z/r), and here r=1. Similar expression for phi
        theta_c = np.arccos(xyz_unit[:,2])
        phi_c = np.arctan2(xyz_unit[:,1], xyz_unit[:,0])

        # Now loop through each point source
        pixel_hist=np.array([])
        outlist=[]
        for i in range(self.n_psf):

            # First load array of delta p vals then convert to dtheta and dphi
            dp = self.gpi.gauss_pdf_draw(self.n_pts_per_psf)
            dangle = np.random.uniform(0,2*np.pi,self.n_pts_per_psf)
            dtheta = dp*np.sin(dangle)
            dphi = dp*np.cos(dangle)/(np.sin(theta_c[i]+dtheta/2))

            # Now combine with position of point source
            theta_temp = theta_c[i] + dtheta
            phi_temp = phi_c[i] + dphi

            # Do the theta mapping carefully
            theta_change_where_up = np.where(theta_temp > np.pi)[0]
            theta_temp[theta_change_where_up] = 2*np.pi-theta_temp[theta_change_where_up]
            theta_change_where_down = np.where(theta_temp < 0)[0]
            theta_temp[theta_change_where_down] = -theta_temp[theta_change_where_down]

            phi_temp[theta_change_where_down] += np.pi
            phi_temp[theta_change_where_up] += np.pi

            phi_temp = np.mod(phi_temp, 2*np.pi)

            # Now as the kings function extends to infty, chance the above still won't
            # get everything within the right theta range. As this happens very rarely
            # better to just cut such pixels out rather than slow code down to deal with
            good_val = np.where((theta_temp <= np.pi) & (theta_temp >= 0))[0]
            theta = theta_temp[good_val]
            phi = phi_temp[good_val]

            # Now map to a healpix pixel
            pixel = hp.ang2pix(self.nside, theta, phi)

            # Now want to determine the flux fraction per pixel to get rho(f)
            mn = np.min(pixel)
            mx = np.max(pixel) + 1
            pixel_hist = np.histogram(pixel, bins=mx-mn, range=(mn, mx), normed=1)[0]
            outlist.append(pixel_hist)

        f_values = np.concatenate(outlist)

        f_values_trunc = f_values[f_values >= self.f_trunc]
        rho_ary, f_bin_edges = np.histogram(f_values_trunc, bins=self.num_f_bins, range=(0.,1.))
        df = f_bin_edges[1] - f_bin_edges[0]
        f_ary = (f_bin_edges[:-1] + f_bin_edges[1:])/2.
        rho_ary = rho_ary / (df * self.n_psf)
        f_ary_edge = f_bin_edges[:-1]
        # NB: can't save all params as string would be too long so use save tag
        #np.savetxt(self.psf_dir+'f_ary-' + self.save_tag + '-'+ str(self.nside) + '-' + str(self.num_f_bins) + '.dat',f_ary)
        #np.savetxt(self.psf_dir+'rho_ary-' + self.save_tag + '-'+ str(self.nside) + '-' + str(self.num_f_bins) + '.dat',rho_ary)
        rho_ary = rho_ary / np.sum(df*f_ary*rho_ary)
        df_rho_div_f_ary = df*rho_ary / f_ary

        return f_ary, df_rho_div_f_ary

    def load_PSF(self):
        f_ary = np.loadtxt(self.psf_dir+'f_ary-' + self.save_tag + '-'+ str(self.nside) + '-' + str(self.num_f_bins) + '.dat')
        df = f_ary[1] - f_ary[0]
        f_ary_edge = f_ary - df/2 + 0.0001
        rho_ary = np.loadtxt(self.psf_dir+'rho_ary-' + self.save_tag + '-'+ str(self.nside) + '-' + str(self.num_f_bins) + '.dat')

        rho_ary = rho_ary / np.sum(df*f_ary*rho_ary)
        df_rho_div_f_ary = df*rho_ary / f_ary

        return f_ary, df_rho_div_f_ary 
