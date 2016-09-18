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

import sys, os
import numpy as np
import healpy as hp
import pulsars.king_pdf as kpi

class psf_king:
    """ A class to rapidly evaluate and save, or load, the Fermi PSF in a form needed for NPTF scans"""
    def __init__(self,fcore,score,gcore,stail,gtail,SpE,frac=[1],psf_dir='',nside=128,
                 num_f_bins=10,n_ps=50000,n_pts_per_king=1000,f_trunc=0.01,
                 king_sampling=10000,save_tag=''):

        # Convert floats to arrays
        fcore = np.array([fcore]).flatten()
        score = np.array([score]).flatten()
        gcore = np.array([gcore]).flatten()
        stail = np.array([stail]).flatten()
        gtail = np.array([gtail]).flatten()
        SpE   = np.array([SpE]).flatten()
        frac  = np.array([frac]).flatten()

        # Store basic variables to self
        self.fcore = fcore
        self.score = score
        self.gcore = gcore
        self.stail = stail
        self.gtail = gtail
        self.SpE = SpE
        self.psf_dir = psf_dir
        self.nside = nside
        self.num_f_bins = num_f_bins
        self.n_ps = n_ps
        self.f_trunc = f_trunc
        self.save_tag = save_tag

        # Setup equal frac if not user defined, and then distribute n_pts_per_king using frac
        n_king = len(self.fcore)
        self.n_king = n_king
        if not len(self.fcore) == len(frac):
            frac = np.empty(n_king)
            frac.fill(1./n_king)

        n_pts_per_king_array = np.empty(n_king)
        for i in range(n_king):
            n_pts_per_king_array[i]=int(round(n_pts_per_king*frac[i]))
        self.n_pts_per_king_array = n_pts_per_king_array
        self.total_n_pts_per_king = int(np.sum(n_pts_per_king_array))

        # Import king function pdfs
        kp = np.array([])
        for i in range(n_king):
            kp = np.append(kp,kpi.king_pdf(fcore[i],score[i],gcore[i],stail[i],gtail[i],SpE[i],king_sampling=king_sampling))
        self.kp = kp

    def make_PSF(self):
        f_values = np.array([])
        f_values_temp = np.array([])

        # First establish an array of n_ps unit vectors
        # By sampling vals from a normal, end up with uniform normed vectors
        xyz = np.random.normal(size=(self.n_ps, 3))
        xyz_unit = np.divide(xyz, np.linalg.norm(xyz, axis=1)[:, None])

        # Convert to array of theta and phi values
        # Recall theta = arccos(z/r), and here r=1. Similar expression for phi
        theta_c = np.arccos(xyz_unit[:,2])
        phi_c = np.arctan2(xyz_unit[:,1], xyz_unit[:,0])

        # Now loop through each point source
        pixel_hist=np.array([])
        outlist=[]
        for i in range(self.n_ps):
            # First load array of delta p vals then convert to dtheta and dphi
            dp = np.array([])
            for j in range(self.n_king):
                dp = np.append(dp,self.kp[j].king_pdf(self.n_pts_per_king_array[j]))
            dangle = np.random.uniform(0,2*np.pi,self.total_n_pts_per_king)
            dtheta = dp*np.sin(dangle)
            dphi = dp*np.cos(dangle)/(np.sin(theta_c[i]+dtheta/2))
            #dphi = dp*np.cos(dangle)

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

            # f_values_temp = np.append(f_values_temp, pixel_hist)
            # if (i + 1) % np.floor(self.n_ps/50.) == 0:
            #     f_values = np.append(f_values, np.array(f_values_temp).ravel())
            #     f_values_temp = np.array([])

        f_values = np.concatenate(outlist)

        f_values_trunc = f_values[f_values >= self.f_trunc]
        rho_ary, f_bin_edges = np.histogram(f_values_trunc, bins=self.num_f_bins, range=(0.,1.))
        df = f_bin_edges[1] - f_bin_edges[0]
        f_ary = (f_bin_edges[:-1] + f_bin_edges[1:])/2.
        rho_ary = rho_ary / (df * self.n_ps)
        f_ary_edge = f_bin_edges[:-1]
        # NB: can't save all params as string would be too long so use save tag
        np.savetxt(self.psf_dir+'f_ary-' + self.save_tag + '-'+ str(self.nside) + '-' + str(self.num_f_bins) + '.dat',f_ary)
        np.savetxt(self.psf_dir+'rho_ary-' + self.save_tag + '-'+ str(self.nside) + '-' + str(self.num_f_bins) + '.dat',rho_ary)
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
