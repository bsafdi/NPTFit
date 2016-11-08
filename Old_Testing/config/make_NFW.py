import sys, os

# current_dir = os.getcwd()
# change_path = ".."
# os.chdir(change_path)

import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt

import healpy as hp
import mpmath as mp
from astropy.io import fits

import logging
from contextlib import contextmanager
import sys, os

from scipy import integrate, interpolate


size = np.size

class make_NFW:
    def __init__(self,en_array = np.array([3.0]),r_0 = 8.5,r_s = 20., gamma_nfw = 1.25, sigma_v = 1.7, rho_0 = 0.3, m_chi = 35.,deg_steps = 0.5/4, GC_lng = 0,nside = 128,data_maps_dir = 'False'):
        #variables for creating NFW template
        self.r_0 = r_0
        self.r_s             = r_s  #NFW scale radius in kpc
        self.gamma_nfw       = gamma_nfw #1.25  #gNFW power, use values from Daylan et al. fits
        self.sigma_v         = sigma_v  #in 10^-26 cm^3 / s
        self.rho_0           = rho_0   #in GeV / cm^3
        self.m_chi           = m_chi   #in GeV
        self.GC_lng = GC_lng #longitude for NFW center (Degrees)
        self.en_array = en_array

        #variables for making intesnsity map
        self.nside = nside
        self.npix = hp.nside2npix(self.nside)
        self.deg_steps = deg_steps

        #variables for saving the map
        self.data_maps_dir = data_maps_dir

        self.make_NFW_intensity_map()


    def photons_per_ann(self,en_array):
        #returns dN/dE per annihilation at values of en_array, using Bergstrom fit to spectrum
        return 2./self.m_chi*0.42*np.exp(-8.*en_array/self.m_chi)/((en_array/self.m_chi)**1.5 + .00014)

    def rho_dimless(self,r):
        #r in kpc
        rho_0_factor = np.power(self.r_0/self.r_s, self.gamma_nfw) * np.power(1. + self.r_0/self.r_s, 3. - self.gamma_nfw)
        return np.power(r/self.r_s, -self.gamma_nfw) * np.power(1. + r/self.r_s, self.gamma_nfw - 3.) * rho_0_factor

    def r_NFW(self,l, psi_deg):
        return np.sqrt(self.r_0**2. + l**2. - 2.*self.r_0*l*np.cos(np.radians(psi_deg)))

    def L_NFW_integral(self,psi_deg):
        return integrate.quad(lambda l: self.rho_dimless(self.r_NFW(l, psi_deg))**2., 0., 100.*self.r_s)[0]

    def make_NFW_intensity_map(self):
        self.map_string = 'NFW_map-En'+str(self.en_array[0])+ '-gamma-' + str(self.gamma_nfw)+ '-nside-' + str(self.nside) + '-GC_lng-'+str(self.GC_lng)

        if self.data_maps_dir == 'False' or not os.path.isfile(self.data_maps_dir + self.map_string + '.npy') :
            print 'Need to make NFW intensity map ...'
            self.do_NFW_intensity_map()
            if self.data_maps_dir != 'False':
                self.save_NFW_intensity_map()
        else:
            print 'Just loading NFW intensity map ...'
            self.load_NFW_intensity_map()

    def do_NFW_intensity_map(self):
        en_array = self.en_array
        psi_deg = np.arange(0., 180.5, self.deg_steps) #BS: put in /4
        intensity_NFW = self.sigma_v * self.rho_0**2. / (8. * np.pi * self.m_chi**2.) * 3.08E-5 * np.vectorize(self.L_NFW_integral)(psi_deg)
        intensity_NFW_interp = interpolate.interp1d(psi_deg, intensity_NFW)
        GC_vec = [np.cos(np.deg2rad(self.GC_lng)), np.sin(np.deg2rad(self.GC_lng)), 0.]
        psi_deg_pixels = np.array([np.degrees(np.arccos(np.dot(GC_vec, hp.pix2vec(self.nside, pix)))) for pix in range(self.npix)])

        self.NFW_profile = np.outer(self.photons_per_ann(en_array), intensity_NFW_interp(psi_deg_pixels))

    def save_NFW_intensity_map(self):
        np.save(self.data_maps_dir + self.map_string, self.NFW_profile)

    def load_NFW_intensity_map(self):
        self.NFW_profile = np.load(self.data_maps_dir + self.map_string + '.npy')
