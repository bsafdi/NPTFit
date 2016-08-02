
# coding: utf-8

# In[1]:

import numpy as np
import scipy
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

import healpy as hp
import mpmath as mp
from astropy.io import fits

# import pulsars.special as spc
# import pulsars.masks as masks
# import pulsars.CTB as CTB
# import pulsars.psf as psf
# import pulsars.diffuse_fermi_not_pdep as df

# import globalVariables as gv
# import make_maps as mm
# import compute_PSF as CPSF

import logging
from contextlib import contextmanager
import sys, os

import config.smooth_with_psf as swp
import config.smooth_with_king as swk

from tempfile import mkdtemp
cachedir = mkdtemp()
from joblib import Memory
memory = Memory(cachedir=cachedir, verbose=0)

size = np.size


class ps_template:
    # Class to create a point source template at a given location
    # Updated to include up and down binning in nside
    # Nick Rodd 14-10-15
    def __init__(self,ell,b,nside,psf):
        self.nside = nside
        self.npix = hp.nside2npix(self.nside)
        self.ell = ell
        self.b = b
        self.psf = psf

        self.ell_and_b_to_pix()
        self.make_ps_template()

    def ell_and_b_to_pix(self):
        # Make a point source at the origin and rotate it later (get better behaviour near poles)
        self.phi = 0
        self.theta = np.pi/2
        # Upbin the nside depending on the value of the psf
        # Idea is to make the sure the pixel size is always smaller than 1/5 of the PSF
        # Bin to a minimum of 256, and a maximum of 1028 (for reasons of speed)
        if (self.psf < 0.15):
            nsideupbin = 512
            if (self.psf < 0.075):
                nsideupbin=1024
        else:
            nsideupbin=256
        self.nsideupbin = nsideupbin
        self.npixupbin = hp.nside2npix(self.nsideupbin)
        self.pix_num = hp.ang2pix(self.nsideupbin,self.theta,self.phi)

    def make_ps_template(self):
        the_map = np.zeros(self.npixupbin)
        the_map[self.pix_num] = 1
        # Smooth using a Gaussian
        self.sgp = swp.smooth_gaussian_psf(self.psf,the_map)
        # Now rotate the map
        thetarot, phirot = hp.pix2ang(self.nsideupbin, np.arange(self.npixupbin))
        r = hp.Rotator(rot = [self.ell,self.b], coord=None, inv=False, deg=False, eulertype='ZYX')
        thetarot, phirot = r(thetarot, phirot)
        pixrot = hp.ang2pix(self.nsideupbin, thetarot, phirot, nest=False)
        smooth_ps_map_rot = self.sgp.smooth_the_pixel(self.pix_num)[..., pixrot]
        # Downbin for the final map
        smooth_ps_map = hp.ud_grade(smooth_ps_map_rot,self.nside,power=-2)
        # Readjust to 1 as gets mucked up by rotation and downbinning sometimes
        self.smooth_ps_map = smooth_ps_map/np.sum(smooth_ps_map)

@memory.cache
def ps_template_king_fast(ell,b,nside,maps_dir,Ebin,eventclass,eventtype):
    # This function is identical to ps_template_king class, but makes use of memory class
    # which means it is much faster if calling multiple times
    # Currently I've kept both as this require joblib to be installed to run
    # Make it a function for easier use with memory
    # As opposed to the above this uses a King function and not a Gaussian
    # Updated to include up and down binning in nside
    # Nick Rodd 14-10-15
    phi = 0
    theta = np.pi/2

    # Upbin the nside depending on the value of the psf
    # Idea is to make the sure the pixel size is always smaller than 1/5 of the PSF
    # Bin to a minimum of 256, and a maximum of 1028 (for reasons of speed)
    if (Ebin >= 11):
        nsideupbin = 512
        if (Ebin >= 15):
            nsideupbin=1024
    else:
        nsideupbin=256
    nsideupbin = nsideupbin
    npixupbin = hp.nside2npix(nsideupbin)
    pix_num = hp.ang2pix(nsideupbin,theta,phi)

    the_map = np.zeros(npixupbin)
    the_map[pix_num] = 1
    # Smooth using a king function
    skf = swk.smooth_king_psf(maps_dir,the_map,Ebin,eventclass,eventtype)
    # Now rotate the map
    thetarot, phirot = hp.pix2ang(nsideupbin, np.arange(npixupbin))
    r = hp.Rotator(rot = [ell,b], coord=None, inv=False, deg=False, eulertype='ZYX')
    thetarot, phirot = r(thetarot, phirot)
    pixrot = hp.ang2pix(nsideupbin, thetarot, phirot, nest=False)
    smooth_ps_map_rot = skf.smooth_the_pixel(pix_num)[..., pixrot]
    # Downbin for the final map
    smooth_ps_map = hp.ud_grade(smooth_ps_map_rot,nside,power=-2)
    # Normalise back to 1 (gets mucked up by rotation and downbinning)
    return smooth_ps_map/np.sum(smooth_ps_map)


class ps_template_king:
    # Class to create a point source template at a given location
    # As opposed to the above this uses a King function and not a Gaussian
    # Updated to include up and down binning in nside
    # Nick Rodd 14-10-15
    def __init__(self,ell,b,nside,maps_dir,Ebin,eventclass,eventtype):
        self.nside = nside
        self.npix = hp.nside2npix(self.nside)
        self.ell = ell
        self.b = b
        self.maps_dir = maps_dir
        self.Ebin = Ebin
        self.eventclass = eventclass
        self.eventtype = eventtype

        self.ell_and_b_to_pix()
        self.make_ps_template()

    def ell_and_b_to_pix(self):
        # Make a point source at the origin and rotate it later (get better behaviour near poles)
        self.phi = 0
        self.theta = np.pi/2

        # Upbin the nside depending on the value of the psf
        # Idea is to make the sure the pixel size is always smaller than 1/5 of the PSF
        # Bin to a minimum of 256, and a maximum of 1028 (for reasons of speed)
        if (self.Ebin >= 11):
            nsideupbin = 512
            if (self.Ebin >= 15):
                nsideupbin=1024
        else:
            nsideupbin=256
        self.nsideupbin = nsideupbin
        self.npixupbin = hp.nside2npix(self.nsideupbin)
        self.pix_num = hp.ang2pix(self.nsideupbin,self.theta,self.phi)

    def make_ps_template(self):
        the_map = np.zeros(self.npixupbin)
        the_map[self.pix_num] = 1
        # Smooth using a king function
        self.skf = swk.smooth_king_psf(self.maps_dir,the_map,self.Ebin,self.eventclass,self.eventtype)
        # Now rotate the map
        thetarot, phirot = hp.pix2ang(self.nsideupbin, np.arange(self.npixupbin))
        r = hp.Rotator(rot = [self.ell,self.b], coord=None, inv=False, deg=False, eulertype='ZYX')
        thetarot, phirot = r(thetarot, phirot)
        pixrot = hp.ang2pix(self.nsideupbin, thetarot, phirot, nest=False)
        smooth_ps_map_rot = self.skf.smooth_the_pixel(self.pix_num)[..., pixrot]
        # Downbin for the final map
        smooth_ps_map = hp.ud_grade(smooth_ps_map_rot,self.nside,power=-2)
        # Normalise back to 1 (gets mucked up by rotation and downbinning)
        self.smooth_ps_map = smooth_ps_map/np.sum(smooth_ps_map)

def photons_per_ann(en_array):
    #returns dN/dE per annihilation at values of en_array, using Bergstrom fit to spectrum
    return 2./m_chi*0.42*np.exp(-8.*en_array/m_chi)/((en_array/m_chi)**1.5 + .00014)

def return_diff_p8_um(fits_file ,CTB_en_bins,CTB_en_min,CTB_en_max,NSIDE=128):
    # Designed to replace function below just for p8 diffuse models
    NPIX = hp.nside2npix(NSIDE)
    diff_fits = fits.open(fits_file)
    diff_master_p8 = diff_fits[0].data
    difftoreturn = [hp.ud_grade(diff_master_p8[:,i],NSIDE)*(CTB_en_bins[i+1-CTB_en_min]-CTB_en_bins[i-CTB_en_min])*4*np.pi/NPIX*10**(3) for i in range(CTB_en_min,CTB_en_max)]
    return difftoreturn

def return_diff_um(fits_file ,CTB_en_bins,CTB_en_min,CTB_en_max,NSIDE=128,mode='p6',is_p8 = False):
    NPIX = hp.nside2npix(NSIDE)
    
    diff_fits = fits.open(fits_file)
    if not is_p8:
                #load diffuse-background model
        diff_master = diff_fits[1].data['COARSEBINS']
        if mode=='p7':
            diffp7 = [hp.ud_grade(diff_master[0,1,::,i],NSIDE)*(CTB_en_bins[i+1-CTB_en_min]-CTB_en_bins[i-CTB_en_min])*4*np.pi/NPIX*10**(3) for i in range(CTB_en_min,CTB_en_max)]
            return diffp7
        else:
            diffp6 = [hp.ud_grade(diff_master[0,0,::,i],NSIDE)*(CTB_en_bins[i+1-CTB_en_min]-CTB_en_bins[i-CTB_en_min])*4*np.pi/NPIX*10**(3) for i in range(CTB_en_min,CTB_en_max)]
            return diffp6
    else:
        diff_master_p8 = diff_fits[0].data
        diffp6 = [hp.ud_grade(diff_master_p8[:,i],NSIDE)*(CTB_en_bins[i+1-CTB_en_min]-CTB_en_bins[i-CTB_en_min])*4*np.pi/NPIX*10**(3) for i in range(CTB_en_min,CTB_en_max)]
        return diffp6

def return_bubbles(bubble_map,CTB_en_min,CTB_en_max,CTB_exposure_maps,NSIDE=128):
    NPIX = hp.nside2npix(NSIDE)

    bubble_maps = [np.load(bubble_map)[0] for i in range(CTB_en_max -CTB_en_min)]#[CTB_en_min-1:CTB_en_max-1]
    bubble_fixed = np.array( [ hp.ud_grade(bubs,NSIDE) for bubs in bubble_maps])
    fermi_bubs_count_maps_CTB_binned = bubble_fixed* CTB_exposure_maps * 4*np.pi/NPIX

    return fermi_bubs_count_maps_CTB_binned

def return_isotropic(CTB_en_min,CTB_en_max,CTB_exposure_maps,NSIDE=128):
    NPIX = hp.nside2npix(NSIDE)
    iso_count_maps_CTB_binned = np.array(CTB_exposure_maps * 4*np.pi/NPIX,'double')

    iso_norm = [map/np.mean(map) for map in iso_count_maps_CTB_binned]

    return iso_norm

def return_template(map_name,CTB_en_min,CTB_en_max,CTB_exposure_maps,NSIDE = 128,already_has_exposure = False):
    #map name should be the file name of the intensity map (not exposure corrected)
    NPIX = hp.nside2npix(NSIDE)
 
    maps = np.load(map_name)[CTB_en_min-1:CTB_en_max-1]
    maps_fixed = np.array( [ hp.ud_grade(a_map,NSIDE,power=-2) for a_map in maps])
    if already_has_exposure:
        count_maps_CTB_binned = maps_fixed[CTB_en_min-1:CTB_en_max-1] 
    else:
        count_maps_CTB_binned = maps_fixed[CTB_en_min-1:CTB_en_max-1]* CTB_exposure_maps * 4*np.pi/NPIX

    return count_maps_CTB_binned


def make_xbg_PSF(A_vec, x_vec):#, A_NFW):
    '''Make masked, integrated, smoothed, exposure-corrected count maps for
        linear combination of background components.'''
    #x_diff_PSF, x_bubs_PSF, x_iso_PSF, X_NFW_PSF = x_PSF_um
    A_diff, A_bubs, A_iso = A_vec
    x_diff_PSF_um, x_bubs_PSF_um, x_iso_PSF_um = x_vec
    return A_diff*x_diff_PSF_um + A_bubs*x_bubs_PSF_um + A_iso*x_iso_PSF_um

def return_masked(the_map,the_mask):
    map_masked = hp.ma(the_map)
    map_masked.mask = the_mask
    return map_masked

def return_compressed(the_map,the_mask):
    return return_masked(the_map,the_mask).compressed()
