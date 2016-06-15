
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

# x_iso_PSF_um = np.sum(np.transpose(iso_array_norm* np.transpose([exp/np.mean(exp) for exp in CTB_exposure_maps])),axis=0)

# def main(CTB_en_bins, CTB_exposure_maps,bin_Number = 'All',normalization = 'Null'):


#     print 'The bin number is ', bin_Number

#     if normalization!='Null':
#         diff_array_norm, bubs_array_norm, iso_array_norm, NFW_array_norm = normalization
#         diff_array_norm, bubs_array_norm, iso_array_norm, NFW_array_norm = np.array(diff_array_norm), np.array(bubs_array_norm), np.array(iso_array_norm), np.array(NFW_array_norm)
#     else:
#         array_norm = np.array([1., 1.,1.,1.,1.,1.,1.,1.])
#         diff_array_norm, bubs_array_norm, iso_array_norm, NFW_array_norm = [array_norm,array_norm,array_norm,array_norm]
    
#     if bin_Number == 'All':
#         #sigma_PSF = gv.sigma_PSF
#         sigma_PSF = CPSF.main(gv.NSIDE,use_generic_sigma = False,just_sigma = True,CTB_en_bins = CTB_en_bins)*np.pi/180
#     else:
        
#         sigma_PSF = CPSF.main(gv.NSIDE,Energy = CTB_en_bins[bin_Number],use_generic_sigma = False,just_sigma = True)*np.pi/180
    

#     #load diffuse-background model
#     diff_fits = fits.open(gv.ferm_diff2_dir+'diffmodelsp15_ultraclean_Q2-front.fits')
#     diff_master = diff_fits[1].data['COARSEBINS']

#     diffp6 = [hp.ud_grade(diff_master[0,0,::,i],gv.NSIDE)*(CTB_en_bins[i+1-gv.CTB_en_min]-CTB_en_bins[i-gv.CTB_en_min])*4*np.pi/gv.NPIX*10**(3) for i in range(gv.CTB_en_min,gv.CTB_en_max)] #was minus 1!
    
#     diffp7 = [hp.ud_grade(diff_master[0,1,::,i],gv.NSIDE)*(CTB_en_bins[i+1-gv.CTB_en_min]-CTB_en_bins[i-gv.CTB_en_min])*4*np.pi/gv.NPIX*10**(3) for i in range(gv.CTB_en_min,gv.CTB_en_max)]
    
#     #this is temporary for \ell = 30!!
#     # array_norm = np.array([1.056, 1.010,0.979,0.967,0.960,1.003,1.024])
#     # array_norm = np.array([1.032, 0.9875,0.9592,0.9626,0.956,1.0006,1.06,1.093])
#     diffp6_norm = [diffp6[i]*diff_array_norm[i] for i in range(np.size(diff_array_norm))]
#     if gv.set_p7_norm_1:
#         diff_array_norm = np.array([1., 1.,1.,1.,1.,1.,1.,1.])
#     diffp7_norm = [diffp7[i]*diff_array_norm[i] for i in range(np.size(diff_array_norm))]

#     if bin_Number == 'All':
#         diffp6T=np.sum(diffp6_norm,axis=0) #summed p6 background model
#         diffp7T=np.sum(diffp7_norm,axis=0) #summed p7 background model
#     else:
#         diffp6T=diffp6[bin_Number]
#         diffp7T=diffp7[bin_Number]

#     if not gv.extra_diff:
#         if gv.boolp6==True:
#             print 'using p6 diffuse model'
#             x_diff_PSF_um = diffp6T
#         else:
#             print 'using p7 diffuse model'
#             x_diff_PSF_um = diffp7T
#     else:
#         diff_fits = fits.open(gv.extra_diff_dir+gv.extra_diff_name)
#         diff_master = diff_fits[1].data['COARSEBINS']

#         diff_c0 = [hp.ud_grade(diff_master[0,0,::,i],gv.NSIDE)*(CTB_en_bins[i+1-gv.CTB_en_min]-CTB_en_bins[i-gv.CTB_en_min])*4*np.pi/gv.NPIX*10**(3) for i in range(gv.CTB_en_min,gv.CTB_en_max)] #was minus 1!

#         diff_c1 = [hp.ud_grade(diff_master[0,1,::,i],gv.NSIDE)*(CTB_en_bins[i+1-gv.CTB_en_min]-CTB_en_bins[i-gv.CTB_en_min])*4*np.pi/gv.NPIX*10**(3) for i in range(gv.CTB_en_min,gv.CTB_en_max)]

#         diff_c2 = [hp.ud_grade(diff_master[0,2,::,i],gv.NSIDE)*(CTB_en_bins[i+1-gv.CTB_en_min]-CTB_en_bins[i-gv.CTB_en_min])*4*np.pi/gv.NPIX*10**(3) for i in range(gv.CTB_en_min,gv.CTB_en_max)]

#         diff_total = np.array(diff_c0) + np.array(diff_c1) + np.array(diff_c2)

#         x_diff_PSF_um = np.sum(diff_total,axis=0)
#         x_c0_PSF_um = np.sum(diff_c0,axis=0)
#         x_c1_PSF_um = np.sum(diff_c1,axis=0)
#         x_c2_PSF_um = np.sum(diff_c2,axis=0)
#         logging.info('%s','using diffuse model ' + gv.extra_diff_name)

#         diff_fits_add = fits.open(gv.extra_diff_dir+gv.extra_diff_add)
#         diff_master_add = diff_fits_add[1].data['COARSEBINS']
#         diff_add = np.sum([[hp.ud_grade(diff_master_add[0,j,::,i],gv.NSIDE)*(CTB_en_bins[i+1-gv.CTB_en_min]-CTB_en_bins[i-gv.CTB_en_min])*4*np.pi/gv.NPIX*10**(3) for i in range(gv.CTB_en_min,gv.CTB_en_max)] for j in range(np.shape(diff_master_add)[1])],axis=0)
#         x_c3_PSF_um = np.sum(diff_add,axis=0)


#     # In[12]:

#     #load Fermi bubbles intensity maps integrated in CTB energy bins (from Wei's code)
# #    print gv.CTB_en_max - gv.CTB_en_min
# #    print np.shape(CTB_exposure_maps)
# #    print CTB_en_bins
# #    print CTB_en_bins[0:gv.CTB_en_max-gv.CTB_en_min]
#     fermi_bubs_count_maps_CTB_binned = np.load(gv.maps_dir + 'bubbles_intensity_maps.npy')[gv.CTB_en_min-1:gv.CTB_en_max-1]* CTB_exposure_maps * 4*np.pi/gv.NPIX

#     # x_iso_PSF_um_total = np.load(gv.maps_dir + 'iso_exposure_maps.npy')
#     # x_iso_pSF_um_range = x_iso_PSF_um_total[gv.CTB_en_min-1:gv.CTB_en_max-1]

#     if bin_Number == 'All':
#         #x_bubs_PSF_um = np.sum(fermi_bubs_count_maps_CTB_binned[gv.CTB_en_min-1:gv.CTB_en_max-2],axis=0)
#         x_bubs_PSF_um = np.sum(np.transpose(np.transpose(np.array(fermi_bubs_count_maps_CTB_binned))*bubs_array_norm),axis=0)
#         # x_iso_PSF_um_total = np.load(gv.maps_dir + 'iso_exposure_maps.npy')
#         # x_iso_pSF_um_range = x_iso_PSF_um_total[gv.CTB_en_min-1:gv.CTB_en_max-1]
#         #print 'x_iso size is ' + str(np.shape(x_iso_pSF_um_range))
#         x_iso_PSF_um = np.sum(np.transpose(iso_array_norm* np.transpose([exp/np.mean(exp) for exp in CTB_exposure_maps])),axis=0)
#         # if normalization == 'Null':
#         #     # x_iso_PSF_um = np.sum(np.transpose(np.transpose(x_iso_pSF_um_range)*iso_array_norm),axis=0)/np.mean(np.sum(x_iso_pSF_um_range*iso_array_norm,axis=0))
#         #     # x_iso_PSF_um = np.sum(np.transpose(iso_array_norm* np.transpose([x_iso/np.mean(x_iso) for x_iso in x_iso_pSF_um_range])),axis=0)
#         #     x_iso_PSF_um = np.sum(np.transpose(iso_array_norm* np.transpose([exp/np.mean(exp) for exp in CTB_exposure_maps])),axis=0)
#         # else:
#         #     x_iso_PSF_um = np.sum(np.transpose(iso_array_norm* np.transpose([x_iso/np.mean(x_iso) for x_iso in x_iso_pSF_um_range])),axis=0)
#         # #x_iso_PSF_um = np.load(gv.maps_dir + 'iso_exposure_map.npy')
#     else:
#         x_bubs_PSF_um = fermi_bubs_count_maps_CTB_binned[bin_Number]
#         x_iso_PSF_um =  CTB_exposure_maps[bin_Number]/np.mean(CTB_exposure_maps[bin_Number])
#     #load isotropic map (normalized to have mean 1)
    
#     #====================================
#     #Now do NFW
    
#     #======================================
#     #construct NFW template: differential intensity maps for NFW components and integrated in CTB energy bins

#     from scipy import integrate, interpolate

#     def photons_per_ann(en_array):
#         #returns dN/dE per annihilation at values of en_array, using Bergstrom fit to spectrum
#         return 2./m_chi*0.42*np.exp(-8.*en_array/m_chi)/((en_array/m_chi)**1.5 + .00014)

#     rho_0_factor = np.power(r_0/r_s, gamma_nfw) * np.power(1. + r_0/r_s, 3. - gamma_nfw)

#     def rho_dimless(r):
#         #r in kpc
#         return np.power(r/r_s, -gamma_nfw) * np.power(1. + r/r_s, gamma_nfw - 3.) * rho_0_factor

#     def r_NFW(l, psi_deg):
#         return np.sqrt(r_0**2. + l**2. - 2.*r_0*l*np.cos(np.radians(psi_deg)))

#     def L_NFW_integral(psi_deg):
#         return integrate.quad(lambda l: rho_dimless(r_NFW(l, psi_deg))**2., 0., 100.*r_s)[0]

#     def make_NFW_intensity_map(en_array):
#         GC_lng=gv.GC_lng
#         psi_deg = np.arange(0., 180.5, 0.5/4) #BS: put in /4
#         intensity_NFW = sigma_v * rho_0**2. / (8. * np.pi * m_chi**2.) * 3.08E-5 * np.vectorize(L_NFW_integral)(psi_deg)
#         intensity_NFW_interp = interpolate.interp1d(psi_deg, intensity_NFW)
#         GC_vec = [np.cos(np.deg2rad(GC_lng)), np.sin(np.deg2rad(GC_lng)), 0.]
#         psi_deg_pixels = np.array([np.degrees(np.arccos(np.dot(GC_vec, hp.pix2vec(gv.NSIDE, pix)))) for pix in range(gv.NPIX)])
#         return np.outer(photons_per_ann(en_array), intensity_NFW_interp(psi_deg_pixels))

#     #make NFW intensity maps integrated in CTB energy bins
#     #use fermi_en (same energies as Fermi diffuse model) to evaluate dN/dE
#     # NFW_int_maps = make_NFW_intensity_map(fermi_en, lng)
#     # NFW_int_maps_CTB_binned, NFW_int_maps_fermi_binned = df.rebin_fermi_int(NFW_int_maps, fermi_en, CTB_en_bins)
#     #use CTB_en_bins_extended to evaluate dN/dE
#     CTB_en_bins_extended = CTB.get_CTB(gv.CTB_dir, gv.NSIDE, gv.CTB_en_min - 1, gv.CTB_en_max + 1)[0]
#     NFW_int_maps = make_NFW_intensity_map(CTB_en_bins_extended)
#     NFW_int_maps_CTB_binned = df.rebin_fermi_int(NFW_int_maps, CTB_en_bins_extended, CTB_en_bins)[0]


#     # In[47]:

#     #create xbg_PSF_um = unmasked healpix map of predicted mean diffuse background counts integrated over energy range
#     #from Fermi diffuse model components and CTBCORE exposure maps.
#     #this is done by multiplying binned intensity maps by binned CTB exposure maps,
#     #summing over bins,
#     #and smoothing using same PSF energy-weighted by excess.
#     #NOTE THAT THIS IS WRONG! SHOULD SMOOTH EACH BINNED INTENSITY MAP WITH PSF AT THAT ENERGY FIRST,
#     #THEN MULTIPLY BY EXPOSURE MAPS AND SUM!!!  TO DO!!!!!
#     #LATER NOTE: I GUESS WE GOT THIS FROM TRACY, FOR THE x_diff COMPONENT!!!

#     #first construct integrated, smoothed, exposure-corrected count maps for each component
#     def make_x_comp_PSF_um(comp_int_maps_CTB_binned, CTB_exposure_maps, sigma_PSF, npix, comp, **kwargs):
#         print 'doing component ', comp
#         if comp == 'NFW':
#             if normalization == 'Null':
#                 #for NFW component, apply mask_no_GC before smoothing
#                 x_comp_um = hp.ma(np.sum(comp_int_maps_CTB_binned * CTB_exposure_maps * 4*np.pi/npix, axis=0))
#                 x_comp_um.mask = kwargs['mask_for_NFW']
#                 if gv.smooth_NFW == False:
#                     print 'Not smoothing the NFW map!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
#                     x_comp_PSF_um = x_comp_um
#                 else:
#                     print 'Smoothing the NFW map and gv.smooth_NFW = ', gv.smooth_NFW
#                     with gv.suppress_stdout():
#                         x_comp_PSF_um = hp.smoothing(x_comp_um, sigma=sigma_PSF, verbose=False)
#             else:
#                 sigma_PSF_vec = [ CPSF.main(gv.NSIDE,Energy = en,use_generic_sigma = False,just_sigma = True)*np.pi/180 for en in CTB_en_bins ]
#                 logging.info('%s', 'The PSF_vec, used for smoothing NFW maps, is ' + str(sigma_PSF_vec))
#                 NFW_maps = np.transpose(NFW_array_norm*np.transpose(comp_int_maps_CTB_binned * CTB_exposure_maps * 4*np.pi/npix))
#                 with gv.suppress_stdout():
#                     NFW_maps_smooth =  [NFW_maps[i] for i in range(np.shape(NFW_maps)[0]) ]
#                     # NFW_maps_smooth = [ hp.smoothing(NFW_maps[i], sigma=sigma_PSF_vec[i], verbose=False) for i in range(np.shape(NFW_maps)[0]) ]
#                 x_comp_PSF_um = np.sum(NFW_maps_smooth,axis=0)
#         else:
#             x_comp_um = np.sum(comp_int_maps_CTB_binned * CTB_exposure_maps * 4*np.pi/npix, axis=0)
#             x_comp_PSF_um = x_comp_um
            
#         return x_comp_PSF_um

#     def make_x_comp_PSF_um_bin(comp_int_maps_CTB_binned, CTB_exposure_maps, sigma_PSF, npix, comp,bin_Number, **kwargs):
#         if comp == 'NFW':
#             #for NFW component, apply mask_no_GC before smoothing
#             x_comp_um = hp.ma(comp_int_maps_CTB_binned[bin_Number] * CTB_exposure_maps[bin_Number] * 4*np.pi/npix)
#             x_comp_um.mask = kwargs['mask_for_NFW']
#         else:
#             x_comp_um = comp_int_maps_CTB_binned[bin_Number] * CTB_exposure_maps[bin_Number] * 4*np.pi/npix
#         if comp == 'bubs' or comp == 'iso':
#             #don't smooth bubbles or isotropic
#             x_comp_PSF_um = x_comp_um
#         else:
#             with gv.suppress_stdout():
#                 #NOTE TO BEN: IF YOU TAKE OUT "REGRESSION" OPTION, MAKE SURE YOU ADD BACK MEAN IF NEEDED!
#                 #x_comp_PSF_um = hp.smoothing(x_comp_um, sigma=sigma_PSF, verbose=False, regression=False)
#                 if gv.smooth_NFW == False:
#                     print 'Not smoothing the NFW map!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
#                     x_comp_PSF_um = x_comp_um
#                 else:
#                     print 'Smoothing the NFW map and gv.smooth_NFW = ', gv.smooth_NFW
#                     x_comp_PSF_um = hp.smoothing(x_comp_um, sigma=sigma_PSF, verbose=False)
#         return x_comp_PSF_um

# #   x_bubs_PSF_um = make_x_comp_PSF_um(fermi_bubs_int_maps_CTB_binned, CTB_exposure_maps, gv.sigma_PSF, gv.NPIX, 'bubs')

# # x_iso_PSF_um = make_x_comp_PSF_um(np.outer(fermi_iso_int_CTB_binned, np.ones(NPIX)), CTB_exposure_maps,
# #                                   sigma_PSF, NPIX, 'iso')
# #
# #x_iso_PSF_um = np.load(maps_dir + 'iso_exposure_map.npy') #normalized to mean = 1, so A_iso parameter
# #now interpreted as counts x_iso --- change priors accordingly!!
#     #mask to remove GC (for removing singularity when smoothing NFW template map)
#     deg_GC = 0.5#0.5
#     mask_no_GC = masks.mask_ring(0, deg_GC, 90, gv.GC_lng, gv.NSIDE)
#     if bin_Number == 'All':
#         x_NFW_PSF_um = make_x_comp_PSF_um(NFW_int_maps_CTB_binned, CTB_exposure_maps, sigma_PSF, gv.NPIX, 'NFW', mask_for_NFW=mask_no_GC)
#     else:
#         x_NFW_PSF_um = make_x_comp_PSF_um_bin(NFW_int_maps_CTB_binned, CTB_exposure_maps, sigma_PSF, gv.NPIX, 'NFW',bin_Number, mask_for_NFW=mask_no_GC)
    
    
    
#     #change this later
#     #  x_NFW_PSF_um = x_iso_PSF_um
#     #full x_PSF
#     if not gv.extra_diff:
#         x_PSF = (x_diff_PSF_um, x_bubs_PSF_um, x_iso_PSF_um, x_NFW_PSF_um)
#     else:
#         if gv.extra_ring:
#             x_PSF = (x_c0_PSF_um,x_c1_PSF_um,x_c2_PSF_um, x_c3_PSF_um,x_bubs_PSF_um, x_iso_PSF_um, x_NFW_PSF_um)
#             hp.mollview(x_c3_PSF_um, min=0, max = 10, title='Inner Ring Extra Profile')
#             plt.savefig(gv.plots_dir + gv.file_tag + '/'  + 'c3_profile.pdf')
#             plt.close()
#         else:
#             x_PSF = (x_c0_PSF_um,x_c1_PSF_um,x_c2_PSF_um, x_bubs_PSF_um, x_iso_PSF_um, x_NFW_PSF_um)
#         hp.mollview(x_c0_PSF_um, min=0, max = 10, title='c0 Profile')
#         plt.savefig(gv.plots_dir + gv.file_tag + '/'  + 'c0_profile.pdf')
#         plt.close()
#         hp.mollview(x_c1_PSF_um, min=0, max = 10, title='c1 Profile')
#         plt.savefig(gv.plots_dir + gv.file_tag + '/'  + 'c1_profile.pdf')
#         plt.close()
#         hp.mollview(x_c2_PSF_um, min=0, max = 10, title='c2 Profile')
#         plt.savefig(gv.plots_dir + gv.file_tag + '/'  + 'c2_profile.pdf')
#         plt.close()


#     hp.mollview(x_NFW_PSF_um, min=0, max = 10, title='NFW Profile')
#     plt.savefig(gv.plots_dir + gv.file_tag + '/'  + 'NFW_profile.pdf')
#     plt.close()

#     hp.mollview(x_NFW_PSF_um, min=0, max = 1, title='NFW Profile (Outer)')
#     plt.savefig(gv.plots_dir + gv.file_tag + '/'  + 'NFW_profile_outer.pdf')
#     plt.close()
    

#     return x_PSF

