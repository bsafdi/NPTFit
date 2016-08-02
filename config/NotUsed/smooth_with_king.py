# This is a function setup to smooth a pixel with a king function
# This can then be used to smooth a map pixel by pixel - but doing this isn't a good idea as it is extremely slow. Thus it is just recommended for use with smoothing a pixel
# Nick Rodd, MIT

import sys, os
import numpy as np
import healpy as hp
import pulsars.masks as masks
import fermi.PSF_king_class as PSFKC
from astropy.io import fits

class smooth_king_psf:
    def __init__(self,maps_dir,the_map,Ebin,eventclass,eventtype):
	self.psf_data_dir = maps_dir + 'psf_data/'
	self.the_map = the_map
	self.Ebin = Ebin
	self.eventclass = eventclass
	self.eventtype = eventtype
	self.npix = np.size(the_map)
	self.nside = hp.npix2nside(self.npix)
	self.pixel_number_array = np.arange(self.npix)
	self.init_theta_phi_array()
	# Define the various parameters needed to load the king function parameters

    def init_theta_phi_array(self):
	self.theta_phi_array = hp.pix2ang(self.nside,self.pixel_number_array) 

    def smooth_the_pixel(self,pix_num_center):
	output_map = self.smooth_a_quartile(pix_num_center,1)
	if ((self.eventtype==4) or (self.eventtype==5) or (self.eventtype==0)):
	    output_map = (output_map + self.smooth_a_quartile(pix_num_center,2))/2.
	    if ((self.eventtype==5) or (self.eventtype==0)):
		output_map = (2*output_map + self.smooth_a_quartile(pix_num_center,3))/3.
		if self.eventtype==0:
		    output_map = (3*output_map + self.smooth_a_quartile(pix_num_center,4))/4.

	return output_map

    def smooth_a_quartile(self,pix_num_center,quartile):
	# Smooth a pixel for a given quartile
	# Load the parameters needed
	if quartile==1:
	    params_index_psf=10
	    rescale_index_psf=11
	    if self.eventclass==2:
		theta_norm_psf=[0.0000000,9.7381019e-06,0.0024811595,0.022328802,0.080147663,0.17148392,0.30634315,0.41720551]
	    if self.eventclass==5:
		theta_norm_psf=[0.0000000,9.5028121e-07,0.00094418357,0.015514370,0.069725775,0.16437751,0.30868705,0.44075016]
	if quartile==2:
	    params_index_psf=7
            rescale_index_psf=8
            if self.eventclass==2:
                theta_norm_psf=[0.0000000,0.00013001938,0.010239333,0.048691643,0.10790632,0.18585539,0.29140913,0.35576811]
            if self.eventclass==5:
                theta_norm_psf=[0.0000000,1.6070284e-05,0.0048551576,0.035358049,0.091767466,0.17568974,0.29916159,0.39315185]
	if quartile==3:
	    params_index_psf=4
            rescale_index_psf=5
            if self.eventclass==2:
                theta_norm_psf=[0.0000000,0.00074299273,0.018672204,0.062317201,0.12894928,0.20150553,0.28339386,0.30441893]
            if self.eventclass==5:
                theta_norm_psf=[0.0000000,0.00015569366,0.010164870,0.048955837,0.11750811,0.19840060,0.29488095,0.32993394]
	if quartile==4:
	    params_index_psf=1
            rescale_index_psf=2
            if self.eventclass==2:
                theta_norm_psf=[4.8923139e-07,0.011167475,0.092594658,0.15382001,0.16862869,0.17309118,0.19837774,0.20231968]
            if self.eventclass==5:
                theta_norm_psf=[0.0000000,0.0036816313,0.062240006,0.14027030,0.17077023,0.18329804,0.21722594,0.22251374]
	if self.eventclass==2:
	    psf_file_name = 'psf_P8R2_SOURCE_V6_PSF.fits'
	elif self.eventclass==5:
	    psf_file_name = 'psf_P8R2_ULTRACLEANVETO_V6_PSF.fits'
	fits_file_name = self.psf_data_dir + psf_file_name
	f = fits.open(fits_file_name)

	# Finally convert from an energy bin to an energy in GeV
	energyvalarray = np.array([0.336605534472,0.423761259986,0.533483821671,0.67161633,0.84551489,1.0644402,1.3400508,1.6870239,2.1238373,2.6737528,3.3660554,4.2376127,5.3348382,6.7161634,8.4551486,10.644402,13.400508,16.870240,21.238374,26.737528,33.660553,42.376126,53.348382,67.161633,84.551489,106.44402,134.00508,168.70240,212.38374,267.37528])
	energyval = energyvalarray[self.Ebin]

	# Now load up the parameters
	kparam=PSFKC.PSF_king(f, theta_norm=theta_norm_psf, rescale_index=rescale_index_psf, params_index=params_index_psf)

	fcore=kparam.return_king_params(energyval,'fcore')
	score=kparam.return_king_params(energyval,'score')
	gcore=kparam.return_king_params(energyval,'gcore')
	stail=kparam.return_king_params(energyval,'stail')
	gtail=kparam.return_king_params(energyval,'gtail')
	SpE=kparam.rescale_factor(energyval)

	# Note using the ring mask we only calculate the PSF in a region where it
	# is appreciably non-zero so we don't have to calculate it for the full sky
	smoothed_pixel_map = np.zeros(self.npix)
	theta_center,phi_center = hp.pix2ang(self.nside,pix_num_center)
	lat = theta_center
	lng = phi_center
	# Note we convert the radius to degrees as this is the format mask_ring takes in
	mask_where = masks.mask_ring(0, 22*SpE*((score+stail)/2.)*180./np.pi,lat*180./np.pi, lng*180./np.pi, self.nside)
	mask_pixel_vals = np.where(mask_where == 1)[0]
	mask_pixel_vals = np.array(list(set(list(mask_pixel_vals) + [pix_num_center])))
	smoothed_pixel_map[mask_pixel_vals] = np.vectorize(self.king_func)(pix_num_center,mask_pixel_vals,fcore,score,gcore,stail,gtail,SpE)

	return self.the_map[pix_num_center]*smoothed_pixel_map/np.sum(smoothed_pixel_map) #normalize smoothing template to the value of the map at the pixel location

    def king_fn_base(self,x,sigma,gamma):
        # A basic king function
        return (1/(2*np.pi*sigma**2))*(1-1/gamma)*(1+x**2/(2*gamma*sigma**2))**(-gamma)

    def king_fn_full(self,r,fcore,score,gcore,stail,gtail,SpE):
        # The combination of two king functions relevant for the Fermi PSF
        return fcore*self.king_fn_base(r/SpE,score,gcore)+(1-fcore)*self.king_fn_base(r/SpE,stail,gtail)

    def king_func(self,pix_num_center,pix_num,fcore,score,gcore,stail,gtail,SpE):
	# Evaluate the king function at every point
	r = self.find_r(pix_num_center,pix_num)
	return self.king_fn_full(r,fcore,score,gcore,stail,gtail,SpE) 

    def find_r(self,pix_num_center,pix_num):
	# This function finds the distance between a given pixel (pix num) and
	# a central pixel (pix_num_center) on a sphere
	theta_center,phi_center = self.theta_phi_array[0][pix_num_center], self.theta_phi_array[1][pix_num_center]
	theta,phi = self.theta_phi_array[0][pix_num], self.theta_phi_array[1][pix_num]
	return self.dist_sphere(theta_center,theta,phi_center,phi) 
	
    def dist_sphere(self,theta1,theta2,phi1,phi2):
	# This is a function to calculate the angular distance on a sphere, assuming that
	# the two points are close
	# For details (and the full formula): https://en.wikipedia.org/wiki/Great-circle_distance
	return np.sqrt( (theta1 - theta2)**2 + self.dist_phi(phi1,phi2)**2*(np.sin((theta1+theta2)/2.))**2 )
	
    def dist_phi(self,phi1,phi2):
	# Determine the shortest distance in phi accounting for the periodicity
        return np.minimum(np.abs(phi1-phi2), 2*np.pi-np.abs((phi1-phi2)))
