import sys, os

import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt

import healpy as hp
import mpmath as mp
from astropy.io import fits

import pulsars.masks as masks

import logging
from contextlib import contextmanager
import sys, os

from scipy import integrate, interpolate

#from numba import autojit
import time

size = np.size

class smooth_gaussian_psf_quick:
	def __init__(self,sigma,the_map):
		self.the_smooth_map = hp.smoothing(the_map, sigma= sigma*2*np.pi/360)
		where_vec = np.where(self.the_smooth_map < 0)[0]
		self.the_smooth_map[where_vec] = 0

# class smooth_save_fits:
# 	def __init__(self,fits_file_path):

# 		f = fits.open(fits_file_path)



class smooth_gaussian_psf:
	def __init__(self,sigma,the_map,mult_sigma_for_smooth = 10):
		self.sigma = sigma*2*np.pi/360 #convert to radians
		self.the_map = the_map
		self.mult_sigma_for_smooth = mult_sigma_for_smooth
		self.npix = size(the_map)
		self.nside = hp.npix2nside(self.npix)

		self.pixel_number_array = np.arange(self.npix)

		self.init_theta_phi_array()

		#self.smooth_the_map()


	def smooth_the_map(self):
		self.the_smooth_map = np.sum( [self.smooth_the_pixel(p) for p in self.pixel_number_array],axis=0)

	def init_theta_phi_array(self):
		self.theta_phi_array = hp.pix2ang(self.nside,self.pixel_number_array) 


	#@autojit
	def smooth_the_pixel(self,pix_num_center):
		smoothed_pixel_map = np.zeros(self.npix)

		# masks.mask_ring(0, 10*self.sigma, lat, lng, self.nside)

		theta_center,phi_center = hp.pix2ang(self.nside,pix_num_center)
		# theta,phi = hp.pix2ang(self.nside,pix_num)
		lat = theta_center*360/(2*np.pi)
		lng = phi_center*360/(2*np.pi)
		mask_where = masks.mask_ring(0, self.mult_sigma_for_smooth*self.sigma*360/(2*np.pi),lat, lng, self.nside)
		mask_pixel_vals = np.where(mask_where == 1)[0]
		mask_pixel_vals = np.array(list(set(list(mask_pixel_vals) + [pix_num_center])))

		# smoothed_pixel_map = np.vectorize(self.gaussian_func)(pix_num_center,self.pixel_number_array)
		
		smoothed_pixel_map[mask_pixel_vals] = np.vectorize(self.gaussian_func)(pix_num_center,mask_pixel_vals)

		return self.the_map[pix_num_center]*smoothed_pixel_map/np.sum(smoothed_pixel_map) #normalize to one

	#@autojit

	def smooth_the_pixel_by_angles(self,angle_vec,counts_vec): #angle vec of the form [ [theta,phi], [theta,phi], ... ]
		smoothed_pixel_map = np.zeros(self.npix)


		#theta_center,phi_center = hp.pix2ang(self.nside,pix_num_center)
		# theta_center=angle_vec[::,0]
		# phi_center=angle_vec[::,1]
		for theta_center,phi_center, count in map(None,angle_vec[0],angle_vec[1],counts_vec):
			#ta = time.time()
			#theta_center,phi_center = angle
			lat = theta_center*360/(2*np.pi)
			lng = phi_center*360/(2*np.pi)

			mask_where = masks.mask_ring(0, self.mult_sigma_for_smooth*self.sigma*360/(2*np.pi),lat, lng, self.nside)
			#mask_where=np.ones(self.npix)
			#print 'start: ', time.time() - ta
			mask_pixel_vals = np.where(mask_where == 1)[0]
			pix_num_center=hp.ang2pix(self.nside,theta_center,phi_center)
			mask_pixel_vals = np.array(list(set(list(mask_pixel_vals) + [pix_num_center])))
			#print 'The mask_pixel_vals are ',mask_pixel_vals

			#t0 = time.time()
			vals = np.vectorize(self.gaussian_func_ang)(theta_center,phi_center,mask_pixel_vals)
			#print 'step vals: ', time.time() - t0

			#print 'The vals are ', vals
			#print 'Sum of vals: ', np.sum(vals)
			if np.sum(vals)==0:
				smoothed_pixel_map[pix_num_center] += count
			else:
				smoothed_pixel_map[mask_pixel_vals] += count*vals/np.sum(vals)

			#print 'total: ', time.time() - ta

		return smoothed_pixel_map #normalize to one


	def gaussian_func(self,pix_num_center, pix_num):
		r = self.find_r(pix_num_center,pix_num)
		val = (1/(2.*np.pi*float(self.sigma)**2))*np.exp(-r**2/(2*self.sigma**2))
		return val

	def gaussian_func_ang(self,theta_center,phi_center, pix_num):
		r = self.find_r_ang(theta_center,phi_center,pix_num)
		val = (1/(2.*np.pi*float(self.sigma)**2))*np.exp(-r**2/(2*self.sigma**2))
		return val

	def total_gaussian_func(self,pix_num):
		return np.sum(self.the_map*self.gaussian_func(self.pixel_number_array,pix_num))
		# return np.sum(map( lambda pix_num_center : self.the_map[pix_num_center]*self.gaussian_func(pix_num_center,pix_num), self.pixel_number_array))

	def return_smoothed_maps(self):
		self.the_smooth_map = np.vectorize(self.total_gaussian_func)(self.pixel_number_array)
		return self.the_smooth_map

	def r_to_gauss(self,r):
		val = np.exp(-np.square(r)/(2*self.sigma**2)) #(1/(2.*np.pi*float(self.sigma)**2))*
		return val


	def make_r_matrix(self,vec):
		theta_center,phi_center = np.array(self.theta_phi_array[0]), np.array(self.theta_phi_array[1])
		theta,phi = np.transpose(np.mat(self.theta_phi_array[0][vec])), np.transpose(np.mat(self.theta_phi_array[1][vec]))
		self.r_matrix = np.sqrt( np.square(theta - theta_center) + np.square(self.dist_phi(phi,phi_center)) )
		self.g_matrix =  self.r_to_gauss( self.r_matrix )
		self.g_matrix_norm = self.g_matrix*1/np.sum(self.g_matrix,axis=1) #*np.transpose(np.mat(self.the_map[0:n])) 
		#print self.g_matrix
		#print np.shape(self.g_matrix)
		self.map_array = np.sum(np.array(self.g_matrix_norm)*self.the_map,axis=0)
		#self.g_matrix*np.transpose(np.mat(self.the_map))
		#self.map_array = np.array(np.transpose(   self.g_matrix*np.transpose(np.mat(self.the_map[0:n]))))[0]
		return self.map_array

	def find_r_ang(self,theta_center,phi_center,pix_num):
		theta,phi = self.theta_phi_array[0][pix_num], self.theta_phi_array[1][pix_num]
		return np.sqrt( (theta - theta_center)**2 + self.dist_phi(phi,phi_center)**2 )

	def find_r(self,pix_num_center,pix_num):
		# theta_center,phi_center = hp.pix2ang(self.nside,pix_num_center)
		# theta,phi = hp.pix2ang(self.nside,pix_num)
		theta_center,phi_center = self.theta_phi_array[0][pix_num_center], self.theta_phi_array[1][pix_num_center]
		theta,phi = self.theta_phi_array[0][pix_num], self.theta_phi_array[1][pix_num]
		return np.sqrt( (theta - theta_center)**2 + self.dist_phi(phi,phi_center)**2*(np.sin((theta+theta_center)/2.))**2 ) 

	# def find_r_vec(self,pix_num_center,pix_num):
	# 	# theta_center,phi_center = hp.pix2ang(self.nside,pix_num_center)
	# 	# theta,phi = hp.pix2ang(self.nside,pix_num)
	# 	theta_center,phi_center = np.array(self.theta_phi_array[0][pix_num_center]), np.array(self.theta_phi_array[1][pix_num_center])
	# 	theta,phi = self.theta_phi_array[0][pix_num], self.theta_phi_array[1][pix_num]
	# 	return np.sqrt( (theta - theta_center)**2 + self.dist_phi(phi,phi_center)**2 )

	def dist_phi(self,phi1,phi2):
		#print "phi1:",phi1
		#print "phi2:",phi2
		#print "dist:",np.minimum(np.abs(phi1-phi2), 2*np.pi-np.abs((phi1-phi2)))
		return np.minimum(np.abs(phi1-phi2), 2*np.pi-np.abs((phi1-phi2)))

		#np.sqrt( (theta - theta_center)**2 + np.sin((theta - theta_center))**2 *(phi - phi_center)**2  )


