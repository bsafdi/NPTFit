import numpy as np
import healpy as hp
import pulsars.masks as masks


class smooth_gaussian_psf_quick:
	def __init__(self,sigma,the_map):
		self.the_smooth_map = hp.smoothing(the_map, sigma= sigma*2*np.pi/360)
		where_vec = np.where(self.the_smooth_map < 0)[0]
		self.the_smooth_map[where_vec] = 0


class smooth_gaussian_psf:
	def __init__(self,sigma,the_map,mult_sigma_for_smooth = 10):
		self.sigma = sigma*2*np.pi/360 #convert to radians
		self.the_map = the_map
		self.mult_sigma_for_smooth = mult_sigma_for_smooth
		self.npix = np.size(the_map)
		self.nside = hp.npix2nside(self.npix)

		self.pixel_number_array = np.arange(self.npix)

		self.theta_phi_array = hp.pix2ang(self.nside,self.pixel_number_array)

	def smooth_the_map(self):
		self.the_smooth_map = np.sum( [self.smooth_the_pixel(p) for p in self.pixel_number_array],axis=0)

	def smooth_the_pixel(self,pix_num_center):
		smoothed_pixel_map = np.zeros(self.npix)

		theta_center,phi_center = hp.pix2ang(self.nside,pix_num_center)
		lat = theta_center*360/(2*np.pi)
		lng = phi_center*360/(2*np.pi)
		mask_where = masks.mask_ring(0, self.mult_sigma_for_smooth*self.sigma*360/(2*np.pi),lat, lng, self.nside)
		mask_pixel_vals = np.where(mask_where == 1)[0]
		mask_pixel_vals = np.array(list(set(list(mask_pixel_vals) + [pix_num_center])))

		smoothed_pixel_map[mask_pixel_vals] = np.vectorize(self.gaussian_func)(pix_num_center,mask_pixel_vals)

		return self.the_map[pix_num_center]*smoothed_pixel_map/np.sum(smoothed_pixel_map) #normalize to one


	def gaussian_func(self,pix_num_center, pix_num):
		r = self.find_r(pix_num_center,pix_num)
		val = (1/(2.*np.pi*float(self.sigma)**2))*np.exp(-r**2/(2*self.sigma**2))
		return val

	def find_r(self,pix_num_center,pix_num):
		theta_center,phi_center = self.theta_phi_array[0][pix_num_center], self.theta_phi_array[1][pix_num_center]
		theta,phi = self.theta_phi_array[0][pix_num], self.theta_phi_array[1][pix_num]
		return np.sqrt( (theta - theta_center)**2 + self.dist_phi(phi,phi_center)**2*(np.sin((theta+theta_center)/2.))**2 ) 

	def dist_phi(self,phi1,phi2):
		return np.minimum(np.abs(phi1-phi2), 2*np.pi-np.abs((phi1-phi2)))
