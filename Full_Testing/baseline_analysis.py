import sys, os
import numpy as np
from mpmath import mp

import matplotlib.pyplot as plt
import matplotlib as mpl
import healpy as hp

import triangle

import pulsars.masks as masks

import config.compute_PSF as CPSF
import config.config_file as cf
import config.bayesian_scan as bs
from config.NPTF import NPTF
import config.bayesian_scan_models as bsm
#from config.bayesian_scan_models import compute_poissonian_spectrum as cps
from astropy.io import fits

from config.analyze_results import analyze_results as ar
from scipy.interpolate import interp1d
from scipy.integrate import quad
import fermi.fermi_plugin as fp
import analysis_classes.mask_for_analysis as mfa

import pandas as pd

import pickle

class mask_for_analysis:
    def __init__(self, nside = 128, plane_mask = True, band_mask_range = [-30,30], lcut = False, lmin = -20, lmax = 20, bcut = False, bmin = -20, bmax = 20, mask_ring = False, lat = 90, lng = 0,inner = 0, outer = 30, input_mask = False, the_input_mask = 'Null'):

        self.nside = nside
        
        self.set_mask_total(plane_mask = plane_mask,band_mask_range = band_mask_range,lcut = lcut, lmin = lmin,lmax = lmax,bcut = bcut,bmin = bmin,bmax = bmax,mask_ring = mask_ring,lat=lat,lng = lng,inner = inner, outer= outer,input_mask = input_mask, the_input_mask = the_input_mask)

    def make_plane_mask(self, band_mask_range):
        plane_mask_array = masks.mask_lat_band(90+band_mask_range[0], 90+band_mask_range[1], self.nside)
        return plane_mask_array

    def make_long_mask(self, lmin, lmax):
        long_mask_array = masks.mask_not_long_band(lmin, lmax, self.nside)
        return long_mask_array

    def make_lat_mask(self, bmin, bmax):
        lat_mask_array = masks.mask_not_lat_band(90+bmin, 90+bmax, self.nside)
        return lat_mask_array

    def make_ring_mask(self, inner, outer, lat, lng):
        ring_mask = np.logical_not(masks.mask_ring(inner, outer, lat, lng, self.nside))
        return ring_mask

    def mask_and_compress(self,the_map,mask):
        map_masked = hp.ma(the_map)
        map_masked.mask = mask
        return map_masked.compressed()

    def set_mask_total(self, nside = 128, plane_mask = True, band_mask_range = [-30,30], lcut = True, lmin = -20, lmax = 20, bcut = True, bmin = -20, bmax = 20, mask_ring = False, lat = 90, lng = 0,inner = 0, outer = 30, input_mask = False, the_input_mask = 'Null'):
        
        self.plane_mask_ar = self.make_plane_mask(band_mask_range)
        self.long_mask_ar = self.make_long_mask(lmin, lmax)
        self.lat_mask_ar = self.make_lat_mask(bmin, bmax)
        self.ring_mask = self.make_ring_mask(inner, outer, lat, lng)
    
        if input_mask:
            self.mask_total = the_input_mask
        else:
            # Need to initialise array, just add the plane mask, if not called then remove it
            mask_array = self.plane_mask_ar

            if not plane_mask:
                mask_array -= self.plane_mask_ar

            if lcut:
                mask_array += self.long_mask_ar

            if bcut:
                mask_array += self.lat_mask_ar

            if mask_ring:
                mask_array += self.ring_mask
                pixel_num_center = hp.ang2pix(self.nside, lat*2*np.pi/360,lng*2*np.pi/360)
                if inner==0 and band_mask_range[0] == 0. and band_mask_range[1] == 0. and mask_array[pixel_num_center] == True :
                    print 'setting inner pixel to False ...'
                    mask_array[pixel_num_center] = False

            mask_geom_total = mask_array

            self.mask_total = mask_geom_total


class compute_spectra(mfa.mask_for_analysis):
    def __init__(self,b,run_tag,En_min,En_max,edep=False,polychord=False,minuit=False,minuit_new=False,ebin=0,*args,**kwargs):
        self.minuit=minuit
        self.run_tag = run_tag
        self.b = b
        if not minuit:
            self.b.load_scan(run_tag = self.run_tag,edep=edep,polychord=polychord,minuit=minuit_new)
            self.samples = self.b.samples
        self.edep = edep
        self.ebin=ebin

        self.nside = self.b.nside

        self.En_min = En_min
        self.En_max = En_max
        self.En_center = 10**((np.log10(self.En_max)+np.log10(self.En_min))/2)
        self.dE = self.En_max - self.En_min

        self.scan_intensity_array_dict = {}

        #self.set_mask_total(*args,**kwargs)
        mfa.mask_for_analysis.__init__(self,nside=self.nside,*args,**kwargs)
        self.make_mask_total_dict()


    def make_mask_total_dict(self):
        b = self.b
        self.mask_total_dict = {comp: self.mask_total for comp in b.poiss_models.keys()}



    ##########################
    ###not automatic

    def return_template_intensity(self,mask,comp):
        b = self.b
        if self.edep:
            summed_template = b.templates_dict_nested[comp]['templates'][self.ebin]
            mean_exposure = [b.exposure_map][self.ebin]
        else:
            summed_template = b.templates_dict_nested[comp]['summed_templates_not_compressed']
            mean_exposure = np.mean([b.exposure_map],axis=0)

        summed_template_masked_compressed = self.mask_and_compress(summed_template,mask)
        mean_exposure_masked_compressed = self.mask_and_compress(mean_exposure,mask)

        self.st = summed_template_masked_compressed
        self.se = mean_exposure_masked_compressed

        template_intensity = np.mean(summed_template_masked_compressed/mean_exposure_masked_compressed/hp.nside2pixarea(self.nside))

        return template_intensity

    def calculate_poiss_flux(self,mask,comp,log_prior=False):
        b = self.b
        template_intensity = self.return_template_intensity(mask,comp)
        if not self.minuit:
            samples=self.samples
        else:
            samples=[b.medians]
        if not log_prior:
            scan_intensity_array = template_intensity*np.array([samples[i][self.b.poiss_comp_numbers[comp]] for i in range(len(samples))])
        else:
            scan_intensity_array = template_intensity*np.array([10**samples[i][self.b.poiss_comp_numbers[comp]] for i in range(len(samples))])
        return scan_intensity_array

    def return_poiss_intensity_array(self,comp):
        b = self.b #intensity, dE, units [photons/cm^2/s/sr]
        mask = self.mask_total_dict[comp]
        scan_intensity_array = self.calculate_poiss_flux(mask,comp,log_prior = self.b.poiss_models[comp]['log_prior'])
        return scan_intensity_array

    def return_poiss_intensity_conf_int(self,comp,quant = [0.16,0.5,0.84]):
        scan_intensity_array = self.return_poiss_intensity_array(comp)
        self.scan_intensity_array_dict[comp] = scan_intensity_array
        print 'at comp ', comp, 'log prior: ', self.b.poiss_models[comp]['log_prior'], 'intensity: ' ,triangle.quantile(scan_intensity_array,quant)
        return triangle.quantile(scan_intensity_array,quant)

    # def return_poiss_intensity_conf_int(self,comp,quant = [0.16,0.5,0.84]):
    #     b = self.b #intensity, dE, units [photons/cm^2/s/sr]
    #     mask = self.mask_total_dict[comp]
    #     scan_intensity_array = self.calculate_poiss_flux(mask,comp,log_prior = self.b.poiss_models[comp]['log_prior'])
    #     print 'at comp ', comp, 'log prior: ', self.b.poiss_models[comp]['log_prior'], 'intensity: ' ,triangle.quantile(scan_intensity_array,quant)
    #     return triangle.quantile(scan_intensity_array,quant)

    def return_spectrum_conf_int(self,comp,quant  = [0.16,0.5,0.84]):
        self.spectrum = self.En_center**2*np.array(self.return_poiss_intensity_conf_int(comp, quant = quant))/self.dE
        return self.spectrum

    def return_spectrum_minuit(self,comp):
        self.spectrum = self.En_center**2*np.array(self.return_poiss_intensity_minuit(comp))/self.dE
        return self.spectrum

    #####################
    ###functions to be called

    def make_spectra_dict(self):
        self.spectra_dict = {comp: self.return_spectrum_conf_int(comp) for comp in self.b.poiss_models.keys()}

    def make_norm_dict(self):
        self.norm_dict = self.b.norms_poiss # {comp: np.mean(np.array([self.samples[i][self.b.poiss_comp_numbers[comp]] for i in range(len(self.samples))])) for comp in self.b.poiss_models.keys() }

    def save_spectra_dict(self,spect_file,emin,emax,over_write= False):
        if os.path.isfile(spect_file + '.npy') and not over_write:
            print 'appending to spectrum ...'
            spect = list(np.load(spect_file + '.npy'))
            spect += [[ self.En_center, self.spectra_dict]]
            np.save(spect_file,np.array(spect))
        else:
            print 'saving spectra for first time'
            spect = [[ self.En_center, self.spectra_dict]]
            np.save(spect_file,np.array(spect))

    def save_norm_dict(self,norm_file,emin,emax,over_write= False):
        if os.path.isfile(norm_file + '.npy') and not over_write:
            print 'appending to norm array ...'
            norm = list(np.load(norm_file + '.npy'))
            norm += [[ self.En_center, self.norm_dict]] 
            print 'The norm array is now, ', np.array(norm) 
            np.save(norm_file,np.array(norm))
        else:
            print 'saving norm array for first time'
            norm = [[ self.En_center, self.norm_dict]]
            np.save(norm_file,np.array(norm))

class compute_edep_spectra(mfa.mask_for_analysis):
    def __init__(self,b,run_tag,En_bins,polychord=False,minuit=False,minuit_new=False,*args,**kwargs):
        self.minuit=minuit
        self.polychord = polychord
        self.minuit_new=minuit_new
        self.args = args
        self.kwargs = kwargs

        self.run_tag = run_tag
        self.b = b
        if not self.minuit:
            self.b.load_scan(run_tag = self.run_tag,edep=True,polychord=polychord,minuit=minuit_new)
            self.samples = self.b.samples

        self.nside = self.b.nside

        self.scan_intensity_array_dict_array = [{} for ens in En_bins]

        self.En_min = []
        self.En_max = []
        self.En_center = []
        self.dE = []
        for i in range(len(En_bins)-1):
            self.En_min.append(En_bins[i])
            self.En_max.append(En_bins[i+1])
            self.En_center.append( 10**((np.log10(self.En_max[i])+np.log10(self.En_min[i]))/2) )
            self.dE.append( self.En_max[i] - self.En_min[i]  )

        #self.set_mask_total(*args,**kwargs)
        mfa.mask_for_analysis.__init__(self,nside=self.nside,*args,**kwargs)
        self.make_compute_spectra_array()

    def make_compute_spectra_array(self):
        self.cs_array = [ compute_spectra(self.b,self.run_tag,self.En_min[i],self.En_max[i],edep=True,polychord=self.polychord,minuit=self.minuit,minuit_new=self.minuit_new,ebin=i,*self.args,**self.kwargs) for i in range(len(self.En_min)) ]
        for cs in self.cs_array:
            cs.make_mask_total_dict()

    def return_spectrum_conf_int(self,comp,quant  = [0.16,0.5,0.84]):
        spect = []
        for en in range(len(self.En_min)):
            spect.append(self.cs_array[en].return_spectrum_conf_int(comp+'-'+str(en),quant=quant))
        return np.array(spect)

    def log_to_normal(self,array,is_log):
        array_2 = []
        for i in range(len(array)):
            if is_log[i]:
                array_2.append(10**array[i])
            else:
                array_2.append(array[i])
        return array_2

    # def return_non_poiss_samples(self,comp,sample,ebin):
    #     b = self.b
    #     self.n_poiss = b.n_poiss
    #     self.model_decompression_non_poiss = np.array(b.model_decompression_key[b.n_poiss:])
    #     model_where = np.where(self.model_decompression_non_poiss[:,0]==comp)[0]+b.n_poiss
    #     is_log_prior = np.array(b.model_decompression_key)[model_where][:,1]
    #     is_log_prior =  list(is_log_prior=='True')
    #     samples_model_not_log = self.log_to_normal(np.array(sample)[model_where],is_log_prior)

    #     #print samples_model_not_log
    #     samples_model_not_log_fixed = samples_model_not_log[0:3] + [samples_model_not_log[3+ebin]] 
    #     samples_model_not_log_fixed[0] = samples_model_not_log_fixed[0]/samples_model_not_log_fixed[-1]

    #     # print self.model_decompression_non_poiss
    #     # print is_log_prior
    #     # print b.model_decompression_key
    #     # print model_where

    #     return samples_model_not_log_fixed

    def return_non_poiss_samples(self,comp,sample,ebin):
        b = self.b
        self.n_poiss = b.n_poiss
        self.model_decompression_non_poiss = np.array(b.model_decompression_key[b.n_poiss:])
        model_where = np.where(self.model_decompression_non_poiss[:,0]==comp)[0]+b.n_poiss
        is_log_prior = np.array(b.model_decompression_key)[model_where][:,1]
        is_log_prior =  list(is_log_prior=='True')
        samples_model_not_log = self.log_to_normal(np.array(sample)[model_where],is_log_prior)
        if self.b.nbreak==1:
            samples_model_not_log_fixed = samples_model_not_log[0:3] + [samples_model_not_log[3+ebin]] 
        elif self.b.nbreak==2:
            samples_model_not_log_fixed = samples_model_not_log[0:5] + [samples_model_not_log[5+ebin]] 
        elif self.b.nbreak==3:
            samples_model_not_log_fixed = samples_model_not_log[0:7] + [samples_model_not_log[7+ebin]]
        samples_model_not_log_fixed[0] = samples_model_not_log_fixed[0]/samples_model_not_log_fixed[-1]

        return samples_model_not_log_fixed

    # def dnds(self,comp,sample,s,ebin): #over full region
    #     # print 'new sample!'
    #     # print sample
    #     samples_reduced = self.return_non_poiss_samples(comp,sample,ebin)
    #     # print samples_reduced

    #     #print start_pos,end_pos,samples_reduced_log,samples_reduced

    #     APS, n1,n2,Sb = samples_reduced
    #     dnds = np.zeros(len(s))
    #     where_vec = np.where(s < Sb)[0]
    #     not_where_vec = np.where(s >= Sb)[0]
    #     dnds[where_vec] = APS*(s[where_vec]/Sb)**(-n2)
    #     dnds[not_where_vec] = APS*(s[not_where_vec]/Sb)**(-n1)
    #     # if s < Sb:
    #     #     dnds = APS*(s/Sb)**(-n2)
    #     # else:
    #     #     dnds = APS*(s/Sb)**(-n1)
    #     # # print APS, n1, n2, Sb, dnds
    #     return dnds

    def dnds(self,comp,sample,s,ebin): #over full region
        samples_reduced = self.return_non_poiss_samples(comp,sample,ebin)

        if self.b.nbreak==1:
            APS, n1,n2,Sb = samples_reduced
            #APS = APS/Sb
            dnds = np.zeros(len(s))
            where_vec =np.where(s < Sb)[0]
            not_where_vec = np.where(s >= Sb)[0]
            dnds[where_vec] = APS*(s[where_vec]/Sb)**(-n2)
            dnds[not_where_vec] = APS*(s[not_where_vec]/Sb)**(-n1)
        elif self.b.nbreak==2:
            APS, n1,n2, n3, lam, Sb2 = samples_reduced
            if self.b.dnds_model_array[0] == 'specify_breaks':
                Sb1 = lam
            else:
                Sb1 = Sb2*lam
            #APS = APS #/Sb2
            dnds = np.zeros(len(s))
            where_vec_1 =np.where(s < Sb1)[0]
            where_vec_2 =np.where( (s >= Sb1) & (s < Sb2) )[0]
            not_where_vec = np.where(s >= Sb2)[0]
            dnds[where_vec_1] = APS*(s[where_vec_1]/Sb1)**(-n3)
            dnds[where_vec_2] = APS*(s[where_vec_2]/Sb1)**(-n2)
            dnds[not_where_vec] = (Sb2 / Sb1)**-n2*APS*(s[not_where_vec]/Sb2)**(-n1)
        elif self.b.nbreak==3:
            APS, n1,n2, n3,n4, lam1,lam2, Sb3 = samples_reduced
            if self.b.dnds_model_array[0] == 'specify_breaks':
                Sb2 = lam2
                Sb1 = lam1
            else:
                Sb2 = Sb3*lam2
                Sb1 = Sb2*lam1
            #APS = APS #/Sb2
            dnds = np.zeros(len(s))
            where_vec_1 =np.where(s < Sb1)[0]
            where_vec_2 =np.where( (s >= Sb1) & (s < Sb2) )[0]
            where_vec_3 =np.where( (s >= Sb2) & (s < Sb3) )[0]
            not_where_vec = np.where(s >= Sb3)[0]
            dnds[where_vec_1] = APS*(s[where_vec_1]/Sb1)**(-n4)
            dnds[where_vec_2] = APS*(s[where_vec_2]/Sb1)**(-n3)
            dnds[where_vec_3] = (Sb1 / Sb2)**n3*APS*(s[where_vec_3]/Sb2)**(-n2)
            dnds[not_where_vec] = (Sb2 / Sb3)**n2*(Sb1 / Sb2)**n3*APS*(s[not_where_vec]/Sb3)**(-n1)


        return dnds


    def return_template_intensity(self,mask,comp,ebin):
        b = self.b
        summed_template = b.templates_dict_nested[comp]['summed_templates_not_compressed']
        mean_exposure = [b.exposure_map][ebin]

        summed_template_masked_compressed = self.mask_and_compress(summed_template,mask)
        mean_exposure_masked_compressed = self.mask_and_compress(mean_exposure,mask)

        template_intensity = np.mean(summed_template_masked_compressed/mean_exposure_masked_compressed/hp.nside2pixarea(self.nside))

        return template_intensity

    def return_intensity_arrays_non_poiss(self,comp,ebin,smin=0.00001,smax=10000,nsteps=10000,quant=[0.16,0.5,0.84]):
        b = self.b
        self.template_instensity = self.return_template_intensity(self.mask_total,comp,ebin)

        #source count
        self.sarray = 10**np.linspace(np.log10(smin),np.log10(smax),nsteps)#np.arange(smin,smax,float((smax-smin)/float(nsteps)))
        self.ds = [self.sarray[i+1]-self.sarray[i] for i in range(len(self.sarray)-1)]
        self.ds = np.array( self.ds + [self.ds[-1]] )
        #self.dnds_array = np.array([self.dnds(comp_number, sam,s) for sam in self.samples])
        if not self.minuit:
            self.intensity_array_non_poiss = map( lambda sam : np.sum(self.template_instensity*self.dnds(comp, sam,self.sarray,ebin)*self.sarray*self.ds ) , b.samples )
        else:
            ###for now, only sample is medians for minuit scan
            self.intensity_array_non_poiss = map( lambda sam : np.sum(self.template_instensity*self.dnds(comp, sam,self.sarray,ebin)*self.sarray*self.ds ) , [b.medians] )

        self.scan_intensity_array_dict_array[ebin][comp] = self.intensity_array_non_poiss

        return triangle.quantile(self.intensity_array_non_poiss,quant)

    def return_non_poiss_spectrum_conf_int(self,comp,smin=0.00001,smax=10000,nsteps=10000, quant  = [0.16,0.5,0.84]):
        self.spectrum = []
        for en in range(len(self.En_min)):
            comp_en = comp #+ '-'+str(en)
            self.spectrum.append( self.En_center[en]**2*np.array(self.return_intensity_arrays_non_poiss(comp_en, en,smin=smin,smax=smax,nsteps=nsteps, quant = quant))/self.dE[en]  )
        return np.array(self.spectrum)

    def make_spectra_dict(self,smin=0.00001,smax=10000,nsteps=10000,quant = [0.16,0.5,0.84]):
        self.spectra_dict = [ [self.En_center[en], {comp: self.return_spectrum_conf_int(comp,quant=quant)[en] for comp in self.b.poiss_models_bare_keys}] for en in range(len(self.En_center))]

        for comp in self.b.non_poiss_models.keys():
            spect_comp_non_poiss = self.return_non_poiss_spectrum_conf_int(comp,smin=smin,smax=smax,nsteps=nsteps, quant  = quant)
            print spect_comp_non_poiss
            for en in range(len(self.En_center)):
                self.spectra_dict[en][1][comp+'-np']=np.array(spect_comp_non_poiss[en])

    def compute_total_spectra(self,quant = [0.16,0.5,0.84]):
        self.total_spectrum_array = [] #np.zeros( len(self.En_center) )
        self.spect_debug_dict = {}
        for en in range(len(self.En_center)):
            total_array = np.zeros( len(self.cs_array[en].scan_intensity_array_dict[self.b.poiss_model_keys[0]])   )
            for comp in self.b.poiss_model_keys:
                self.spect_debug_dict[comp] = self.En_center[en]**2*self.cs_array[en].scan_intensity_array_dict[comp]/self.dE[en]
                total_array += self.En_center[en]**2*self.cs_array[en].scan_intensity_array_dict[comp]/self.dE[en]
            for comp in self.b.non_poiss_models.keys():
                total_array += self.En_center[en]**2*np.array(self.scan_intensity_array_dict_array[en][comp])/self.dE[en]
            self.total_spectrum_array += [triangle.quantile( total_array ,quant)]

    def compute_EGB_spectra(self,quant = [0.16,0.5,0.84]):
        self.EGB_spectrum_array = [] #np.zeros( len(self.En_center) )
        for en in range(len(self.En_center)):
            total_array = np.zeros( len(self.cs_array[en].scan_intensity_array_dict[self.b.poiss_model_keys[0]])   )
            comp = 'iso-0'
            total_array += self.En_center[en]**2*self.cs_array[en].scan_intensity_array_dict[comp]/self.dE[en]
            comp = 'iso'
            total_array += self.En_center[en]**2*np.array(self.scan_intensity_array_dict_array[en][comp])/self.dE[en]
            self.EGB_spectrum_array += [triangle.quantile( total_array ,quant)]

    def compute_PS_over_EGB_spectra(self,quant = [0.16,0.5,0.84]):
        self.PS_over_EGB_spectrum_array = [] #np.zeros( len(self.En_center) )
        for en in range(len(self.En_center)):
            total_array = np.zeros( len(self.cs_array[en].scan_intensity_array_dict[self.b.poiss_model_keys[0]])   )
            total_ps_array = np.zeros( len(self.cs_array[en].scan_intensity_array_dict[self.b.poiss_model_keys[0]])   )
            comp = 'iso-0'
            total_array += self.En_center[en]**2*self.cs_array[en].scan_intensity_array_dict[comp]/self.dE[en]
            comp = 'iso'
            total_array += self.En_center[en]**2*np.array(self.scan_intensity_array_dict_array[en][comp])/self.dE[en]
            total_ps_array += self.En_center[en]**2*np.array(self.scan_intensity_array_dict_array[en][comp])/self.dE[en]
            self.PS_over_EGB_spectrum_array += [triangle.quantile( np.array(total_ps_array)/np.array(total_array) ,quant)]

    def save_spectra_dict(self,spect_file):
        np.save(spect_file,np.array(self.spectra_dict))


class save_ps_template_info:
    def __init__(self,b,run_tag,edep=False,polychord=False,minuit=False,minuit_new=False):
        self.b = b
        if not minuit:
            self.b.load_scan(run_tag = run_tag,edep=edep,polychord=polychord,minuit=minuit_new)

    def save_ps_template(self,ps_tag,save_tag,**kwargs):
        b = self.b
        save_vec=[]
        save_vec.append(b.norms_non_poiss[ps_tag])
        save_vec.append(b.templates_dict_nested[ps_tag]['summed_templates_not_compressed'])
        # ps_results_dir = b.data_dir_base + 'ps_results/'
        # b.make_dirs([ps_results_dir])
        # with open(ps_results_dir + save_tag,'w') as thefile:
        #     for item in save_vec:
        #         thefile.write("%s\n" % item)

        if 'save_path' in kwargs.keys():
            save_path = kwargs['save_path']
        else:
            save_path = b.data_ps_dir + save_tag


        with open(save_path, 'wb') as f:
            pickle.dump(save_vec, f)
        #np.savetxt(ps_results_dir + save_tag,save_vec)




class make_triangle:
    def __init__(self,b,run_tag,edep=False,polychord=False,minuit=False,minuit_new=False):
        self.run_tag = run_tag
        self.b = b
        self.b.load_scan(run_tag = self.run_tag, edep=edep, polychord=polychord, minuit=minuit_new)

    def make_triangle(self,plot_name='triangle.pdf',return_fig = False,save_fig = True,the_range=None):
        plot_dir_for_triangle = self.b.plots_dir_for_run

        if the_range==None:
            self.the_range = [1 for i in range(self.b.n_params)]  #0.95 
        else:
            self.the_range = the_range

        self.tri_fig = triangle.corner(self.b.samples, labels=self.b.params,smooth=1.5,smooth1d=1,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_fmt='.2f',
        title_args={'fontsize':14},
        range=self.the_range,
        plot_datapoints=False,
        verbose=False)
        if save_fig:
            plt.savefig(plot_dir_for_triangle + plot_name)
        plt.close()

        if return_fig:
            return self.tri_fig

class plot_2FHL(mfa.mask_for_analysis):

    def __init__(self,ps_data_dir,nside,*args,**kwargs):

        self.nside = nside
        self.ps_data_dir = ps_data_dir

        self.load_2FHL_data()
        self.CTB_hi_en_bins = 10**np.linspace(np.log10(50), np.log10(2000),11)
        self.CTB_hi_en_bin_centers = 10**np.array([(np.log10(self.CTB_hi_en_bins[i])+np.log10(self.CTB_hi_en_bins[i+1]))/2. for i in range(len(self.CTB_hi_en_bins)-1)]) 


        mfa.mask_for_analysis.__init__(self,nside=self.nside,*args,**kwargs)

    def load_2FHL_data(self):

        # Load in the 2FHL catalog as downloaded
        # These have l and b in equatorial coordinates
        source_2fhl_df = pd.read_csv(self.ps_data_dir + '2fhl.dat', sep='|', comment='#')

        # Remove whitespace from appropriate columns
        source_2fhl_df.rename(columns=lambda x: x.strip(),inplace=True ) 
        for col in source_2fhl_df.columns.values:
            try: 
                source_2fhl_df[col] = source_2fhl_df[col].map(str.strip)
            except TypeError:
                continue

        # Convert to numeric data
        self.source_2fhl_df=source_2fhl_df.convert_objects(convert_numeric=True);

    def lb2pix(self, nside,l,b,nest=False):
        '''
        Function to convert right ascension and descent to galactic
        coordinates, then get index corresponding to healpix array
        '''
        return hp.ang2pix(nside,np.deg2rad(90-b),np.deg2rad(l),nest=nest) 

    def return_counts(self, CTB_hi_en_min = 0, CTB_hi_en_max = 9, flux_min = 7e-12, flux_max = 4e-10, flux_bins = 14,error_range = 0.68,cumul=False):
        """
        Return histogrammed source counts from 2FHL data
        """
        self.error_range = error_range
        self.area_factor = 1.
        if not self.mask_total == None: # Reduce dataframe to unmasked region
            to_include = []
            for i in range(len(self.source_2fhl_df)):
                source_pix = (self.lb2pix(self.nside, self.source_2fhl_df['lii'].values[i],self.source_2fhl_df['bii'].values[i]))
                if self.mask_total[source_pix] == 0: 
                    to_include.append(i)
            self.source_2fhl_df = self.source_2fhl_df.iloc[to_include]
            self.area_factor = 1. - np.sum(self.mask_total)/float(hp.nside2npix(self.nside))

        # If one of the standard 2FHL bins, use that flux, otherwise use provided power-law fit
        # if CTB_hi_en_min == 0 and CTB_hi_en_max == 3:
        #     flux_values_reduced = self.source_2fhl_df.flux50_171_gev.values
        # elif CTB_hi_en_min == 3 and CTB_hi_en_max == 6:
        #     flux_values_reduced = self.source_2fhl_df.flux171_585_gev.values
        # elif CTB_hi_en_min == 6 and CTB_hi_en_max == 10:
        #     flux_values_reduced = self.source_2fhl_df.flux585_2000_gev.values
        # else:
        flux_values = self.source_2fhl_df.flux.values
        gamma_values = self.source_2fhl_df.spectral_index.values
        flux_values_reduced = [flux*(self.CTB_hi_en_bins[CTB_hi_en_max]**(-gamma+1) - self.CTB_hi_en_bins[CTB_hi_en_min]**(-gamma+1))/(self.CTB_hi_en_bins[-1]**(-gamma+1) - self.CTB_hi_en_bins[0]**(-gamma+1)) for flux, gamma in zip(flux_values, gamma_values)]

        deg = 180/np.pi
        sr = 4*np.pi
        srdeg2 = sr*deg**2 #(180/np.pi)**2#

        counts, bin_edges = np.histogram(flux_values_reduced, bins=np.logspace(int(np.log10(flux_min)),int(np.log10(flux_max)), flux_bins))
        #plt.gca().set_xscale("log") #, patches

        bin_centres = 10**((np.log10(bin_edges[:-1]) + np.log10(bin_edges[1:]))/2.)
        bin_centers = bin_centres #british to american
        bin_width = bin_edges[1:]-bin_edges[:-1]
        #plt.close()

        x_counts = bin_centres
        y_counts = np.array(counts/(self.area_factor*bin_width*srdeg2))
        y_counts_err = np.array(np.sqrt(counts)/(self.area_factor*bin_width*srdeg2))
        if cumul:
            counts_cum = np.array( [np.sum(counts[i:]) for i in range(len(counts))])
            counts = counts_cum
            y_counts = np.array(counts_cum/(self.area_factor*srdeg2))
            y_counts_cum_err = np.array(np.sqrt(counts_cum)/(self.area_factor*bin_width*srdeg2))

        if self.error_range == 0.68:
            errors_sigma = np.loadtxt(self.ps_data_dir + 'ErrorTab68.dat')
        else:
            errors_sigma = np.loadtxt(self.ps_data_dir + 'ErrorTab95.dat')
        errors = errors_sigma[counts]
        error_L = []
        error_H = []
        for i in range(len(counts)):
            error_L.append(counts[i]-errors[i][0]-10**-8)
            error_H.append(errors[i][1]-counts[i])
        if cumul:
            self.error_L = np.array(error_L)/(self.area_factor*srdeg2)
            self.error_H = np.array(error_H)/(self.area_factor*srdeg2)
        else:
            self.error_L = np.array(error_L)/(self.area_factor*bin_width*srdeg2)
            self.error_H = np.array(error_H)/(self.area_factor*bin_width*srdeg2)

        self.x_errors_L = np.array([ bin_centers[i]-bin_edges[i] for i in range(np.size(bin_centers))])
        self.x_errors_H = np.array([ bin_edges[i+1]-bin_centers[i] for i in range(np.size(bin_centers))])
        
        # if cumul:
        #     return x_counts, y_counts_cum, self.error_L, self.error_H#y_counts_err
        # else:
        return x_counts, y_counts, self.error_L, self.error_H, self.x_errors_L, self.x_errors_H



class plot_3FGL(mfa.mask_for_analysis):
    def __init__(self,f,only_gal = False,CTB_en_min = 'False',CTB_en_max = 'False',*args,**kwargs):
        self.f = f
        self.nside = f.nside
        self.npix = hp.nside2npix(self.nside)
        energies = 10**np.linspace(np.log10(0.3),np.log10(300),31)
        if CTB_en_min != 'False' and CTB_en_max != 'False':
            self.CTB_en_min, self.CTB_en_max = np.where(np.round(energies,3)==np.round(CTB_en_min,3))[0][0], np.where(np.round(energies,3)==np.round(CTB_en_max,3))[0][0]
        else:
            self.CTB_en_min = f.CTB_en_min
            self.CTB_en_max = f.CTB_en_max
        print self.CTB_en_min, self.CTB_en_max
        # self.CTB_en_min_bin_number = np.where(f.CTB_en_bins == self.CTB_en_min)[0][0]
        # self.CTB_en_max_bin_number = np.where(f.CTB_en_bins == self.CTB_en_max)[0][0]
        self.CTB_en_min_bin_number = self.CTB_en_min-f.CTB_en_min
        self.CTB_en_max_bin_number = self.CTB_en_max-f.CTB_en_min
        print self.CTB_en_min_bin_number,self.CTB_en_max_bin_number


        self.ps_data_dir = f.maps_dir + '3FGL/'
        self.only_gal = only_gal

        mfa.mask_for_analysis.__init__(self,nside=self.nside,*args,**kwargs)
        self.total_exposure_map = np.mean(f.exposure_map[self.CTB_en_min_bin_number:self.CTB_en_max_bin_number],axis=0) ###this should be masked!
        # self.total_exposure_map = []
        # for exp in f.exposure_map:
        #     exp_masked = hp.ma(exp)
        #     exp_masked.mask = self.mask_total
        #     exp_masked_compressed = exp_masked.compressed()
        #     self.total_exposure_map += [np.mean(exp_masked_compressed)]
        # self.total_exposure_map = np.array(self.total_exposure_map)

        self.load_3FGL_data()
        self.make_3FGL_predicted_counts_map()
        self.mask_3FGL_predicted_counts_map()
        self.mask_compress_exposure_map()
        self.make_3FGL_predicted_flux_map()


    def load_3FGL_data(self):
        if self.only_gal:
            psc_tag = 'sphCoordsGAL.dat'
            psc_flux_tag = 'fluxTabGALBin.dat'
        else:
            psc_tag = 'sphCoords.dat'
            psc_flux_tag = 'fluxTabBin.dat'

        psc_3FGL_name = self.ps_data_dir+ psc_tag
        psc_3FGL_flux_name = self.ps_data_dir + psc_flux_tag

        self.psc_3FGL = np.loadtxt(psc_3FGL_name)
        self.psc_3FGL_flux = np.loadtxt(psc_3FGL_flux_name)

        self.psc_3FGL_flux_bin = np.sum(self.psc_3FGL_flux[::,self.CTB_en_min:self.CTB_en_max],axis=1)


    def make_3FGL_predicted_counts_map(self):
        f = self.f
        self.psArray = np.zeros(self.npix)
        self.pix_num_array = [hp.ang2pix(self.nside, self.psc_3FGL[i,0], self.psc_3FGL[i,1]) for i in range(np.shape(self.psc_3FGL)[0]) ]
        self.counts_vec = np.array(self.psc_3FGL_flux_bin)*f.total_exposure_map[np.array(self.pix_num_array)]
        self.psArray[self.pix_num_array] = self.counts_vec


    def mask_3FGL_predicted_counts_map(self):
        self.ps_array_masked = hp.ma(self.psArray)
        self.ps_array_masked.mask = self.mask_total
        self.ps_array_compressed = self.ps_array_masked.compressed()
        self.ps_array_non_zero = self.ps_array_compressed[np.where(self.ps_array_compressed)] 

    def mask_compress_exposure_map(self):
        f = self.f
        self.area_mask = np.sum(1-self.mask_total)*hp.nside2pixarea(self.nside)*(360/(2.*np.pi))**2 #in deg^2
        self.exp_masked = hp.ma(self.total_exposure_map)
        self.exp_masked.mask = self.mask_total
        self.exp_masked_compressed = self.exp_masked.compressed()
        self.exp_masked_mean = np.mean(self.exp_masked_compressed)

    def make_3FGL_predicted_flux_map(self):
        self.mask_compress_exposure_map()
        self.ps_flux_masked_compressed = self.ps_array_non_zero/self.exp_masked_compressed[np.where(self.ps_array_compressed)]


    def make_3FGL_flux_histogram(self,use_counts=False):
        if use_counts:
            yvals, bins = np.histogram(self.ps_array_non_zero,bins=10**np.linspace(np.log10(self.flux_min), np.log10(self.flux_max), self.n_flux_bins));
            bins = np.array(bins) / np.mean(self.exp_masked_compressed[np.where(self.ps_array_compressed)])
        else:
            yvals, bins = np.histogram(self.ps_flux_masked_compressed,bins=10**np.linspace(np.log10(self.flux_min), np.log10(self.flux_max), self.n_flux_bins));
        bin_centers = np.array([10**((np.log10(bins[i])+np.log10(bins[i+1]))/2) for i in range(np.shape(bins)[0]-1)])
        bin_width = np.array([bins[i+1]-bins[i] for i in range(np.size(bins)-1)])

        rescale = 1/self.area_mask
        self.yvals_rescale = yvals*rescale/bin_width
        self.bin_centers_rescale = bin_centers

        if self.error_range == 0.68:
            errors_sigma = np.loadtxt(self.ps_data_dir + 'ErrorTab68.dat')
        else:
            errors_sigma = np.loadtxt(self.ps_data_dir + 'ErrorTab95.dat')
        errors = errors_sigma[yvals]
        error_L = []
        error_H = []
        for i in range(np.shape(yvals)[0]):
            error_L.append(yvals[i]-errors[i][0]-10**-8)
            error_H.append(errors[i][1]-yvals[i])
        self.error_L = np.array(error_L)*rescale/bin_width
        self.error_H = np.array(error_H)*rescale/bin_width
            
        self.x_errors_L = np.array([ bin_centers[i]-bins[i] for i in range(np.size(bins)-1)])
        self.x_errors_H = np.array([ bins[i+1]-bin_centers[i] for i in range(np.size(bins)-1)])

        self.yvals = yvals
        self.bins = bins

        
    #user executes these functions

    def configure_3FGL(self,flux_min=10**-11,flux_max=10**-8,n_flux_bins = 25,error_range = 0.68,**kwargs):
        if 'use_counts' in kwargs.keys():
            use_counts = kwargs['use_counts']
        else:
            use_counts=False
        if 'only_gal' in kwargs.keys():
            self.only_gal = kwargs['only_gal']
        self.flux_min = flux_min
        self.flux_max = flux_max
        self.n_flux_bins = n_flux_bins
        self.error_range = error_range


        self.load_3FGL_data()
        self.make_3FGL_predicted_counts_map()
        self.mask_3FGL_predicted_counts_map()
        self.make_3FGL_predicted_flux_map()
        self.make_3FGL_flux_histogram(use_counts=use_counts)


    def plot_3FGL(self,spow=0,*args,**kwargs):
        # plt.errorbar(self.bin_centers_rescale,self.yvals_rescale*self.bin_centers_rescale**spow,xerr=[self.x_errors_L,self.x_errors_H],yerr=[self.error_L*self.bin_centers_rescale**spow,self.error_H*self.bin_centers_rescale**spow],*args,**kwargs)
        flux_factor = ((self.bin_centers_rescale-self.x_errors_L)*(self.bin_centers_rescale+self.x_errors_H))**(spow/2.)
        plt.errorbar(self.bin_centers_rescale,self.yvals_rescale*flux_factor,xerr=[self.x_errors_L,self.x_errors_H],yerr=[self.error_L*flux_factor,self.error_H*flux_factor],*args,**kwargs)




# class Memoize:
#     def __init__(self, f):
#         self.f = f
#         self.memo = {}
#     def __call__(self, *args):
#         if not args in self.memo:
#             self.memo[args] = self.f(*args)
#         return self.memo[args]
def memoize(f):
    """ Memoization decorator for functions taking one or more arguments. """
    class memodict(dict):
        def __init__(self, f):
            self.f = f
        def __call__(self, *args):
            return self[args]
        def __missing__(self, key):
            ret = self[key] = self.f(*key)
            return ret
    return memodict(f)


class plot_dnds(mfa.mask_for_analysis):
    def __init__(self,b,*args,**kwargs):
        self.b = b
        self.nside = b.nside
        self.npix = hp.nside2npix(self.nside)

        mfa.mask_for_analysis.__init__(self,nside=self.nside,*args,**kwargs)

    def mask_compress_exposure_map(self):
        b = self.b
        self.area_mask = np.sum(1-self.mask_total)*hp.nside2pixarea(self.nside)*(360/(2.*np.pi))**2 #in deg^2
        self.exp_masked = hp.ma(b.total_exposure_map)
        self.exp_masked.mask = self.mask_total
        self.exp_masked_compressed = self.exp_masked.compressed()
        self.exp_masked_mean = np.mean(self.exp_masked_compressed)

    def log_to_normal(self,array,is_log):
        array_2 = []
        for i in range(len(array)):
            if is_log[i]:
                array_2.append(10**array[i])
            else:
                array_2.append(array[i])
        return array_2

    def return_non_poiss_samples(self,comp,sample):
        b = self.b
        self.n_poiss = b.n_poiss
        self.model_decompression_non_poiss = np.array(b.model_decompression_key[b.n_poiss:])
        model_where = np.where(self.model_decompression_non_poiss[:,0]==comp)[0]+b.n_poiss
        is_log_prior = np.array(b.model_decompression_key)[model_where][:,1]
        is_log_prior =  list(is_log_prior=='True')
        samples_model_not_log = self.log_to_normal(np.array(sample)[model_where],is_log_prior) 

        # print self.model_decompression_non_poiss
        # print is_log_prior
        # print b.model_decompression_key
        # print model_where

        return samples_model_not_log


    def return_non_poiss_samples_edep(self,comp,sample,ebin):
        b = self.b
        self.n_poiss = b.n_poiss
        self.model_decompression_non_poiss = np.array(b.model_decompression_key[b.n_poiss:])
        model_where = np.where(self.model_decompression_non_poiss[:,0]==comp)[0]+b.n_poiss
        is_log_prior = np.array(b.model_decompression_key)[model_where][:,1]
        is_log_prior =  list(is_log_prior=='True')
        samples_model_not_log = self.log_to_normal(np.array(sample)[model_where],is_log_prior)
        if self.b.nbreak==1:
            samples_model_not_log_fixed = samples_model_not_log[0:3] + [samples_model_not_log[3+ebin]] 
        elif self.b.nbreak==2:
            samples_model_not_log_fixed = samples_model_not_log[0:5] + [samples_model_not_log[5+ebin]] 
        elif self.b.nbreak==3:
            samples_model_not_log_fixed = samples_model_not_log[0:7] + [samples_model_not_log[7+ebin]]

        samples_model_not_log_fixed[0] = samples_model_not_log_fixed[0]/samples_model_not_log_fixed[-1]

        return samples_model_not_log_fixed


    @memoize
    def dnds_calc(APS,n1,n2,Sb,s):
        if s < Sb:
            dnds = APS*(s/Sb)**(-n2)
        else:
            dnds = APS*(s/Sb)**(-n1)
        return dnds

    def dnds(self,comp,sample,s): #over full region
        # print 'new sample!'
        # print sample
        samples_reduced = self.return_non_poiss_samples(comp,sample)
        # print samples_reduced

        #print start_pos,end_pos,samples_reduced_log,samples_reduced

        APS, n1,n2,Sb = samples_reduced
        if s < Sb:
            dnds = APS*(s/Sb)**(-n2)
        else:
            dnds = APS*(s/Sb)**(-n1)
        # print APS, n1, n2, Sb, dnds
        return dnds

    def dnds_edep_reduced(self,comp,sample_reduced,s,ebin):
        if self.b.nbreak==1:
            APS, n1,n2,Sb = sample_reduced[0:3] + [sample_reduced[3+ebin]]
            APS = APS/Sb
            dnds = np.zeros(len(s))
            where_vec =np.where(s < Sb)[0]
            not_where_vec = np.where(s >= Sb)[0]
            dnds[where_vec] = APS*(s[where_vec]/Sb)**(-n2)
            dnds[not_where_vec] = APS*(s[not_where_vec]/Sb)**(-n1)
        elif self.b.nbreak==2:
            APS, n1,n2, n3, lam, Sb2 = sample_reduced[0:5] + [sample_reduced[5+ebin]]
            if self.b.dnds_model_array[0] == 'specify_breaks':
                Sb1 = lam
            else:
                Sb1 = Sb2*lam
            APS = APS/Sb2
            dnds = np.zeros(len(s))
            where_vec_1 =np.where(s < Sb1)[0]
            where_vec_2 =np.where( (s >= Sb1) & (s < Sb2) )[0]
            not_where_vec = np.where(s >= Sb2)[0]
            dnds[where_vec_1] = APS*(s[where_vec_1]/Sb1)**(-n3)
            dnds[where_vec_2] = APS*(s[where_vec_2]/Sb1)**(-n2)
            dnds[not_where_vec] = (Sb2 / Sb1)**-n2*APS*(s[not_where_vec]/Sb2)**(-n1)
        elif self.b.nbreak==3:
            APS, n1,n2, n3,n4, lam1,lam2, Sb3 = sample_reduced[0:7] + [sample_reduced[7+ebin]]
            if self.b.dnds_model_array[0] == 'specify_breaks':
                Sb2 = lam2
                Sb1 = lam1
            else:
                Sb2 = Sb3*lam2
                Sb1 = Sb2*lam1
            #APS = APS #/Sb2
            dnds = np.zeros(len(s))
            where_vec_1 =np.where(s < Sb1)[0]
            where_vec_2 =np.where( (s >= Sb1) & (s < Sb2) )[0]
            where_vec_3 =np.where( (s >= Sb2) & (s < Sb3) )[0]
            not_where_vec = np.where(s >= Sb3)[0]
            dnds[where_vec_1] = APS*(s[where_vec_1]/Sb1)**(-n4)
            dnds[where_vec_2] = APS*(s[where_vec_2]/Sb1)**(-n3)
            dnds[where_vec_3] = (Sb1 / Sb2)**n3*APS*(s[where_vec_3]/Sb2)**(-n2)
            dnds[not_where_vec] = (Sb2 / Sb3)**n2*(Sb1 / Sb2)**n3*APS*(s[not_where_vec]/Sb3)**(-n1)
        return dnds

    def dnds_edep(self,comp,sample,s,ebin): #over full region
        samples_reduced = self.return_non_poiss_samples_edep(comp,sample,ebin)

        # APS, n1,n2,Sb = samples_reduced

        # #return self.dnds_calc(APS,n1,n2,Sb,s)
        # dnds = np.zeros(len(s))
        # where_vec =np.where(s < Sb)[0]
        # # print s[where_vec]
        # # print APS*(s[where_vec]/Sb)
        # not_where_vec = np.where(s >= Sb)[0]
        # dnds[where_vec] = APS*(s[where_vec]/Sb)**(-n2)
        # dnds[not_where_vec] = APS*(s[not_where_vec]/Sb)**(-n1)
        # # if s < Sb:
        # #     dnds = APS*(s/Sb)**(-n2)
        # # else:
        # #     dnds = APS*(s/Sb)**(-n1)
        # # # print APS, n1, n2, Sb, dnds
        if self.b.nbreak==1:
            APS, n1,n2,Sb = samples_reduced
            #APS = APS/Sb
            dnds = np.zeros(len(s))
            where_vec =np.where(s < Sb)[0]
            not_where_vec = np.where(s >= Sb)[0]
            dnds[where_vec] = APS*(s[where_vec]/Sb)**(-n2)
            dnds[not_where_vec] = APS*(s[not_where_vec]/Sb)**(-n1)
        elif self.b.nbreak==2:
            APS, n1,n2, n3, lam, Sb2 = samples_reduced
            if self.b.dnds_model_array[0] == 'specify_breaks':
                Sb1 = lam
            else:
                Sb1 = Sb2*lam
            #APS = APS #/Sb2
            dnds = np.zeros(len(s))
            where_vec_1 =np.where(s < Sb1)[0]
            where_vec_2 =np.where( (s >= Sb1) & (s < Sb2) )[0]
            not_where_vec = np.where(s >= Sb2)[0]
            dnds[where_vec_1] = APS*(s[where_vec_1]/Sb1)**(-n3)
            dnds[where_vec_2] = APS*(s[where_vec_2]/Sb1)**(-n2)
            dnds[not_where_vec] = (Sb2 / Sb1)**-n2*APS*(s[not_where_vec]/Sb2)**(-n1)
        elif self.b.nbreak==3:
            APS, n1,n2, n3,n4, lam1,lam2, Sb3 = samples_reduced
            if self.b.dnds_model_array[0] == 'specify_breaks':
                Sb2 = lam2
                Sb1 = lam1
            else:
                Sb2 = Sb3*lam2
                Sb1 = Sb2*lam1
            #APS = APS #/Sb2
            dnds = np.zeros(len(s))
            where_vec_1 =np.where(s < Sb1)[0]
            where_vec_2 =np.where( (s >= Sb1) & (s < Sb2) )[0]
            where_vec_3 =np.where( (s >= Sb2) & (s < Sb3) )[0]
            not_where_vec = np.where(s >= Sb3)[0]
            dnds[where_vec_1] = APS*(s[where_vec_1]/Sb1)**(-n4)
            dnds[where_vec_2] = APS*(s[where_vec_2]/Sb1)**(-n3)
            dnds[where_vec_3] = (Sb1 / Sb2)**n3*APS*(s[where_vec_3]/Sb2)**(-n2)
            dnds[not_where_vec] = (Sb2 / Sb3)**n2*(Sb1 / Sb2)**n3*APS*(s[not_where_vec]/Sb3)**(-n1)


        return dnds

    def calculate_dnds_arrays(self,comp,smin=0.01,smax=1000,nsteps=1000,qs =[0.16, 0.5, 0.84],edep=False,ebin=0,spow=0 ,plot_sum=False):
        b = self.b
        self.mask_compress_exposure_map()

        #deal with the template
        # if edep:
        #     summed_template = b.templates_dict_nested[comp]['templates'][ebin]
        # else:
        summed_template = b.templates_dict_nested[comp]['summed_templates_not_compressed']
        #mean_exposure = np.mean(self.the_dict['exposure_map'],axis=0)
        summed_template_masked_compressed = self.mask_and_compress(summed_template,self.mask_total)
        self.template_sum = np.sum(summed_template_masked_compressed)

        #source count
        self.sarray = 10**np.linspace(np.log10(smin),np.log10(smax),nsteps) #np.arange(smin,smax,float((smax-smin)/float(nsteps)))
        #self.dnds_array = np.array([self.dnds(comp_number, sam,s) for sam in self.samples])
        self.flux_array = self.sarray/self.exp_masked_mean

        rf = self.template_sum*self.exp_masked_mean/self.area_mask
        if edep != True:
            qArray = [triangle.quantile(np.array([self.dnds(comp, sam,s) for sam in b.samples]),qs) for s in self.sarray]
            self.qmean = rf*np.array([np.mean(np.array([self.dnds(comp, sam,s) for sam in b.samples])) for s in self.sarray])
        else:
            # qArray = [triangle.quantile(np.array([self.dnds_edep(comp, sam,s,ebin) for sam in b.samples]),qs) for s in self.sarray]
            # self.qmean = rf*np.array([np.mean(np.array([self.dnds_edep(comp, sam,s,ebin) for sam in b.samples])) for s in self.sarray])
            data_array = np.array([self.flux_array**spow*self.dnds_edep(comp, sam,self.sarray,ebin) for sam in b.samples ])
            self.data_array = data_array
            #print np.shape(data_array)
            qArray = [triangle.quantile(data_array[::,i],qs) for i in range(len(self.sarray)) ]
            self.qmean = rf*np.array([np.mean(data_array[::,i]) for i in range(len(self.sarray)) ])

            if plot_sum:
                self.dfarray = [self.flux_array[i+1] - self.flux_array[i] for i in range(len(self.flux_array)-1)]
                self.dfarray += [self.dfarray[-1]]
                self.dfarray = np.array(self.dfarray)
                self.sum_array = np.array([  np.array( [ np.sum(data_array[sam,i:]*self.dfarray[i:]) for i in range(len(self.sarray)) ] ) for sam in range(len(b.samples))])
                #*self.flux_array[i:]

                qArray_sum = [triangle.quantile(self.sum_array[::,i],qs) for i in range(len(self.sarray)) ]
                self.qmean_sum = rf*np.array([np.mean(self.sum_array[::,i]) for i in range(len(self.sarray)) ])
            # qArray = [triangle.quantile(np.array([self.dnds_edep(comp, sam,s,ebin) for sam in b.samples]),qs) for s in self.sarray]
            # self.qmean = rf*np.array([np.mean(np.array([self.dnds_edep(comp, sam,s,ebin) for sam in b.samples])) for s in self.sarray])

        self.q16 = rf*np.array([q[0] for q in qArray]) #units:photons^-1 cm^2 s deg^-2
        self.q5 = rf*np.array([q[1] for q in qArray])
        self.q84 = rf*np.array([q[2] for q in qArray])

        if plot_sum:
            rf2 = self.template_sum/self.area_mask
            self.q16_sum = rf*np.array([q[0] for q in qArray_sum]) #units:photons^-1 cm^2 s deg^-2
            self.q5_sum = rf*np.array([q[1] for q in qArray_sum])
            self.q84_sum = rf*np.array([q[2] for q in qArray_sum])
            # Saving these so we can read them off the CTA limit
            np.save("q16_sum.npy", self.q16_sum)
            np.save("q5_sum.npy", self.q5_sum)
            np.save("q84_sum.npy", self.q84_sum)
            np.save("flux_array.npy", self.flux_array)

    def calculate_fixed_dnds_arrays(self,comp,smin=0.01,smax=1000,nsteps=1000,ebin=0,spow=0):
        medians=self.b.fixed_ps_model_dict[comp]['medians']
        summed_template = self.b.fixed_ps_model_dict[comp]['template']

        self.calculate_static_dnds_arrays(comp,summed_template,medians,smin=smin,smax=smax,nsteps=nsteps,ebin=ebin,spow=spow)

    def calculate_minuit_dnds_arrays(self,comp,smin=0.01,smax=1000,nsteps=1000,ebin=0,spow=0):
        medians=self.b.norms_non_poiss[comp]
        summed_template = self.b.templates_dict_nested[comp]['summed_templates_not_compressed']
        self.calculate_static_dnds_arrays(comp,summed_template,medians,smin=smin,smax=smax,nsteps=nsteps,ebin=ebin,spow=spow)

    def calculate_static_dnds_arrays(self,comp,summed_template,medians,smin=0.01,smax=1000,nsteps=1000,ebin=0,spow=0):
        b = self.b
        self.mask_compress_exposure_map()

        summed_template_masked_compressed = self.mask_and_compress(summed_template,self.mask_total)
        self.template_sum = np.sum(summed_template_masked_compressed)

        #source count
        self.sarray = 10**np.linspace(np.log10(smin),np.log10(smax),nsteps) #np.arange(smin,smax,float((smax-smin)/float(nsteps)))

        self.flux_array = self.sarray/self.exp_masked_mean

        rf = self.template_sum*self.exp_masked_mean/self.area_mask
        data = self.dnds_edep_reduced(comp, medians,self.sarray,ebin)*self.flux_array**spow
        self.data = data
        self.q = self.data*rf
        

    def plot_source_count_band(self,comp,smin=0.01,smax=1000,nsteps=1000,qs=[0.16,0.5,0.84],calculate_dnds=False,edep=False,ebin=0,spow=0,plot_sum=False,*args,**kwargs):
        if calculate_dnds:
            self.calculate_dnds_arrays(comp,smin=smin,smax=smax,nsteps=nsteps,qs=qs,edep=edep,ebin=ebin,spow=spow,plot_sum=plot_sum)
        if plot_sum:
            plt.fill_between(self.flux_array,self.q16_sum,self.q84_sum,*args,**kwargs)#
        else:        
            plt.fill_between(self.flux_array,self.q16,self.q84,*args,**kwargs)#facecolor='Chartreuse',interpolate=True, alpha=1,linewidth=0)
        #plt.plot(self.flux_array ,self.q5,'k--',label="NFW PS")

    def plot_source_count_median(self,comp,smin=0.01,smax=1000,nsteps=1000,qs=[0.16,0.5,0.84],calculate_dnds=False,edep=False,ebin=0,spow=0,plot_sum=False,*args,**kwargs):
        if calculate_dnds:
            self.calculate_dnds_arrays(comp,smin=smin,smax=smax,nsteps=nsteps,qs=qs,edep=edep,ebin=ebin,spow=spow,plot_sum=plot_sum)
        if plot_sum:
            plt.plot(self.flux_array ,self.q5_sum,*args,**kwargs)
        else:
            plt.plot(self.flux_array ,self.q5,*args,**kwargs)

    def plot_fixed_source_count(self,comp,smin=0.01,smax=1000,nsteps=1000,calculate_dnds=True,ebin=0,spow=0,*args,**kwargs):
        if calculate_dnds:
            self.calculate_fixed_dnds_arrays(comp,smin=smin,smax=smax,nsteps=nsteps,ebin=ebin,spow=spow)
        plt.plot(self.flux_array ,self.q,*args,**kwargs)

    def plot_minuit_source_count(self,comp,smin=0.01,smax=1000,nsteps=1000,calculate_dnds=True,ebin=0,spow=0,*args,**kwargs):
        if calculate_dnds:
            self.calculate_minuit_dnds_arrays(comp,smin=smin,smax=smax,nsteps=nsteps,ebin=ebin,spow=spow)
        plt.plot(self.flux_array ,self.q,*args,**kwargs)

    # def plot_source_count_mean(self,comp,smin=0.01,smax=1000,nsteps=1000,qs=[0.16,0.5,0.84],calculate_dnds=False,edep=False,ebin=0,spow=0,*args,**kwargs):
    #     if calculate_dnds:
    #         self.calculate_dnds_arrays(comp,smin=smin,smax=smax,nsteps=nsteps,qs=qs,edep=edep,ebin=ebin,spow=spow)
    #     plt.plot(self.flux_array ,self.qmean,*args,**kwargs)
    def plot_source_count_mean(self,comp,smin=0.01,smax=1000,nsteps=1000,qs=[0.16,0.5,0.84],calculate_dnds=False,edep=False,ebin=0,spow=0,plot_sum=False,*args,**kwargs):
        if calculate_dnds:
            self.calculate_dnds_arrays(comp,smin=smin,smax=smax,nsteps=nsteps,qs=qs,edep=edep,ebin=ebin,spow=spow,plot_sum=plot_sum)
        if plot_sum:
            plt.plot(self.flux_array ,self.qmean_sum,*args,**kwargs)
        else:
            plt.plot(self.flux_array ,self.qmean,*args,**kwargs)

class save_log_evidence:
    def __init__(self,b,edep=False,polychord=False):
        self.b = b
        self.b.load_scan(run_tag = b.run_tag,edep=edep,polychord=polychord)
        self.lge = self.b.s['nested sampling global log-evidence']
        self.lge_err = self.b.s['nested sampling global log-evidence error']
        print 'Log-global evidence ', [self.lge, self.lge_err]
        np.savetxt(b.plots_dir_for_run + 'lge.txt', [self.lge, self.lge_err])


class plot_flux_fraction(mfa.mask_for_analysis):
    def __init__(self,b,*args,**kwargs):
        self.b = b
        self.nside = b.nside
        self.npix = hp.nside2npix(self.nside)

        mfa.mask_for_analysis.__init__(self,nside=self.nside,*args,**kwargs)

    def mask_compress_exposure_map(self):
        b = self.b
        self.area_mask = np.sum(1-self.mask_total)*hp.nside2pixarea(self.nside)*(360/(2.*np.pi))**2 #in deg^2
        self.exp_masked = hp.ma(b.total_exposure_map)
        self.exp_masked.mask = self.mask_total
        self.exp_masked_compressed = self.exp_masked.compressed()
        self.exp_masked_mean = np.mean(self.exp_masked_compressed)

    def log_to_normal(self,array,is_log):
        array_2 = []
        for i in range(len(array)):
            if is_log[i]:
                array_2.append(10**array[i])
            else:
                array_2.append(array[i])
        return array_2

    # def return_non_poiss_samples(self,comp,sample):
    #     b = self.b
    #     self.n_poiss = b.n_poiss
    #     self.model_decompression_non_poiss = np.array(b.model_decompression_key[b.n_poiss:])
    #     model_where = np.where(self.model_decompression_non_poiss[:,0]==comp)[0]+b.n_poiss
    #     is_log_prior = np.array(b.model_decompression_key)[model_where][:,1]
    #     is_log_prior =  list(is_log_prior=='True')
    #     samples_model_not_log = self.log_to_normal(np.array(sample)[model_where],is_log_prior) 

    #     # print self.model_decompression_non_poiss
    #     # print is_log_prior
    #     # print b.model_decompression_key
    #     # print model_where

    #     return samples_model_not_log


    def return_non_poiss_samples(self,comp,sample,ebin):
        b = self.b
        self.n_poiss = b.n_poiss
        self.model_decompression_non_poiss = np.array(b.model_decompression_key[b.n_poiss:])
        model_where = np.where(self.model_decompression_non_poiss[:,0]==comp)[0]+b.n_poiss
        is_log_prior = np.array(b.model_decompression_key)[model_where][:,1]
        is_log_prior =  list(is_log_prior=='True')
        samples_model_not_log = self.log_to_normal(np.array(sample)[model_where],is_log_prior)

        samples_model_not_log_fixed = samples_model_not_log[0:3] + [samples_model_not_log[3+ebin]] 
        samples_model_not_log_fixed[0] = samples_model_not_log_fixed[0]/samples_model_not_log_fixed[-1]

        return samples_model_not_log_fixed


    # def return_poiss_samples(self,comp,sample):
    #     b = self.b
    #     self.n_poiss = b.n_poiss
    #     self.model_decompression_poiss = np.array(b.model_decompression_key[0:b.n_poiss])
    #     model_where = np.where(self.model_decompression_poiss[:,0]==comp)[0]
    #     is_log_prior = np.array(b.model_decompression_key)[model_where][:,1]
    #     is_log_prior =  list(is_log_prior=='True')
    #     samples_model_not_log = self.log_to_normal(np.array(sample)[model_where],is_log_prior)[0]

    #     return samples_model_not_log



    def return_poiss_samples(self,comp,sample):
        b = self.b
        self.n_poiss = b.n_poiss
        self.model_decompression_poiss = np.array(b.model_decompression_key[0:b.n_poiss])
        model_where = np.where(self.model_decompression_poiss[:,0]==comp)[0]
        is_log_prior = self.model_decompression_poiss[model_where][0][1]
        samples_model_not_log = self.log_to_normal(np.array(sample)[model_where],is_log_prior)[0]

        return samples_model_not_log




    # def dnds(self,comp,sample,s): #over full region
    #     # print 'new sample!'
    #     # print sample
    #     samples_reduced = self.return_non_poiss_samples(comp,sample)
    #     # print samples_reduced

    #     #print start_pos,end_pos,samples_reduced_log,samples_reduced

    #     APS, n1,n2,Sb = samples_reduced
    #     if s < Sb:
    #         dnds = APS*(s/Sb)**(-n2)
    #     else:
    #         dnds = APS*(s/Sb)**(-n1)
    #     # print APS, n1, n2, Sb, dnds
    #     return dnds



    def dnds(self,comp,sample,s,ebin=0): #over full region
        samples_reduced = self.return_non_poiss_samples(comp,sample,ebin)
        APS, n1,n2,Sb = samples_reduced
        dnds = np.zeros(len(s))
        where_vec = np.where(s < Sb)[0]
        not_where_vec = np.where(s >= Sb)[0]
        dnds[where_vec] = APS*(s[where_vec]/Sb)**(-n2)
        dnds[not_where_vec] = APS*(s[not_where_vec]/Sb)**(-n1)
        return dnds


    def mask_compress_data(self):
        b= self.b
        self.CTB_masked = [hp.ma(the_data) for the_data in b.count_map]
        for the_map in self.CTB_masked:
            the_map.mask = self.mask_total
        self.CTB_masked_compressed = [the_map.compressed() for the_map in self.CTB_masked]

        self.CTB_masked_compressed_summed = np.sum(self.CTB_masked_compressed,axis=0)
        self.total_counts = float(np.sum(self.CTB_masked_compressed_summed))
        print 'The total number of counts in this region is ', self.total_counts

    def calculate_flux_arrays_non_poiss(self,comp,smin=0,smax=1000,nsteps=10000):
        b = self.b
        self.mask_compress_exposure_map()
        self.mask_compress_data()

        #deal with the template
        summed_template = b.templates_dict_nested[comp]['summed_templates_not_compressed']
        #mean_exposure = np.mean(self.the_dict['exposure_map'],axis=0)
        summed_template_masked_compressed = self.mask_and_compress(summed_template,self.mask_total)
        self.template_sum = np.sum(summed_template_masked_compressed)

        #source count
        # self.sarray = np.arange(smin,smax,float((smax-smin)/float(nsteps)))
        self.sarray = 10**np.linspace( np.log10(smin),np.log10(smax),nsteps)
        self.ds_array = [self.sarray[i+1] - self.sarray[i] for i in range(nsteps-1)]
        self.ds_array = np.array( self.ds_array + [ self.ds_array[-1] ])
        self.ds = self.sarray[1]-self.sarray[0]
        #self.dnds_array = np.array([self.dnds(comp_number, sam,s) for sam in self.samples])
        # self.flux_array_non_poiss = map( lambda sam : np.sum(map(lambda s: self.template_sum*self.dnds(comp, sam,s)*s*self.ds, self.sarray) ) , b.samples )

        self.flux_array_non_poiss = np.sum( np.array([self.template_sum*self.dnds(comp, sam,self.sarray)*self.sarray*self.ds_array for sam in b.samples]),axis=1 )


        self.flux_fraction_array_non_poiss = np.array(self.flux_array_non_poiss) / self.total_counts

    def calculate_flux_arrays_poiss(self,comp):
        b = self.b
        self.mask_compress_exposure_map()
        self.mask_compress_data()

        #deal with the template
        summed_template = b.templates_dict_nested[comp]['summed_templates_not_compressed']
        #mean_exposure = np.mean(self.the_dict['exposure_map'],axis=0)
        summed_template_masked_compressed = self.mask_and_compress(summed_template,self.mask_total)
        self.template_sum = np.sum(summed_template_masked_compressed)

        self.flux_array_poiss = map( lambda sam : self.template_sum*self.return_poiss_samples(comp,sam) , b.samples )


        self.flux_fraction_array_poiss = np.array(self.flux_array_poiss) / self.total_counts
        

    def plot_flux_fraction_non_poiss(self,comp,smin=0.00001,smax=1000,nsteps=1000,bins=200,bin_range=(0,100),color_vert = 'b',ls_vert = 'dashed',*args,**kwargs):
        self.calculate_flux_arrays_non_poiss(comp,smin=smin,smax=smax,nsteps=nsteps)
        frac_hist_comp, bin_edges_comp = np.histogram(100*np.array(self.flux_fraction_array_non_poiss),bins=bins,range=bin_range) #bins were 600 #bins=200, range=(0,100)
        qvalues_comp = triangle.quantile(100*np.array(self.flux_fraction_array_non_poiss), [0.16, 0.5, 0.84])
        plt.plot(bin_edges_comp[:-1], frac_hist_comp/float(sum(frac_hist_comp)),*args,**kwargs ) #'b', drawstyle='steps-mid',label = 'NFW PS'
        for q in qvalues_comp:
            plt.axvline(q, ls=ls_vert, color=color_vert)
        self.qvalues_comp = qvalues_comp

    def plot_flux_fraction_poiss(self,comp,bins=200,bin_range=(0,100),color_vert = 'b',ls_vert = 'dashed',*args,**kwargs):
        self.calculate_flux_arrays_poiss(comp)
        frac_hist_comp, bin_edges_comp = np.histogram(100*np.array(self.flux_fraction_array_poiss),bins=bins,range=bin_range) #bins were 600 #bins=200, range=(0,100)
        qvalues_comp = triangle.quantile(100*np.array(self.flux_fraction_array_poiss), [0.16, 0.5, 0.84])
        plt.plot(bin_edges_comp[:-1], frac_hist_comp/float(sum(frac_hist_comp)),*args,**kwargs ) #'b', drawstyle='steps-mid',label = 'NFW PS'
        for q in qvalues_comp:
            plt.axvline(q, ls=ls_vert, color=color_vert)
        self.qvalues_comp = qvalues_comp




