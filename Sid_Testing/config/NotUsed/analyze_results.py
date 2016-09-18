import os

# current_dir = os.getcwd()
# change_path = ".."
# os.chdir(change_path)

import healpy as hp
import numpy as np
import config.set_dirs as sd
import config.make_templates as mt
import config.compute_PSF as CPSF

import pulsars.masks as masks

import pymultinest
import triangle
import matplotlib.pyplot as plt

import sys
import inspect

import logging


class analyze_results:
    def __init__(self,**kwargs):
        if "dict_dir" in kwargs:
            self.dict_dir = kwargs['dict_dir']
            self.run_tag = kwargs['run_tag']

            #self.the_dict = np.load(self.dict_dir+ self.run_tag + '.npy')

            self.keys = np.load(self.dict_dir+ self.run_tag + '-keys' + '.npz')['arr_0']
            self.values = np.load(self.dict_dir+ self.run_tag + '-values' + '.npz')['arr_0']

            self.the_dict = {self.keys[i]: self.values[i] for i in range(np.shape(self.keys)[0])}
        else:
            self.the_dict = kwargs['the_dict'] 

        self.samples = self.the_dict['samples']
        self.nside = self.the_dict['nside']
        self.npix = self.the_dict['npix']

        #Lina: Running the mask a. set total mask
        #set_mask_total()

    def delete_dict(self):
        os.remove(self.dict_dir+ self.run_tag + '-keys' + '.npz')
        os.remove(self.dict_dir+ self.run_tag + '-values' + '.npz')
        print 'deleted the dictionary corresponding to the file ' + self.dict_dir+ self.run_tag + '-keys' + '.npz' 

    def make_triangle(self,plot_name='triangle.pdf',return_fig = False,save_fig = True,**kwargs):
        if 'plot_dir_for_triangle' in kwargs.keys():
            plot_dir_for_triangle = kwargs['plot_dir_for_triangle']
        else:
            plot_dir_for_triangle = self.the_dict['plots_dir_for_run']

        self.the_range = [0.95 for i in range(self.the_dict['n_params'])]
        self.tri_fig = triangle.corner(self.the_dict['samples'], labels=self.the_dict['params'],smooth=1.5,smooth1d=1,
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

    def counts_to_intensity(self,the_map):
        return the_map / self.the_dict['CTB_exposure_maps'] / hp.nside2pixarea(self.nside)

    def mask_and_compress(self,the_map,mask):
        map_masked = hp.ma(the_map)
        map_masked.mask = mask
        return map_masked.compressed()

    def return_template_intensity(self,mask,comp):
        summed_template = self.the_dict['templates_dict_nested'][comp]['summed_templates_not_compressed']
        mean_exposure = np.mean(self.the_dict['CTB_exposure_maps'],axis=0)

        summed_template_masked_compressed = self.mask_and_compress(summed_template,mask)
        mean_exposure_masked_compressed = self.mask_and_compress(mean_exposure,mask)

        template_intensity = np.mean(summed_template_masked_compressed/mean_exposure_masked_compressed/hp.nside2pixarea(self.nside))

        return template_intensity
      #   template_array = self.the_dict['templates'][comp_number]
      #   template_intensity_array_um = self.counts_to_intensity(template_array)
      # #  mask_array = [mask for i in range(np.shape(template_array)[0])]

      #   template_intensity_array_masked_compressed = np.array( [self.mask_and_compress(tf_um,mask) for tf_um in template_intensity_array_um ] )
      #   len_of_array = len(template_intensity_array_masked_compressed[0])

      #   template_intensity_array = np.sum(template_intensity_array_masked_compressed,axis=1)/len_of_array

      #   return np.sum(template_intensity_array,axis=0)

    def calculate_poiss_flux(self,mask,comp,log_prior=False):
        template_intensity = self.return_template_intensity(mask,comp)
        if not log_prior:
            scan_intensity_array = template_intensity*np.array([self.samples[i][self.the_dict['poiss_comp_numbers'][comp]] for i in range(len(self.samples))])
        else:
            scan_intensity_array = template_intensity*np.array([10**self.samples[i][self.the_dict['poiss_comp_numbers'][comp]] for i in range(len(self.samples))])
        return scan_intensity_array

    def calculate_number_non_ps(self):
        self.n_single_ps = len(np.where(np.array(self.the_dict['comp_array'])=='ps')[0])
        self.n_non_ps = self.the_dict['n_poiss'] - self.n_single_ps
        #print 'number of non_ps is ', self.n_non_ps

    def calculate_single_ps_flux(self,ps_number,log_prior=False):
        self.calculate_number_non_ps()
        if not log_prior:
            scan_intensity_array = np.array([self.samples[i][ps_number + self.n_non_ps] for i in range(len(self.samples))])
        else:
            scan_intensity_array = np.array([10**self.samples[i][ps_number + self.n_non_ps] for i in range(len(self.samples))])
        return scan_intensity_array

    def calculate_ps_norm_factor(self,sample,other_comp_numbers):
        n1, n2, Sb = [sample[i] for i in other_comp_numbers]
        if self.the_dict['Sc'] != 'False':
            max_photon_cut = self.the_dict['Sc']
            norm_factor = (Sb**2)*( 1/(2.-n2) - (1 - (Sb**(n1-2))*(max_photon_cut**(2-n1)))/(2.-n1))
        else:
            norm_factor = (n2 - n1) * Sb**2 / ((n1 - 2.) * (n2 - 2))
        return norm_factor

    def log_to_normal(self,array,is_log):
        array_2 = []
        for i in range(len(array)):
            if is_log[i]:
                array_2.append(10**array[i])
            else:
                array_2.append(array[i])
        return array_2

    def calculate_ps_flux(self,mask,comp_number):

        self.n_poiss = self.the_dict['n_poiss']
        start_pos = self.n_poiss+np.sum(self.the_dict['non_poiss_list_num_model_params'][0:comp_number-1])
        end_pos = start_pos + self.the_dict['non_poiss_list_num_model_params'][comp_number]
        samples_reduced_log = [sam[start_pos:end_pos] for sam in self.samples]

        samples_reduced = [self.log_to_normal(sam,self.the_dict['non_poiss_list_is_log_prior_uncompressed'][comp_number] ) for sam in samples_reduced_log] 

        template_intensity = self.return_template_intensity(mask,self.the_dict['non_poiss_list_template_number'][comp_number])

        scan_intensity_array = template_intensity*np.array([self.calculate_ps_norm_factor(sam,[1,2,3])*sam[0] for sam in samples_reduced])
        
        return scan_intensity_array

    def return_poiss_intensity_conf_int(self,comp,quant = [0.16,0.5,0.84]): #intensity, dE, units [photons/cm^2/s/sr]
        mask = self.mask_total
        scan_intensity_array = self.calculate_poiss_flux(mask,comp,log_prior = self.the_dict['poiss_models'][comp]['log_prior'])
       
        return triangle.quantile(scan_intensity_array,quant)

    def return_single_ps_intensity_conf_int(self,ps_number,quant = [0.16,0.5,0.84]): #intensity, dE, units [photons/cm^2/s/sr]
        scan_intensity_array = self.calculate_single_ps_flux(ps_number,log_prior = self.the_dict['poiss_list_is_log_prior'][-1]) #all ps's are same, so just use last
       
        return triangle.quantile(scan_intensity_array,quant)

    def set_energies(self):
        self.En_min = self.the_dict['CTB_en_bins'][0]
        self.En_max = self.the_dict['CTB_en_bins'][-1]
        self.En_center = 10**((np.log10(self.En_max)+np.log10(self.En_min))/2)
        self.dE = self.En_max - self.En_min


    def return_spectrum_conf_int(self,comp,quant  = [0.16,0.5,0.84]):
        self.set_energies()
        self.spectrum = self.En_center**2*np.array(self.return_poiss_intensity_conf_int(comp, quant = quant))/self.dE
        return self.spectrum

    def return_single_ps_spectrum_conf_int(self,ps_number,quant  = [0.16,0.5,0.84]):
        self.set_energies()
        self.spectrum = self.En_center**2*np.array(self.return_single_ps_intensity_conf_int(ps_number, quant = quant))/self.dE
        return self.spectrum

    def return_single_ps_median_counts(self,ps_number):
        return self.return_single_ps_intensity_conf_int(ps_number, quant = [0.5])[0]


    def return_poiss_medians(self):
        self.poiss_medians = [self.the_dict['s']['marginals'][i]['median'] for i in range(self.the_dict['n_poiss']) ]
        return self.poiss_medians

    def return_non_poiss_intensity_conf_int(self,comp_number,quant = [0.16,0.5,0.84]): #intensity, dE, units [photons/cm^2/s/sr]
        mask = self.mask_total
        scan_intensity_array = self.calculate_ps_flux(mask,comp_number)
       
        return triangle.quantile(scan_intensity_array,quant)

    def return_non_poiss_spectrum_conf_int(self,comp_number,quant  = [0.16,0.5,0.84]):
        self.set_energies()
        self.spectrum = self.En_center**2*np.array(self.return_non_poiss_intensity_conf_int(comp_number, quant = quant))/self.dE
        return self.spectrum

    def make_new_ps_model(self,ps_model_dir,model_tag):
        self.calculate_number_non_ps()
        self.new_ps_model = np.sum([self.return_single_ps_median_counts(ps_number)*self.the_dict['summed_templates_not_compressed'][ps_number + self.n_non_ps] for ps_number in range(self.n_single_ps)],axis=0)

        if self.the_dict['use_fixed_template']:
            self.new_plus_old_ps_model = self.new_ps_model + self.the_dict['fixed_templates'][0]
            np.save(ps_model_dir + model_tag+ '-just_bin',self.new_ps_model )
            np.save(ps_model_dir + model_tag,self.new_plus_old_ps_model )
        else:
            np.save(ps_model_dir + model_tag,self.new_ps_model)


    def make_band_mask(self,band_mask_range):
        band_mask = masks.mask_lat_band(90+band_mask_range[0],90+ band_mask_range[1], self.nside)
        return band_mask

    def make_ring_mask(self,inner,outer,lat,lng):
        ring_mask = np.logical_not(masks.mask_ring(inner, outer, lat, lng, self.nside))
        return ring_mask


    def set_mask_total(self,band_mask_range = [-30,30],mask_ring = False,lat= 90,lng = 0,inner = 0, outer= 30,input_mask = False, the_input_mask = 'Null'):
        band_mask = self.make_band_mask(band_mask_range)

        if input_mask:
            self.mask_total = the_input_mask
        else:
            if not mask_ring:
                mask_geom_total = band_mask #fix this!
            else:
                ring_mask = self.make_ring_mask(inner,outer,lat,lng)
                mask_geom_total = ring_mask + band_mask

            self.mask_total = mask_geom_total

    def load_3FGL_data(self):
        if self.only_gal:
            psc_tag = 'sphCoordsGAL.dat'
            psc_flux_tag = 'fluxTabGALBin.dat'
        else:
            psc_tag = 'sphCoords.dat'
            psc_flux_tag = 'fluxTabBin.dat'

            psc_3FGL_name = self.the_dict['ps_flux_dir'] + psc_tag
            psc_3FGL_flux_name = self.the_dict['ps_flux_dir'] + psc_flux_tag

            self.psc_3FGL = np.loadtxt(psc_3FGL_name)
            self.psc_3FGL_flux = np.loadtxt(psc_3FGL_flux_name)

            self.psc_3FGL_flux_bin = np.sum(self.psc_3FGL_flux[::,self.the_dict['CTB_en_min']:self.the_dict['CTB_en_max']],axis=1)

    def make_3FGL_predicted_counts_map(self):
        self.psArray = np.zeros(self.npix)
        self.pix_num_array = [hp.ang2pix(self.nside, self.psc_3FGL[i,0], self.psc_3FGL[i,1]) for i in range(np.shape(self.psc_3FGL)[0]) ]
        self.counts_vec = np.array(self.psc_3FGL_flux_bin)*self.the_dict['total_exposure_map'][np.array(self.pix_num_array)]
        self.psArray[self.pix_num_array] = self.counts_vec

    def mask_3FGL_predicted_counts_map(self):
        self.ps_array_masked = hp.ma(self.psArray)
        self.ps_array_masked.mask = self.mask_total
        self.ps_array_compressed = self.ps_array_masked.compressed()
        self.ps_array_non_zero = self.ps_array_compressed[np.where(self.ps_array_compressed)] 

        #self.area_mask = len(self.mask_total)*hp.nside2pixarea(self.nside)*(360/(2.*np.pi))**2 #in deg^2

    def mask_compress_exposure_map(self):
        self.area_mask = len(self.mask_total)*hp.nside2pixarea(self.nside)*(360/(2.*np.pi))**2 #in deg^2
        self.exp_masked = hp.ma(self.the_dict['total_exposure_map'])
        self.exp_masked.mask = self.mask_total
        self.exp_masked_compressed = self.exp_masked.compressed()
        self.exp_masked_mean = np.mean(self.exp_masked_compressed)

    def make_3FGL_predicted_flux_map(self):
        self.mask_compress_exposure_map()
        self.ps_flux_masked_compressed = self.ps_array_non_zero/self.exp_masked_compressed[np.where(self.ps_array_compressed)]

    def make_3FGL_flux_histogram(self):
        yvals, bins = np.histogram(self.ps_flux_masked_compressed,bins=10**np.linspace(np.log10(self.flux_min), np.log10(self.flux_max), self.n_flux_bins));
        bin_centers = np.array([10**((np.log10(bins[i])+np.log10(bins[i+1]))/2) for i in range(np.shape(bins)[0]-1)])
        bin_width = np.array([bins[i+1]-bins[i] for i in range(np.size(bins)-1)])

        rescale = 1/self.area_mask
        self.yvals_rescale = yvals*rescale/bin_width
        self.bin_centers_rescale = bin_centers

        if self.error_range == 0.68:
            errors_sigma = np.loadtxt(self.the_dict['ps_flux_dir'] + 'ErrorTab68.dat')
        else:
            errors_sigma = np.loadtxt(self.the_dict['ps_flux_dir'] + 'ErrorTab95.dat')
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

        
    def configure_3FGL(self,flux_min=10**-11,flux_max=10**-8,n_flux_bins = 25,only_gal = False,error_range = 0.68):
        self.only_gal = only_gal
        self.flux_min = flux_min
        self.flux_max = flux_max
        self.n_flux_bins = n_flux_bins
        self.error_range = error_range


        self.load_3FGL_data()
        self.make_3FGL_predicted_counts_map()
        self.mask_3FGL_predicted_counts_map()
        self.make_3FGL_predicted_flux_map()
        self.make_3FGL_flux_histogram()


    def plot_3FGL(self,*args,**kwargs):
        plt.errorbar(self.bin_centers_rescale,self.yvals_rescale,xerr=[self.x_errors_L,self.x_errors_H],yerr=[self.error_L,self.error_H],*args,**kwargs)
        #fmt=fmt,color=color,markersize=markersize,label = label
        #fmt = 'o', color='black',markersize=7,label='3FGL PS'

    def dnds(self,comp_number,sample,s): #over full region

        self.n_poiss = self.the_dict['n_poiss']
        start_pos = self.n_poiss+np.sum(self.the_dict['non_poiss_list_num_model_params'][0:comp_number-1])
        end_pos = start_pos + self.the_dict['non_poiss_list_num_model_params'][comp_number]
        samples_reduced_log = sample[start_pos:end_pos] 
        samples_reduced = self.log_to_normal(samples_reduced_log,self.the_dict['non_poiss_list_is_log_prior_uncompressed'][comp_number] ) 

        #print start_pos,end_pos,samples_reduced_log,samples_reduced

        APS, n1,n2,Sb = samples_reduced
        if s < Sb:
            dnds = APS*(s/Sb)**(-n2)
        else:
            dnds = APS*(s/Sb)**(-n1)
        return dnds

    def calculate_dnds_arrays(self,comp_number,smin=0,smax=1000,nsteps=10000,qs =[0.16, 0.5, 0.84] ):
        self.mask_compress_exposure_map()

        #deal with the template
        summed_template = self.the_dict['summed_templates_not_compressed'][self.the_dict['non_poiss_list_template_number'][comp_number]]
        mean_exposure = np.mean(self.the_dict['CTB_exposure_maps'],axis=0)
        summed_template_masked_compressed = self.mask_and_compress(summed_template,self.mask_total)
        self.template_sum = np.sum(summed_template_masked_compressed)

        #source count
        self.sarray = np.arange(smin,smax,float((smax-smin)/float(nsteps)))
        #self.dnds_array = np.array([self.dnds(comp_number, sam,s) for sam in self.samples])
        qArray = [triangle.quantile(np.array([self.dnds(comp_number, sam,s) for sam in self.samples]),qs) for s in self.sarray]

        rf = self.template_sum*self.exp_masked_mean/self.area_mask
        self.qmean = rf*np.array([np.mean(np.array([self.dnds(comp_number, sam,s) for sam in self.samples])) for s in self.sarray])
        self.q16 = rf*np.array([q[0] for q in qArray]) #units:photons^-1 cm^2 s deg^-2
        self.q5 = rf*np.array([q[1] for q in qArray])
        self.q84 = rf*np.array([q[2] for q in qArray])

        self.flux_array = self.sarray/self.exp_masked_mean
        
        # rf = self.exp_masked_mean/self.area_mask #units:photons^-1 cm^2 s deg^-2
        # return sarray/self.exp_masked_mean, q16*rf, q5*rf, q84*rf, qmean*rf

    # class plot_source_count:
    #     def __init__(self,comp_number,smin=0,smax=1000,nsteps=10000,qs=[0.16,0.5,0.84],calculate_dnds=False,*args,**kwargs):
    #         if calculate_dnds:
    #             self.calculate_dnds_arrays(comp_number,smin=smin,smax=smax,nsteps=nsteps,qs=qs)
    #         plt.fill_between(self.flux_array,self.q16,self.q84,*args,**kwargs)
    #         #plt.plot(self.flux_array ,self.q5,'k--',label="NFW PS")


    def plot_source_count_band(self,comp_number,smin=0,smax=1000,nsteps=10000,qs=[0.16,0.5,0.84],calculate_dnds=False,*args,**kwargs):
        if calculate_dnds:
            self.calculate_dnds_arrays(comp_number,smin=smin,smax=smax,nsteps=nsteps,qs=qs)

        plt.fill_between(self.flux_array,self.q16,self.q84,*args,**kwargs)#facecolor='Chartreuse',interpolate=True, alpha=1,linewidth=0)
        #plt.plot(self.flux_array ,self.q5,'k--',label="NFW PS")

    def plot_source_count_median(self,comp_number,smin=0,smax=1000,nsteps=10000,qs=[0.16,0.5,0.84],calculate_dnds=False,*args,**kwargs):
        if calculate_dnds:
            self.calculate_dnds_arrays(comp_number,smin=smin,smax=smax,nsteps=nsteps,qs=qs)
        plt.plot(self.flux_array ,self.q5,*args,**kwargs)

    def plot_source_count_mean(self,comp_number,smin=0,smax=1000,nsteps=10000,qs=[0.16,0.5,0.84],calculate_dnds=False,*args,**kwargs):
        if calculate_dnds:
            self.calculate_dnds_arrays(comp_number,smin=smin,smax=smax,nsteps=nsteps,qs=qs)
        plt.plot(self.flux_array ,self.qmean,*args,**kwargs)








