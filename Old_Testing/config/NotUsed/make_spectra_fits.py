"""
Local version of make_spectra_plot.py that plots functional fit, 
fit parameters and 3FGL source spectra
"""

import sys, os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import healpy as hp
from scipy.optimize import curve_fit
#from sympy.functions.special.delta_functions import Heaviside
# from astropy.modeling import powerlaws, fitting
import fermi.fermi_plugin as fp
import pulsars.CTB as CTB

class make_spectra_plot:
    def __init__(self,spect_file):
        self.spect_file = spect_file
        self.spect = np.load(self.spect_file)
        self.extract_energies()
        self.extract_keys()
        self.make_spectra_arrays()

    def extract_energies(self):
        self.energies = np.array( [ sp[0] for sp in self.spect] )

    def extract_keys(self):
        self.keys = np.array( self.spect[0][1].keys() )

    def make_spectra_comp(self,comp):
        spectra_array_lower = []
        spectra_array_mid = []
        spectra_array_upper = []
        for sp in self.spect:
            spectrum_lower, spectrum_mid, spectrum_upper = list(sp[1][comp])
            spectra_array_lower.append(spectrum_lower)
            spectra_array_mid.append(spectrum_mid)
            spectra_array_upper.append(spectrum_upper)
        return spectra_array_lower, spectra_array_mid, spectra_array_upper

    def make_spectra_arrays(self):
        self.total_spectra_array_lower = []
        self.total_spectra_array_mid = []
        self.total_spectra_array_upper = []
        for comp in self.keys:
            comps = self.make_spectra_comp(comp)
            self.total_spectra_array_lower += [comps[0]]
            self.total_spectra_array_mid += [comps[1]]
            self.total_spectra_array_upper += [comps[2]]

    def plot_spectra_band(self,model,*args,**kwargs):
        
        from scipy.interpolate import spline

        model_number=np.where(self.keys==model)[0][0]
        enew = np.linspace(np.log10(self.energies.min()),np.log10(self.energies.max()),12)
        lower_smooth = spline(np.log10(self.energies),np.log10(self.total_spectra_array_lower[model_number]),enew)        
        upper_smooth = spline(np.log10(self.energies),np.log10(self.total_spectra_array_upper[model_number]),enew)        

        plt.fill_between(self.energies, self.total_spectra_array_lower[model_number], self.total_spectra_array_upper[model_number],*args,**kwargs)

    def plot_spectra_median(self,model,*args,**kwargs):
        from scipy.interpolate import spline
        model_number=np.where(self.keys==model)[0][0]
        enew = np.linspace(np.log10(self.energies.min()),np.log10(self.energies.max()),5)
        medians_smooth = spline(np.log10(self.energies),np.log10(self.total_spectra_array_mid[model_number]),enew)        
        plt.plot(self.energies, self.total_spectra_array_mid[model_number], *args,**kwargs)

    def plot_spectra_bars(self,model,scaling=1,*args,**kwargs):

        model_number=np.where(self.keys==model)[0][0]
        # print model_number, model
        # print self.total_spectra_array_mid[model_number]
        print model, 'spectra are'
        print 'energies', list(self.energies)
        print 'lower', self.total_spectra_array_lower[model_number]
        print 'mid', self.total_spectra_array_mid[model_number]
        print 'upper', self.total_spectra_array_upper[model_number]
        plt.errorbar(self.energies, scaling*np.array(self.total_spectra_array_mid[model_number]), yerr=[np.array(scaling*np.array(self.total_spectra_array_mid[model_number]))-scaling*np.array(self.total_spectra_array_lower[model_number]), scaling*np.array(self.total_spectra_array_upper[model_number])-scaling*np.array(self.total_spectra_array_mid[model_number])], markersize=3, *args,**kwargs)    

    def plot_spectra_err(self,model,scaling=1,*args,**kwargs):

        model_number=np.where(self.keys==model)[0][0]
        # print model_number, model
        # print self.total_spectra_array_mid[model_number]
        print model, 'spectra are'
        print 'energies', list(self.energies)
        print 'lower', self.total_spectra_array_lower[model_number]
        print 'mid', self.total_spectra_array_mid[model_number]
        print 'upper', self.total_spectra_array_upper[model_number]
        if len(self.energies) > 1:
            xerr1 = (self.energies-[10**((np.log10(self.energies[i]*self.energies[0]/self.energies[1])+np.log10(self.energies[i]))/2) for i in range(len(self.energies))])
            xerr2 = ([10**((np.log10(self.energies[i]*self.energies[1]/self.energies[0])+np.log10(self.energies[i]))/2) for i in range(len(self.energies))])-self.energies
            print xerr2
            plt.errorbar(self.energies, scaling*np.array(self.total_spectra_array_mid[model_number]), yerr=[np.array(scaling*np.array(self.total_spectra_array_mid[model_number]))-scaling*np.array(self.total_spectra_array_lower[model_number]), scaling*np.array(self.total_spectra_array_upper[model_number])-scaling*np.array(self.total_spectra_array_mid[model_number])], xerr=[xerr1,xerr2], markersize=3, *args,**kwargs)     
        else:
            plt.errorbar(self.energies, scaling*np.array(self.total_spectra_array_mid[model_number]), yerr=[np.array(scaling*np.array(self.total_spectra_array_mid[model_number]))-scaling*np.array(self.total_spectra_array_lower[model_number]), scaling*np.array(self.total_spectra_array_upper[model_number])-scaling*np.array(self.total_spectra_array_mid[model_number])], markersize=3, *args,**kwargs) 
        # plt.plot(self.energies, scaling*np.array(self.total_spectra_array_mid[model_number]),marker='o',linestyle='none',markersize=6,color='darkgrey')

    def return_spectra_err(self,model,scaling=1,*args,**kwargs):

        model_number=np.where(self.keys==model)[0][0]
        # print model_number, model
        # print self.total_spectra_array_mid[model_number]
        print model, 'spectra are'
        print 'energies', list(self.energies)
        print 'lower', self.total_spectra_array_lower[model_number]
        print 'mid', self.total_spectra_array_mid[model_number]
        print 'upper', self.total_spectra_array_upper[model_number]
        xerr1 = (self.energies-[10**((np.log10(self.energies[i]*self.energies[0]/self.energies[1])+np.log10(self.energies[i]))/2) for i in range(len(self.energies))])
        xerr2 = ([10**((np.log10(self.energies[i]*self.energies[1]/self.energies[0])+np.log10(self.energies[i]))/2) for i in range(len(self.energies))])-self.energies
        return self.total_spectra_array_mid[model_number], self.energies, xerr1, xerr2

    def plot_spectra_err_2models(self,model1,model2,scaling=1,*args,**kwargs):

        model_number1=np.where(self.keys==model1)[0][0]
        model_number2=np.where(self.keys==model2)[0][0]

        total_spectra_array_lower = [m1+m2 for m1,m2 in zip(self.total_spectra_array_lower[model_number1],self.total_spectra_array_lower[model_number2])]
        total_spectra_array_mid = [m1+m2 for m1,m2 in zip(self.total_spectra_array_mid[model_number1],self.total_spectra_array_mid[model_number2])]
        total_spectra_array_upper = [m1+m2 for m1,m2 in zip(self.total_spectra_array_upper[model_number1],self.total_spectra_array_upper[model_number2])]


        # print model_number, model
        # print self.total_spectra_array_mid[model_number]
        # plt.plot(self.energies, scaling*np.array(sum_temp),*args,**kwargs)   

        # model_number=np.where(self.keys==model)[0][0]
        # print model_number, model
        # print self.total_spectra_array_mid[model_number]
        # print model, 'spectra are'
        # print 'energies', list(self.energies)
        # print 'lower', self.total_spectra_array_lower[model_number]
        # print 'mid', self.total_spectra_array_mid[model_number]
        # print 'upper', self.total_spectra_array_upper[model_number]
        if len(self.energies) > 1 :
            xerr1 = (self.energies-[10**((np.log10(self.energies[i]*self.energies[0]/self.energies[1])+np.log10(self.energies[i]))/2) for i in range(len(self.energies))])
            xerr2 = ([10**((np.log10(self.energies[i]*self.energies[1]/self.energies[0])+np.log10(self.energies[i]))/2) for i in range(len(self.energies))])-self.energies
            print xerr2
            plt.errorbar(self.energies, scaling*np.array(total_spectra_array_mid), yerr=[np.array(scaling*np.array(total_spectra_array_mid))-scaling*np.array(total_spectra_array_lower), scaling*np.array(total_spectra_array_upper)-scaling*np.array(total_spectra_array_mid)], xerr=[xerr1,xerr2], markersize=3, *args,**kwargs)     
        else:
            plt.errorbar(self.energies, scaling*np.array(total_spectra_array_mid), yerr=[np.array(scaling*np.array(total_spectra_array_mid))-scaling*np.array(total_spectra_array_lower), scaling*np.array(total_spectra_array_upper)-scaling*np.array(total_spectra_array_mid)], markersize=3, *args,**kwargs)  
        # plt.plot(self.energies, scaling*np.array(total_spectra_array_mid),marker='o',linestyle='none',markersize=6,color='black')

    def return_spectra_err_2models(self,model1,model2,scaling=1,*args,**kwargs):

        model_number1=np.where(self.keys==model1)[0][0]
        model_number2=np.where(self.keys==model2)[0][0]

        total_spectra_array_lower = [m1+m2 for m1,m2 in zip(self.total_spectra_array_lower[model_number1],self.total_spectra_array_lower[model_number2])]
        total_spectra_array_mid = [m1+m2 for m1,m2 in zip(self.total_spectra_array_mid[model_number1],self.total_spectra_array_mid[model_number2])]
        total_spectra_array_upper = [m1+m2 for m1,m2 in zip(self.total_spectra_array_upper[model_number1],self.total_spectra_array_upper[model_number2])]


        # print model_number, model
        # print self.total_spectra_array_mid[model_number]
        # plt.plot(self.energies, scaling*np.array(sum_temp),*args,**kwargs)   

        # model_number=np.where(self.keys==model)[0][0]
        # print model_number, model
        # print self.total_spectra_array_mid[model_number]
        # print model, 'spectra are'
        # print 'energies', list(self.energies)
        # print 'lower', self.total_spectra_array_lower[model_number]
        # print 'mid', self.total_spectra_array_mid[model_number]
        # print 'upper', self.total_spectra_array_upper[model_number]
        xerr1 = (self.energies-[10**((np.log10(self.energies[i]*self.energies[0]/self.energies[1])+np.log10(self.energies[i]))/2) for i in range(len(self.energies))])
        xerr2 = ([10**((np.log10(self.energies[i]*self.energies[1]/self.energies[0])+np.log10(self.energies[i]))/2) for i in range(len(self.energies))])-self.energies
        return total_spectra_array_mid, self.energies, xerr1, xerr2

    def plot_spectra_bars_2models(self,model1,model2,scaling=1,*args,**kwargs):

        model_number1=np.where(self.keys==model1)[0][0]
        model_number2=np.where(self.keys==model2)[0][0]

        total_spectra_array_lower = [m1+m2 for m1,m2 in zip(self.total_spectra_array_lower[model_number1],self.total_spectra_array_lower[model_number2])]
        total_spectra_array_mid = [m1+m2 for m1,m2 in zip(self.total_spectra_array_mid[model_number1],self.total_spectra_array_mid[model_number2])]
        total_spectra_array_upper = [m1+m2 for m1,m2 in zip(self.total_spectra_array_upper[model_number1],self.total_spectra_array_upper[model_number2])]

        # print model_number, model
        # print self.total_spectra_array_mid[model_number]
        plt.errorbar(self.energies, scaling*np.array(total_spectra_array_mid), yerr=[np.array(scaling*np.array(total_spectra_array_mid))-scaling*np.array(total_spectra_array_lower), scaling*np.array(total_spectra_array_upper)-scaling*np.array(total_spectra_array_mid)], markersize=3, *args,**kwargs)    

        # plt.plot(self.energies, scaling*np.array(sum_temp),*args,**kwargs)    

    def get_spectra_median(self,model,*args,**kwargs):
        model_number=np.where(self.keys==model)[0][0]
        return self.energies, self.total_spectra_array_mid[model_number]

class extract_fit_coeffs:
    def __init__(self,spect_files,comp='iso',fitfunc='PLE'):
        
        self.spect_files = spect_files
        self.spects = np.array([spect_file for spect_file in self.spect_files])
        self.comp = comp
        
        self.extract_energies()
        self.extract_keys()
        self.make_spectra_arrays()

        self.total_spectra_array_lower=np.ravel(self.total_spectra_array_lower)
        self.total_spectra_array_mid=np.ravel(self.total_spectra_array_mid)
        self.total_spectra_array_upper=np.ravel(self.total_spectra_array_upper)

        self.energies = np.array([val for sublist in self.energies for val in sublist])
        self.total_spectra_array_mid = np.array([val for sublist in self.total_spectra_array_mid for val in sublist])

        if fitfunc == 'PLE':
            fitfunc = self.power_law_exponent
            initguess = [0.95*10**-4,2.32,80]
        elif fitfunc == 'PL':
            fitfunc = self.power_law
            initguess = [0.95*10**-4,2.32]
        elif fitfunc == 'BPL':
            fitfunc = self.broken_power_law
            initguess = [3*10**-5,2.32,2.32,20]

        self.E2dNdE = [(E2dNdE) for E2dNdE,E in zip(self.total_spectra_array_mid, self.energies)]

        self.E2dNdE = [x for (y,x) in sorted(zip(self.energies, self.E2dNdE))]
        self.energies = sorted(self.energies)

        self.popt, self.pcov = curve_fit(fitfunc, np.array(self.energies), np.array(self.E2dNdE),maxfev = 100000,p0=initguess)
        
        self.E2dNdEfit = [fitfunc(E,*self.popt) for E in  self.energies]

    def theta(self,x):
        return 0.5 * (np.sign(x) + 1)

    def return_fit_spec(self):
        return self.energies, self.E2dNdEfit

    def extract_energies(self):
        self.energies = np.array( [ [ sp[0] for sp in sps ] for sps in self.spects])

    def extract_keys(self):
        self.keys = [[] for i in range(len(self.spects))]
        for i in range(len(self.spects)):
            self.keys[i] = np.array( self.spects[i][0][1].keys() )

    def make_spectra_comp(self,comp,spect):
        spectra_array_lower = []
        spectra_array_mid = []
        spectra_array_upper = []
        for sp in spect:
            spectrum_lower, spectrum_mid, spectrum_upper = list(sp[1][self.comp])
            spectra_array_lower.append(spectrum_lower)
            spectra_array_mid.append(spectrum_mid)
            spectra_array_upper.append(spectrum_upper)
        return spectra_array_lower, spectra_array_mid, spectra_array_upper

    def make_spectra_arrays(self):
        self.total_spectra_array_lower = [[] for i in range(len(self.spects))]
        self.total_spectra_array_mid = [[] for i in range(len(self.spects))]
        self.total_spectra_array_upper = [[] for i in range(len(self.spects))]
        for i in range(len(self.spects)):
            for comp in self.keys[i]:
                if comp == self.comp:
                    comps = self.make_spectra_comp(comp,self.spects[i])
                    self.total_spectra_array_lower[i] += [comps[0]]
                    self.total_spectra_array_mid[i] += [comps[1]]
                    self.total_spectra_array_upper[i] += [comps[2]]
        self.total_spectra_array_lower = np.array(self.total_spectra_array_lower)
        self.total_spectra_array_mid = np.array(self.total_spectra_array_mid)
        self.total_spectra_array_upper = np.array(self.total_spectra_array_upper)

    def plot_PLE(self, I100,gamma,Ecut,*args,**kwargs):
        EgammaList = np.logspace(-1,3,100)
        print 'actually plotting PLE...'
        dNdEgamma = [self.power_law_exponent(Egamma,I100,gamma,Ecut) for Egamma in EgammaList]
        E2dNdEGamma = [E2dNdE for E,E2dNdE in zip(EgammaList,dNdEgamma)]
        # print 'points are', EgammaList,E2dNdEGamma
        # plt.plot(EgammaList,E2dNdEGamma,*args,**kwargs)
        return EgammaList,np.array(E2dNdEGamma)

    ########################## Fit functions ##########################

    def power_law_exponent(self, E, I100, gamma, Ecut):
        return E**2*(I100*(E/.1)**-gamma)*np.exp(-E/Ecut)

    def power_law(self, E, I100, gamma):
        return E**2*I100*(E/.1)**-gamma

    def broken_power_law(self,E, IEc, a1, a2, Ecut):
        if np.shape(E) == ():
            return self.fit_function(E, IEc, a1, a2, Ecut)
        else:
            return [self.fit_function(val, IEc, a1, a2, Ecut) for val in E]  

    def fit_function(self, E, IEc, a1, a2, Ecut):
        if E < Ecut:
            y = E**2*IEc*(E/Ecut)**-a1
        elif E > Ecut:
            y = E**2*IEc*(E/Ecut)**-a2
        else:
            y = 0
        return y  

    ####################################################################

# This doesn't work yet 

class plot_3FGL:
    def __init__(self, bin_min, bin_max):
        self.fluxes_3FGL = np.loadtxt("/mnt/hepheno/CTBCORE/3FGL/fluxTabBin.dat")
        self.fluxes_summed = np.sum(self.fluxes_3FGL, axis=0)
        self.CTB_energies = 10**np.linspace(np.log10(0.3),np.log10(300),31)
        self.CTB_bin_centers = [10**((np.log10(self.CTB_energies[i])+np.log10(self.CTB_energies[i+1]))/2) for i in range(len(self.CTB_energies)-1)]
        self.CTB_bin_widths = [self.CTB_energies[i+1]-self.CTB_energies[i] for i in range(len(self.CTB_energies)-1)]

        self.E2dNdE = [self.CTB_bin_centers[i]**2*self.fluxes_summed[i]/self.CTB_bin_widths[i] for i in range(len(self.CTB_energies)-1)]
        # maps_dir = "/mnt/hepheno/CTBCORE/"
        # work_dir = "/group/hepheno/smsharma/Fermi_High_Latitudes/Fermi_2MASS_Edep"
        # f = fp.fermi_plugin(maps_dir,work_dir=work_dir,CTB_en_min = 0,CTB_en_max=30,nside=128)
        plt.plot(self.CTB_bin_centers, np.array(self.E2dNdE))
        # self.CTB_dir = "/mnt/hepheno/CTBCORE/PASS8_Jun15_UltracleanVeto_specbin/"
        # self.nside = 128
        # self.CTB_en_min = bin_min
        # self.CTB_en_max = bin_max
        # self.CTB_en_bins, self.CTB_count_maps, self.CTB_exposure_maps, self.CTB_psc_masks_temp = CTB.get_CTB(self.CTB_dir, self.nside, self.CTB_en_min, self.CTB_en_max,is_p8 = True)




        








       
