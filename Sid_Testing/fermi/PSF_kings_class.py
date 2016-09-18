import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib
from scipy.interpolate import interp1d as interp1d


class PSF_kings:
    """ A class to return various Kings function parameters """
    def __init__(self,f,theta_norm = [0.00000  ,
9.50266e-07  ,
0.000944168   , 
0.0155147    ,
0.0697261    , 
0.164378     ,
0.308687,
0.440748], rescale_index = 11, params_index = 10):

        print "All input energies must be in GeV"
        self.f = f
        self.rescale_index = rescale_index
        self.params_index = params_index
        self.fill_rescale() 
        self.fill_PSF_params()
        if theta_norm == 'False':
        	self.theta_norm = np.transpose([np.array([1 for i in range(8)]) for j in range(23)])
    	else:
    		self.theta_norm = np.transpose([theta_norm for i in range(23)])
        self.interp_bool = False
        self.interpolate_R()
    
    def fill_rescale(self):
        #rescale_index = 11 #2
        self.rescale_array = self.f[self.rescale_index].data[0][0]
        
    def fill_PSF_params(self):
        #params_index = 10 #1
        self.E_min = self.f[self.params_index].data[0][0] #size is 23
        self.E_max = self.f[self.params_index].data[0][1]
        self.theta_min = self.f[self.params_index].data[0][2] #size is 8
        self.theta_max = self.f[self.params_index].data[0][3]
        self.NCORE = np.array(self.f[self.params_index].data[0][4]) #shape is (8, 23)
        self.NTAIL = np.array(self.f[self.params_index].data[0][5])
        self.SCORE = np.array(self.f[self.params_index].data[0][6])
        self.STAIL = np.array(self.f[self.params_index].data[0][7])
        self.GCORE = np.array(self.f[self.params_index].data[0][8])
        self.GTAIL = np.array(self.f[self.params_index].data[0][9])
        # Now create fcore from the definition
        self.FCORE = np.array([[1/(1+self.NTAIL[i,j]*self.STAIL[i,j]**2/self.SCORE[i,j]**2) for j in range(np.shape(self.NCORE)[1])] for i in range(np.shape(self.NCORE)[0])])
         
    def rescale_factor(self,E): #E in GeV
        SpE = np.sqrt((self.rescale_array[0]*(E*10**3/100)**(self.rescale_array[2]))**2 + self.rescale_array[1]**2)
        return SpE
    
    def interpolate_R(self):
        self.FCORE_int = interp1d((self.E_max+self.E_min)/2.*10**-3, np.sum(self.theta_norm*self.FCORE,axis=0))
        self.SCORE_int = interp1d((self.E_max+self.E_min)/2.*10**-3, np.sum(self.theta_norm*self.SCORE,axis=0))
        self.STAIL_int = interp1d((self.E_max+self.E_min)/2.*10**-3, np.sum(self.theta_norm*self.STAIL,axis=0))
        self.GCORE_int = interp1d((self.E_max+self.E_min)/2.*10**-3, np.sum(self.theta_norm*self.GCORE,axis=0))
        self.GTAIL_int = interp1d((self.E_max+self.E_min)/2.*10**-3, np.sum(self.theta_norm*self.GTAIL,axis=0))
        self.interp_bool = True

    def return_kings_params(self,energies,param): # Put E in in GeV
        if not self.interp_bool:
            self.interpolate_R()
        if param=='fcore':
            return self.FCORE_int(energies)
        elif param=='score':
            return self.SCORE_int(energies)
        elif param=='gcore':
            return self.GCORE_int(energies)
        elif param=='stail':
            return self.STAIL_int(energies)
        elif param=='gtail':
            return self.GTAIL_int(energies)
        else:
            print "Param must be fcore, score, gcore, stail or gtail"

