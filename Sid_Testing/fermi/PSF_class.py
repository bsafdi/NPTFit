import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib
from scipy.interpolate import interp1d as interp1d


class PSF:
    def __init__(self,f,theta_norm = [0.00000  ,
9.50266e-07  ,
0.000944168   , 
0.0155147    ,
0.0697261    , 
0.164378     ,
0.308687,
0.440748], rescale_index = 11, params_index = 10):
        self.f = f
        self.rescale_index = rescale_index
        self.params_index = params_index
        self.fill_rescale() 
        self.fill_PSF_params()
        self.fill_Kings()
        self.fill_68_and_95()
        if theta_norm == 'False':
            self.theta_norm = np.transpose([np.array([1 for i in range(8)]) for j in range(23)])
        else:
            self.theta_norm = np.transpose([theta_norm for i in range(23)])
            #print self.theta_norm

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
        
    def Kings(self,sigma,gamma,x):
        return (1/(2*np.pi*sigma**2))*(1-1/gamma)*(1+(1/(2.*gamma))*(x**2/sigma**2))**-gamma
    
    def fill_Kings(self):
        self.x_var = np.linspace(0,10,10000)
        self.fcore = np.array([[1/(1+self.NTAIL[i,j]*self.STAIL[i,j]**2/self.SCORE[i,j]**2) for j in range(np.shape(self.NCORE)[1])] for i in range(np.shape(self.NCORE)[0])])
        self.Kcore = np.array([[self.Kings(self.SCORE[i,j],self.GCORE[i,j],self.x_var) for j in range(np.shape(self.NCORE)[1])] for i in range(np.shape(self.NCORE)[0])])
        self.Ktail = np.array([[self.Kings(self.STAIL[i,j],self.GTAIL[i,j],self.x_var) for j in range(np.shape(self.NCORE)[1])] for i in range(np.shape(self.NCORE)[0])])
          
        self.PSF = np.array([[self.fcore[i,j]*self.Kcore[i,j] + (1-self.fcore[i,j])*self.Ktail[i,j] for j in range(np.shape(self.NCORE)[1])] for i in range(np.shape(self.NCORE)[0])])
    
    def find_perc(self,xvar,yvar,perc):
        norm = np.sum(xvar*yvar)
        for i in range(len(xvar)):
            sum_t = np.sum(xvar[0:i]*yvar[0:i])/norm
            if sum_t > perc:
                break
        return xvar[i]
    
    def fill_68_and_95(self):
        self.R68_x = np.array([[self.find_perc(self.x_var,2*np.pi*(self.x_var[1]-self.x_var[0])*self.PSF[i,j],0.68) for j in range(np.shape(self.NCORE)[1])] for i in range(np.shape(self.NCORE)[0])])
        self.R95_x = np.array([[self.find_perc(self.x_var,2*np.pi*(self.x_var[1]-self.x_var[0])*self.PSF[i,j],0.95) for j in range(np.shape(self.NCORE)[1])] for i in range(np.shape(self.NCORE)[0])])
        
        self.R68 = (360/(2*np.pi))*np.array([[self.R68_x[i,j]*self.rescale_factor((self.E_max[j]+self.E_min[j])/2.,self.rescale_array) for j in range(np.shape(self.NCORE)[1])] for i in range(np.shape(self.NCORE)[0])])
        self.R95 = (360/(2*np.pi))*np.array([[self.R95_x[i,j]*self.rescale_factor((self.E_max[j]+self.E_min[j])/2.,self.rescale_array) for j in range(np.shape(self.NCORE)[1])] for i in range(np.shape(self.NCORE)[0])])
        
    def rescale_factor(self,E,rescale): #E in MeV
        SpE = np.sqrt((rescale[0]*(E/100)**(rescale[2]))**2 + rescale[1]**2)
        return SpE
    
    def return_C68(self,E):
        
        return self.rescale_factor(E,self.rescale_array)
        
    def make_rescale_array(self,Emin,Emax,nbins):
        self.E_array = np.linspace(Emin,Emax,nbins)
        self.rescale_E_array = []
        self.rescale_E_array = np.array([self.return_C68(E) for E in self.E_array])*360/(2*np.pi) #in degrees
#         for rescale in self.rescale_array:
#             self.rescale_E_array.append(np.array([self.rescale_factor(E,rescale) for E in self.E_array]))
    
    def plot_rescale(self):
        self.rescale_fig = plt.figure(figsize=(8,6))
        plt.plot(self.E_array*10**-3,self.rescale_E_array,'r',label = 'P8R2_ULTRACLEANVETO_V6  (68%)')
        plt.legend(fontsize=13)
        plt.xlabel('energy [GeV]', fontsize=18)
        plt.ylabel('containment radius [degrees]', fontsize=18)
        plt.tick_params(axis='x', length=5,width=2,labelsize=18)
        plt.tick_params(axis='y',length=5,width=2,labelsize=18)
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(b=True, which='major', color='k', linestyle='-',linewidth=0.5)
        plt.grid(b=True, which='minor', color='k', linestyle='-',linewidth=0.25)
        plt.xlim([self.E_array[0]*10**-3,self.E_array[-1]*10**-3])
        #plt.title('PSF',fontsize=18)
        
    def plot_R(self):
        self.R_fig = plt.figure(figsize=(8,6))
        plt.plot((self.E_max+self.E_min)/2.*10**-3,np.sum(self.theta_norm*self.R68,axis=0),'r-',label = 'P8R2_ULTRACLEANVETO_V6  (68%)')
        plt.plot((self.E_max+self.E_min)/2.*10**-3,np.sum(self.theta_norm*self.R95,axis=0),'r--',label = 'P8R2_ULTRACLEANVETO_V6  (95%)')
        plt.legend(fontsize=13)
        plt.xlabel('energy [GeV]', fontsize=18)
        plt.ylabel('containment radius [degrees]', fontsize=18)
        plt.tick_params(axis='x', length=5,width=2,labelsize=18)
        plt.tick_params(axis='y',length=5,width=2,labelsize=18)
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(b=True, which='major', color='k', linestyle='-',linewidth=0.5)
        plt.grid(b=True, which='minor', color='k', linestyle='-',linewidth=0.25)

        #plt.title('PSF',fontsize=18)

    def interpolate_R(self):
        self.R68_int = interp1d((self.E_max+self.E_min)/2.*10**-3, np.sum(self.theta_norm*self.R68,axis=0))
        self.R95_int = interp1d((self.E_max+self.E_min)/2.*10**-3, np.sum(self.theta_norm*self.R95,axis=0))
        self.interp_bool = True

    def return_sigma_gaussian(self,energies,rad = '0.68'):
        if not self.interp_bool:
            self.interpolate_R()
        if rad=='0.68':
            return self.R68_int(energies) / 1.50959
        else:
            return self.R95_int(energies) / 2.44775
        
    def save_rescale_plot(self,dir,name):
        self.rescale_fig.savefig(dir+name)
        self.rescale_fig.close()
        
    def save_R_plot(self,dir,name):
        self.R_fig.savefig(dir+name)
        self.R_fig.close()
        
        
    def return_polyfit(self):
        self.fit_params = np.polyfit(np.log10((self.E_max-self.E_min)/2.*10**-3), np.log10(np.mean(self.R68,axis=0)), 1)
        self.lin_space_fit_params = [10**self.fit_params[0],self.fit_params[1]]
        return self.lin_space_fit_params

# class PSF:
#     def __init__(self,f):
#         self.f = f
#         self.fill_rescale() 
    
#     def fill_rescale(self):
#         rescale_index = [11,11,11,11]#[2,5,8,11]
#         self.rescale_array = []
#         for index in rescale_index:
#             self.rescale_array.append(self.f[index].data[0][0])
            
#     def rescale_factor(self,E,rescale): #E in MeV
#         SpE = np.sqrt((rescale[0]*(E/100)**(rescale[2]))**2 + rescale[1]**2)
#         return SpE
    
#     def return_C68(self,E):
#         self.log_E_array = np.linspace(0.75,6.5,4) #E in MeV
        
#         logE = np.log10(E)
#         if logE < self.log_E_array[0] or logE > self.log_E_array[-1]:
#             return 'energy out of range!'
#         else:
#             for i in range(len(self.log_E_array)):
#                 if logE < self.log_E_array[i]:
#                     bin_number = i-1
#             return self.rescale_factor(E,self.rescale_array[bin_number])

        
#     def make_rescale_array(self,Emin,Emax,nbins):
#         self.E_array = np.linspace(Emin,Emax,nbins)
#         self.rescale_E_array = []
#         self.rescale_E_array = np.array([self.return_C68(E) for E in self.E_array])*360/(2*np.pi) #in degrees
# #         for rescale in self.rescale_array:
# #             self.rescale_E_array.append(np.array([self.rescale_factor(E,rescale) for E in self.E_array]))
    
#     def plot_rescale(self):
#         self.rescale_fig = plt.figure(figsize=(8,6))
#         plt.plot(self.E_array*10**-3,self.rescale_E_array,'r',label = 'P8R2_ULTRACLEANVETO_V6  (68%)')
#         plt.legend(fontsize=13)
#         plt.xlabel('energy [GeV]', fontsize=18)
#         plt.ylabel('containment radius [degrees]', fontsize=18)
#         plt.tick_params(axis='x', length=5,width=2,labelsize=18)
#         plt.tick_params(axis='y',length=5,width=2,labelsize=18)
#         plt.yscale('log')
#         plt.xscale('log')
#         plt.grid(b=True, which='major', color='k', linestyle='-',linewidth=0.5)
#         plt.grid(b=True, which='minor', color='k', linestyle='-',linewidth=0.25)
#         plt.xlim([self.E_array[0]*10**-3,self.E_array[-1]*10**-3])
#         plt.close()
#         #plt.title('PSF',fontsize=18)
        
#     def save_rescale_plot(self,dir,name):
#         self.rescale_fig.savefig(dir+name)


#     def return_polyfit(self):
#         self.fit_params = np.polyfit(np.log10(self.E_array*10**-3), np.log10(self.rescale_E_array), 1)
#         self.lin_space_fit_params = [10**self.fit_params[1],self.fit_params[0]]
#         return self.lin_space_fit_params
