import os
import numpy as np
import pulsars.psf_king as psfk



class make_f_ary_kings:
    def __init__(self,nside,psf_dir,fcore,score,gcore,stail,gtail,SpE,frac=[1],save_tag='temp',**kwargs):
        print 'The save_tag is ', save_tag
        self.fcore = fcore
        self.score = score
        self.gcore = gcore
        self.stail = stail
        self.gtail = gtail
        self.SpE = SpE
        self.frac = frac

        self.save_tag=save_tag
        self.psf_dir = psf_dir
        self.nside = nside
        self.kwargs=kwargs
        if 'num_f_bins' in kwargs.keys():
            self.num_f_bins = kwargs['num_f_bins']
        else:
            self.num_f_bins = 10
            kwargs['num_f_bins'] = 10

        self.make_f_ary_kings()


    def make_f_ary_kings(self):
        self.psfk_inst = psfk.psf_king(self.fcore,self.score,self.gcore,self.stail,self.gtail,self.SpE,frac=self.frac,nside=self.nside,psf_dir = self.psf_dir,save_tag=self.save_tag,**self.kwargs) #kwargs example: psf_dir=self.psf_dir,nside=self.nside,num_f_bins=self.num_f_bins,n_ps=self.n_ps,n_pts_per_king=self.n_pts_per_psf,save_tag=self.psf_save_tag


        self.f_ary_file = self.psf_dir+'f_ary-' + self.save_tag + '-'+ str(self.nside) + '-' + str(self.num_f_bins) + '.dat'
        if not os.path.exists(self.f_ary_file):
            print 'we have to make PSF ...'
            self.f_ary, self.df_rho_div_f_ary = self.psfk_inst.make_PSF()
            print 'finished making PSF ...'
        else:
            print 'just loading PSF ...'
            self.f_ary, self.df_rho_div_f_ary = self.psfk_inst.load_PSF()