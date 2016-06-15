import sys, os

# current_dir = os.getcwd()
# change_path = ".."
# os.chdir(change_path)

import logging
from contextlib import contextmanager

bens_CTB_dir = '''/Users/bsafdi/Dropbox/Edep-NPTF/NPTF-shared/CTBCORE/'#'/Users/bsafdi/Dropbox/pulsars_local_new/CTBCORE/'''
erebus_CTB_dir = '/zfs/bsafdi/data/'
linas_CTB_dir = '/home/lina/Dropbox (MIT)/github/NPTF/CTBCORE/'
#gilly_CTB_dir = '/Users/gelor/Desktop/NPTF-shared/CTBCORE/'
gilly_CTB_dir = '/Users/gelor/Dropbox (MIT)/CTBCORE/'

class set_dirs:
    def __init__(self,tag = 'test',make_plots_data_tag = True,work_dir = '',psf_dir = 'False',**kwargs):
        self.psf_dir = psf_dir
        self.work_dir = work_dir
        self.tag = tag
        #self.data_tag = data_tag
        self.make_plots_data_tag = make_plots_data_tag


       # self.set_initial_dirs(self.tag,self.data_tag,**kwargs)
        self.define_additional_dirs(self.tag)
        self.make_dirs(self.dirs)

    # def set_initial_dirs(self,tag,data_tag,**kwargs):
    #     if 'maps_dir' in kwargs.keys():
    #         self.maps_dir = kwargs['maps_dir']
    #     elif os.path.exists(bens_CTB_dir):
    #         self.maps_dir    = bens_CTB_dir
    #     elif os.path.exists(erebus_CTB_dir):
    #         self.maps_dir    = erebus_CTB_dir
    #     elif os.path.exists(linas_CTB_dir):
    #         self.maps_dir    = linas_CTB_dir
    #     elif os.path.exists(gilly_CTB_dir):
    #         self.maps_dir    = gilly_CTB_dir

    #     self.CTB_dir     = self.maps_dir + self.data_tag + '/' 
    #     self.diff_dir    = self.maps_dir + 'diffuse/'
    #     self.ps_dir = self.maps_dir + '3FGL_masks_by_energy/'
    #     self.ps_flux_dir = self.maps_dir + '3FGL/'
    #     self.template_dir = self.maps_dir + 'additional_templates/'
    #     self.psf_data_dir = self.maps_dir + 'psf_data/'


        # plots_base_dir = dirs[2]
        # chains_base_dir = dirs[5]
        # psf_dir = dirs[0]
        # return maps_dir, CTB_dir, plots_base_dir, diff_dir, chains_base_dir, psf_dir

    def define_additional_dirs(self,tag):
        if self.psf_dir == 'False':
            self.psf_dir = self.work_dir+'psf/'
        #self.sim_dir = self.work_dir+'sim/'   #before, there was /full-sky/
        self.plots_base_dir = self.work_dir+'plots/'#full-sky/'   #before, there was /full-sky/
        self.chains_base_dir = self.work_dir+'''chains/'''#full-sky/' 
        #before, there was /full-sky/
        self.dict_base_dir = self.work_dir+'dict/'
        self.chains_dir = self.chains_base_dir + tag + '/'
        self.plots_dir = self.plots_base_dir + tag + '/'
        self.dict_dir = self.dict_base_dir + tag + '/'
        self.data_dir_base = self.work_dir+'data/'
        self.data_maps_dir = self.data_dir_base + 'maps/'
        #self.data_dir = self.data_dir_base + tag + '/'
        self.data_ps_dir = self.data_dir_base + 'ps_results/'


        self.dirs = [self.psf_dir, self.plots_base_dir, self.plots_dir, self.chains_base_dir, self.chains_dir,self.dict_base_dir, self.dict_dir, self.data_dir_base, self.data_maps_dir,self.data_ps_dir] #self.data_dir

    def make_dirs_for_run(self,run_tag,polychord=False):
        self.chains_dir_for_run = self.chains_dir + run_tag + '/'
        if self.make_plots_data_tag:
            self.plots_dir_for_run = self.plots_dir + run_tag + '/'
            self.data_dir_for_run  = self.chains_dir + run_tag + '/'
        else:
            self.plots_dir_for_run = self.plots_dir
            #self.data_dir_for_run = self.data_dir
        #self.dict_dir_for_run = self.dict_dir + run_tag + '/'
        self.dirs_for_run = [self.plots_dir_for_run, self.chains_dir_for_run] #self.data_dir_for_run  #,self.dict_dir_for_run]
        if polychord:
            self.dirs_for_run+=[self.chains_dir_for_run + 'clusters/']

        self.make_dirs(self.dirs_for_run)

    def make_dirs_for_ps_data(self):
        self.ps_data_dir = self.data_dir_base + 'ps_data/'
        self.make_dirs(self,[self.ps_data_dir])

    # def make_dirs(self,dirs):
    #     for d in dirs:
    #         if not os.path.exists(d):
    #             os.mkdir(d)


    def make_dirs(self,dirs):
        for d in dirs:
            if not os.path.exists(d):
                while True:
                    try:
                        os.mkdir(d)
                        break
                    except OSError, e:
                        if e.errno != 17:
                            raise   
                        pass


    def make_logging(self,plots_dir2):
        #======================================
        #set up logging to print messages both to screen and log
        #logger = logging.getLogger()  #uncomment for running ss iPython notebook
        #logger.setLevel(logging.INFO) #uncomment for running as iPython notebook
        #umcomment below for running as script
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(message)s',
                            datefmt='%m-%d %H:%M:%S',
                            filename=plots_dir2 + 'log.txt',
                            filemode='w')
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)
