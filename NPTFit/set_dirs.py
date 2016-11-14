import os

class set_dirs:
    def __init__(self,tag = 'test',work_dir = '',psf_dir = 'False',**kwargs):
        self.psf_dir = psf_dir
        self.work_dir = work_dir
        self.tag = tag

        self.define_additional_dirs(self.tag)
        self.make_dirs(self.dirs)

    def define_additional_dirs(self,tag):
        if self.psf_dir == 'False':
            self.psf_dir = self.work_dir+'psf/'
        self.chains_base_dir = self.work_dir+'''chains/'''
        self.chains_dir = self.chains_base_dir + tag + '/'


        self.dirs = [self.psf_dir, self.chains_base_dir, self.chains_dir] 

    def make_dirs_for_run(self,run_tag):
        self.chains_dir_for_run = self.chains_dir + run_tag + '/'
        self.data_dir_for_run  = self.chains_dir + run_tag + '/'
        self.dirs_for_run = [self.chains_dir_for_run] 

        self.make_dirs(self.dirs_for_run)

    def make_dirs(self,dirs):
        for d in dirs:
            if not os.path.exists(d):
                try:
                    os.mkdir(d)
                except OSError, e:
                    if e.errno != 17:
                        raise   
