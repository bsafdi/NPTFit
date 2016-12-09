###############################################################################
# set_dirs.py
###############################################################################
#
# Define and create the directories required for the scan.
#
###############################################################################

import os


class SetDirs:  # pragma: no cover
    def __init__(self, tag='Untagged', work_dir=None):
        """ :param tag: label associated with the scan
            :param work_dir: location where all output from the run is stored
        """
        
        self.tag = tag

        # If unset, work_dir defaults to the current working directory, whilst
        # psf_dir is placed inside work_dir
        if work_dir is None:
            self.work_dir = os.getcwd() + '/'
        else:
            self.work_dir = work_dir

        # Chains directories, where the output of 
        self.chains_base_dir = self.work_dir + 'chains/'
        self.chains_dir = self.chains_base_dir + self.tag + '/'

        self.make_dirs([self.chains_base_dir, self.chains_dir])
        
    def make_dirs_for_run(self, run_tag=None):
        """ Set up the directory for the run itself, called at the time of 
            perform_scan
        """
        
        # If no run_tag specified, write directly into the chains directory
        if run_tag is None:
            self.chains_dir_for_run = self.chains_dir
        else:
            self.chains_dir_for_run = self.chains_dir + run_tag + '/'
            self.make_dirs([self.chains_dir_for_run])

    @staticmethod
    def make_dirs(dirs):
        """ Creates directories if they do not already exist 
        """

        for d in dirs:
            if not os.path.exists(d):
                try:
                    os.mkdir(d)
                except OSError as e:
                    if e.errno != 17:
                        raise   
