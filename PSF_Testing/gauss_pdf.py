#
# NAME:
#  king_pdf.py
#
# PURPOSE:
#  Take in a set of king function parameters, normalise the function and allow
#  sampling from a normalised king function PDF
#  For details of the Fermi PSF, see: http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_LAT_IRFs/IRF_PSF.html
#  Note the inverse sampling function doesn't require the input pdf to be normalised
#  Thus in this version of the code we removed this step
#
# REVISION HISTORY:
#  2015-Oct-02 Created by by Nick Rodd, MIT
#  2015-Oct-06 Removed the normalisation 

import numpy as np
import pdf_sampler as pdfs

class gauss_pdf:
    def __init__(self,sigma_PSF_deg,sampling=10000):
        self.sigma_PSF_deg = sigma_PSF_deg
        self.sampling = sampling

    def gauss_fn(self,r=0.1):
        # The combination of two king functions relevant for the Fermi PSF
        # NB: the x out the front is required as we treat this as a radius running from
        # 0 to infty, and this appears as a Jacobian factor
        # Similarly the 2pi appears from the angular integral
        return 2*np.pi*r*np.exp(-r**2/(2.0*(self.sigma_PSF_deg*np.pi/180.)**2))

    def gauss_pdf_draw(self,n):
        # Calcs an array of r values from the king pdf, then converts
        # them to deltap values by multiplying by SpE
        rvals = np.linspace(0,5*self.sigma_PSF_deg,self.sampling)
        pofr = self.gauss_fn(rvals)
        dist = pdfs.pdf_sampler(rvals,pofr)
        return dist(n)
