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

import sys, os
import numpy as np
from scipy import integrate

class king_pdf:
    """ A class to return a normalised king function """ 
    def __init__(self,fcore,score,gcore,stail,gtail,SpE,king_sampling=10000):
        self.fcore = fcore
        self.score = score
        self.gcore = gcore
        self.stail = stail
        self.gtail = gtail
        self.SpE = SpE
        self.king_sampling = king_sampling
        self.king_fn_base()
        self.king_fn_full()

    def king_fn_base(self,x=0.1,sigma=1,gamma=2):
        # A basic king function
        return (1/(2*np.pi*sigma**2))*(1-1/gamma)*(1+x**2/(2*gamma*sigma**2))**(-gamma)

    def king_fn_full(self,r=0.1):
        # The combination of two king functions relevant for the Fermi PSF
        # NB: the x out the front is required as we treat this as a radius running from
        # 0 to infty, and this appears as a Jacobian factor
        # Similarly the 2pi appears from the angular integral
        return 2*np.pi*r*(self.fcore*self.king_fn_base(r/self.SpE,self.score,self.gcore)+(1-self.fcore)*self.king_fn_base(r/self.SpE,self.stail,self.gtail))

    def king_pdf(self,n):
        # Calcs an array of r values from the king pdf, then converts
        # them to deltap values by multiplying by SpE
        rarr = np.linspace(0,10*self.SpE*(self.score+self.stail)/2.,self.king_sampling)
        pdf = self.king_fn_full(rarr)
        dist = Distribution(pdf, interpolation=False, transform=lambda i:[rarr[idr] for idr in np.array(i)])
        return dist(n)[0]

# class for distribution sampling taken from stack exchange
class Distribution(object):
    """
    draws samples from a one dimensional probability distribution,
    by means of inversion of a discrete inverstion of a cumulative density function

    the pdf can be sorted first to prevent numerical error in the cumulative sum
    this is set as default; for big density functions with high contrast,
    it is absolutely necessary, and for small density functions,
    the overhead is minimal

    a call to this distibution object returns indices into density array
    """
    def __init__(self, pdf, sort = True, interpolation = True, transform = lambda x: x):
        self.shape          = pdf.shape
        self.pdf            = pdf.ravel()
        self.sort           = sort
        self.interpolation  = interpolation
        self.transform      = transform

        #a pdf can not be negative
        assert(np.all(pdf>=0))

        #sort the pdf by magnitude
        if self.sort:
            self.sortindex = np.argsort(self.pdf, axis=None)
            self.pdf = self.pdf[self.sortindex]
        #construct the cumulative distribution function
        self.cdf = np.cumsum(self.pdf)
    @property
    def ndim(self):
        return len(self.shape)
    @property
    def sum(self):
        """cached sum of all pdf values; the pdf need not sum to one, and is imlpicitly normalized"""
        return self.cdf[-1]
    def __call__(self, N):
        """draw """
        #pick numbers which are uniformly random over the cumulative distribution function
        choice = np.random.uniform(high = self.sum, size = N)
        #find the indices corresponding to this point on the CDF
        index = np.searchsorted(self.cdf, choice)
        #if necessary, map the indices back to their original ordering
        if self.sort:
            index = self.sortindex[index]
        #map back to multi-dimensional indexing
        index = np.unravel_index(index, self.shape)
        index = np.vstack(index)
        #is this a discrete or piecewise continuous distribution?
        if self.interpolation:
            index = index + np.random.uniform(size=index.shape)
        return self.transform(index)
