#import sys, os

#current_dir = os.getcwd()
#change_path = ".."
#os.chdir(change_path)

import numpy as np
import scipy
# import matplotlib
# #matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib as mpl

import healpy as hp
# import mpmath as mp
# from astropy.io import fits
# import time

import pulsars.special as spc
# import pulsars.masks as masks
# import pulsars.CTB as CTB
# import pulsars.psf as psf
# import pulsars.diffuse_fermi as df
# import pulsars.likelihood_psf_pdep as llpsf

#import numpy as np
cimport numpy as np
cimport cython
#from libcpp cimport bool

import pulsars.special as spc

import logging
logger = logging.getLogger(__name__)

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef extern from "math.h":
    double log(double x) nogil
    double exp(double x) nogil
    double pow(double x, double y) nogil


# from mpmath import mp

# sp_max = 5000 #was 200

# def log_factorial(k):
#     #return log(k!)]
#     return np.sum(np.log(np.arange(1., k + 1., dtype=np.float128)))

# #precalculate m! and log(m!) arrays for m in [0, sp_max - 1]
# #need to use mp.factorial, since np.math.factorial overflows for m ~ 170
# factorial_ary = np.array([np.float128(mp.factorial(m)) for m in np.arange(sp_max)])

sp_max = 5000 

def log_factorial(k):
    return np.sum(np.log(np.arange(1., k + 1., dtype=np.float)))

# cdef double[::1] factorial_ary = np.array([np.float(mp.factorial(m)) for m in np.arange(sp_max)])
cdef double[::1] log_factorial_ary = np.vectorize(log_factorial)(np.arange(sp_max))

#cdef double[::1] log_factorial_ary = spc.log_factorial_ary.astype(dtype='float')

def log_python(mean, k):
    return np.sum( -mean + k*np.log(mean) )

def log_cython_vect(double[::1] mean, int[::1] k):
    return log_cython_vect_int(mean, k)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_cython_int(double mean, int k) nogil:
    return -mean + k*log(mean) - log_factorial_ary[k]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_cython_vect_int(double[::1] mean, int[::1] k):
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(mean)
    cdef double ll=0.0
    for i in range(n):
        ll+=log_cython_int(mean[i],k[i])
    return ll

def log_poisson(double mean, int k):
    return log_poisson_int(mean,k)
    #return -mean + k*np.log(mean) - spc.log_factorial_ary[k.astype(np.int)]

def log_poisson_old(mean, k):
    #return log_poisson_int(mean,k)
    return -mean + k*np.log(mean) - spc.log_factorial_ary[k.astype(np.int)]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_poisson_int(double mean, int k) nogil:
    return -mean + k*log(mean) - log_factorial_ary[k]


def log_likelihood_poissonian(double[::1] back, int[::1] data):
    return log_likelihood_poissonian_int(back,data)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_likelihood_poissonian_int(double[::1] back, int[::1] data):
    cdef double ll = 0.0
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(back)
    for i in range(n):
        ll += log_poisson_int(back[i], data[i])
    return ll



