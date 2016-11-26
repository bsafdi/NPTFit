###############################################################################
# pll.pyx
###############################################################################
#
# Log like for the purely Poissonian case.
#
###############################################################################

import numpy as np
cimport numpy as np
cimport cython

# Type used for all non-integer functions
DTYPE = np.float

# Setup cython functions
cdef extern from "math.h":
    double log(double x) nogil

# Precalculate log(k!) for m in [0, fct_max - 1]
# If data map has more than fct_max counts per pixel fct_max must be increased
fct_max = 5000 

def log_factorial(k):
    return np.sum(np.log(np.arange(1., k + 1., dtype=np.float128)))

log_factorial_ary = np.vectorize(log_factorial)(np.arange(fct_max))

def log_like_poissonian(double[::1] pt_sum_compressed, int[::1] data):
    """ Python wrapper for the log likelihood

    Args:
        pt_sum_compressed: pixel-wise sum of Poissonian model templates
        data: The pixel-wise data

    Returns:
        double log likelihood

    """

    return log_like_poissonian_int(pt_sum_compressed, data)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_like_poissonian_int(double[::1] pt_sum_compressed,
                                    int[::1] data):
    """ Calculation of the log likelihood
        
        Full log likelihood is the sum of the value in all pixels

        Returns:
            double log likelihood
    """
    
    cdef Py_ssize_t p
    cdef double ll = 0.0
    cdef int npix_roi = len(pt_sum_compressed)
    for p in range(npix_roi):
        ll += log_poisson_int(pt_sum_compressed[p], data[p])
    return ll

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_poisson_int(double mean, int k):
    """ Poisson likelihood evaluation for an observed number of counts k, given
        an expected counts mean
    """

    return -mean + k*log(mean) - log_factorial_ary[k]
