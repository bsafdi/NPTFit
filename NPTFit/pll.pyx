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

# Setup cython functions
cdef extern from "math.h":
    double log(double x) nogil
    double lgamma(double x) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cpdef double log_like_poissonian(double[::1] pt_sum_compressed,
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
        ll += -pt_sum_compressed[p] + data[p]*log(pt_sum_compressed[p]) \
              - lgamma(data[p] + 1.)
        
    return ll
