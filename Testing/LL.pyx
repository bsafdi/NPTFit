#------------------------------------------------
# NPTF Log Likelihood
#------------------------------------------------
##
## Calculation of Log Likelihood for a template
## fit involving non-poissonian templates
##
## Calculation broken into two parts:
##  1. Determination of x_m and x_m_sum, which
##     depends on the number of non-poissonian
##     templates & the number of breaks each one
##     has; and
##  2. Determination of LL in terms of these
##     using a recurrence relation
##
#------------------------------------------------

import numpy as np
cimport numpy as np
cimport cython
import incgamma_fct as igf

# Type used for all non-integer functions
DTYPE = np.float

# Setup cython functions
cdef extern from "math.h":
    double log(double x) nogil
    double exp(double x) nogil
    double pow(double x, double y) nogil


# Wrap the cython function as a python function that can be called externally 
def log_like(double[::1] xbg_PSF_compressed, 
             double[::1] theta, 
             double[::1] f_ary, 
             double[::1] df_rho_div_f_ary, 
             double[::1] PS_dist_compressed, 
             int[::1] data):

## SETUP SO IT PASSES SUM OF STUFF ON
    return log_like_internal(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data)


## Implement the recurrence relation and return LL
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_like_internal(double[::1] xbg_PSF_compressed, 
                              int[::1] data):

    ## Determine log likelihood for model to see data
    # This amounts to calculating p_k in every pixel
    # and adding their log likelihoods

    cdef int k_max = np.max(data) + 1
    cdef int npixROI = len(xbg_PSF_compressed)
    cdef size_t p, k
    cdef double ll = 0.
    cdef double f0_ary, f1_ary
    cdef double[:] pk = np.zeros((k_max+1), dtype=DTYPE)

    # Loop over pixels
    for p in range(npixROI):
        # Define p_0 (pk[0]) and p_1 (pk[1])
        # Then the remaining p_k are determined recursively
        # up to the value of k = data in that pixel
        f0_ary = -(xbg_PSF_compressed[p] + x_m_sum[p])
        f1_ary = (xbg_PSF_compressed[p] + x_m_ary[p,1])
        pk[0] = exp(f0_ary)
        pk[1] = pk[0] * f1_ary

        for k in range(2,data[p]+1):
            pk[k] = 0.
            for n in range(0, k - 1):
                pk[k] += (k-n)/ float(k) * x_m_ary[p,k-n] * pk[n]
            pk[k] += f1_ary * pk[k-1] / float(k)

        # Need prob > 0 to define a LL
        # If a very bad fit can get prob = 0,
        # if so then penalise to avoid this
        # region
        if pk[data[p]] > 0:
            ll += log(pk[data[p]])
        else:
            ll += -10.1**10.

    return ll
