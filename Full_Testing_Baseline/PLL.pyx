#------------------------------------------------------------------------------
# Poissonian Log Likelihood
#------------------------------------------------------------------------------
##
## Log likelihood for the purely poissonian case.
##
#------------------------------------------------------------------------------

# TODO: currently written as an edep version (hence sum of PT_sum_compressed and data indices).
# Shoudl simplify once gotten rid of edep in main code. 

import numpy as np
cimport numpy as np
cimport cython

# Type used for all non-integer functions
DTYPE = np.float

# Setup cython functions
cdef extern from "math.h":
    double log(double x) nogil

sp_max = 5000 

def log_factorial(k):
    return np.sum(np.log(np.arange(1., k + 1., dtype=DTYPE)))

cdef double[::1] log_factorial_ary = np.vectorize(log_factorial)(np.arange(sp_max))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_poisson_int(double mean, int k) nogil:
    return -mean + k*log(mean) - log_factorial_ary[k]

def log_likelihood_poissonian(double[::1] PT_sum_compressed, int[::1] data):
    return log_likelihood_poissonian_int(PT_sum_compressed,data)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_likelihood_poissonian_int(double[::1] PT_sum_compressed, int[::1] data):
    cdef double ll = 0.0
    cdef Py_ssize_t i
    cdef Py_ssize_t n = len(PT_sum_compressed)
    for i in range(n):
        ll += log_poisson_int(PT_sum_compressed[i], data[i])
    return ll
