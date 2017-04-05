###############################################################################
# incgamma_fct_p.pyx
###############################################################################
#
# Calculation of the two following arrays:
#  - Gamma(z+m,a)/m! - upper incomplete gamma
#  - Gamma(z+m,0,a)/m! - lower incomplete gamma
#
# Determine these using gamma function recursion relations. The use of 128
# point precision is important, and is why part of the code is python rather
# than cython.
#
# This is the python version of the gamma-functions and uses mpmath.
#
###############################################################################

import numpy as np
cimport numpy as np
from mpmath import mp

## Precalculate log(m!) for m in [0, fct_max - 1]
# If data map has more than fct_max counts per pixel fct_max must be increased
fct_max = 30000

def log_factorial(m):
    return np.sum(np.log(np.arange(1., m + 1., dtype=np.float128)))

log_factorial_ary = np.vectorize(log_factorial)(np.arange(fct_max))

def incgamma_up_fct_ary(np.int m_max, np.float z, np.float a):
    """ Upper incomplete gamma function / m!

    Calculation performed using gamma function recursion relation.

    Returns:
        float64 numpy array of Gamma(z+m,a)/m! for m in [0,m_max]

    """

    cdef np.ndarray igf_list = np.zeros(m_max + 1, dtype=np.float128)
    igf_list[0] = np.float128(mp.gammainc(z, a))

    # i_array[m] = m+1
    cdef np.ndarray i_array = np.arange(1., m_max + 1., dtype=np.float128)
    
    cdef Py_ssize_t i 
    for i in np.arange(1,m_max+1):
        igf_list[i] = ((i_array[i-1] - 1. + z)*igf_list[i-1]/i 
                      + np.exp((i_array[i-1] - 1. + z)*np.log(a) - a 
                      - log_factorial_ary[i]))

    return np.asarray(igf_list,dtype=np.float64)

def incgamma_lo_fct_ary(np.int m_max, np.float z, np.float a):
    """ Lower incomplete gamma function / m!

    Calculation performed using gamma function recursion relation.
    Starts at m_max and moves backwards for greater precision.

    Returns:
        float64 numpy array of Gamma(z+m,0,a)/m! for m in [0,m_max]

    """

    cdef np.ndarray igf_list = np.zeros(m_max + 1,dtype=np.float128)
    igf_list[m_max] = np.float128(mp.gammainc(z + m_max, 0, a) / mp.factorial(m_max))
    # NB: use mp.factorial, since np.math.factorial overflows for m ~ 170

    # i_array[m] = m_max-m
    cdef np.ndarray i_array = np.arange(m_max, 0, -1, dtype=np.float128)

    cdef Py_ssize_t i 
    for i in np.arange(m_max, 0, -1):
        igf_list[i-1] = ((igf_list[i] 
                        + np.exp((i_array[-i] - 1. + np.float128(z))
                        *np.log(np.float128(a)) - np.float128(a) 
                        - log_factorial_ary[i]))*i_array[-i]
                        /(i_array[-i] - 1. + np.float128(z)))

    return np.asarray(igf_list,dtype=np.float64)
