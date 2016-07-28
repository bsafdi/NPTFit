#------------------------------------------------
# Incomplete Gamma Function / Factorial Arrays
#------------------------------------------------
##
## Calculation of Log Likelihood for a template
## fit involving non-poissonian templates
##
#------------------------------------------------

import numpy as np
cimport numpy as np
from mpmath import mp

sp_max = 5000

def log_factorial(k):
    #return log(k!)
    return np.sum(np.log(np.arange(1., k + 1., dtype=np.float128)))

#precalculate m! and log(m!) arrays for m in [0, sp_max - 1]
#need to use mp.factorial, since np.math.factorial overflows for m ~ 170
log_factorial_ary = np.vectorize(log_factorial)(np.arange(sp_max))

gammainc_mem = mp.gammainc
   
def gammainc_up_fct_ary_log(np.int m_max, np.float z, np.float a):
    #calculate in log-space an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
    cdef np.ndarray gi_list = np.zeros(m_max + 1, dtype=np.float128)
    gi_list[0] = np.float128(gammainc_mem(z, a))

    cdef np.ndarray i_array = np.arange(1., m_max + 1., dtype=np.float128)
    
    cdef Py_ssize_t i 
    for i in np.arange(1, m_max + 1):
        gi_list[i] = (i_array[i-1] - 1. + z)*gi_list[i-1]/i + np.exp((i_array[i-1] - 1. + z)*np.log(a) - a - log_factorial_ary[i])
    return np.asarray(gi_list,dtype=np.float64)

def gammainc_lo_fct_ary_back(np.int m_max, np.float z, np.float  b):
    #calculate in log-space an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
    #use backwards recursion to avoid loss of precision from subtraction
    cdef np.ndarray gi_list = np.zeros(m_max + 1,dtype=np.float128)
    gi_list[m_max] = np.float128(gammainc_mem(z + m_max, 0, b) / mp.factorial(m_max))

    cdef np.ndarray i_array = np.arange(m_max, 0, -1, dtype=np.float128)

    cdef Py_ssize_t i 
    for i in np.arange(m_max, 0, -1):
        gi_list[i-1] = (gi_list[i] + np.exp((i_array[-i] - 1. + np.float128(z))*np.log(np.float128(b)) - np.float128(b) - log_factorial_ary[i]))*i_array[-i]/(i_array[-i] - 1. + np.float128(z))
    return np.asarray(gi_list,dtype=np.float64)
