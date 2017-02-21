import numpy as np
cimport numpy as np
cimport cython

cdef extern from "incgamma_fct.h":
    void inc_gamma_upper_ratio_ary(int x_max, double y, double a, double* result) nogil
    void inc_gamma_lower_ratio_ary(int x_max, double y, double a, double* result) nogil

cdef np.ndarray incgamma_up_fct_ary(int m_max, double z, double a):
    """ Upper incomplete gamma function / m!
    Calculation performed using gamma function recursion relation.
    Returns:
        float64 numpy array of Gamma(z+m,a)/m! for m in [0,m_max]
    """
    cdef np.ndarray[double, ndim=1, mode="c"] result = np.zeros(m_max + 1, dtype=np.double, order='c')
    inc_gamma_upper_ratio_ary(m_max+1, z-1, a, &result[0])
    return result

def incgamma_up_fct_ary_py(np.int m_max, np.float z, np.float a):
    """ Python interface, don't call this from Cython code """
    return incgamma_up_fct_ary(m_max, z, a)

cdef np.ndarray incgamma_lo_fct_ary(int m_max, double z, double a):
    """ Lower incomplete gamma function / m!
    Calculation performed using gamma function recursion relation.
    Starts at m_max and moves backwards for greater precision.
    Returns:
        float64 numpy array of Gamma(z+m,0,a)/m! for m in [0,m_max]
    """
    cdef np.ndarray[double, ndim=1, mode="c"] result = np.zeros(m_max + 1, dtype=np.double, order='c')
    inc_gamma_lower_ratio_ary(m_max+1, z-1, a, &result[0])
    return result

def incgamma_lo_fct_ary_py(np.int m_max, np.float z, np.float a):
    """ Python interface, don't call this from Cython code """
    return incgamma_lo_fct_ary(m_max, z, a)
