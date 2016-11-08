import numpy as np
cimport numpy as np
from mpmath import mp

sp_max = 5000 #was 200

def log_factorial(k):
    #return log(k!)]
    return np.sum(np.log(np.arange(1., k + 1., dtype=np.float128)))

#precalculate m! and log(m!) arrays for m in [0, sp_max - 1]
#need to use mp.factorial, since np.math.factorial overflows for m ~ 170
factorial_ary = np.array([np.float128(mp.factorial(m)) for m in np.arange(sp_max)])
log_factorial_ary = np.vectorize(log_factorial)(np.arange(sp_max))

#memoize gamma functions from mpmath to increase speed
# gamma_mem = mp.memoize(mp.gamma)
# gammainc_mem = mp.memoize(mp.gammainc)
gamma_mem = mp.gamma
gammainc_mem = mp.gammainc
   
def gammainc_up_fct_ary(m_max, z, a):
    #calculate naively an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
    #for typical PL parameters, overflows for m ~ 250
    gi_list = np.zeros(m_max + 1, dtype=np.float128)
#    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)
    gi_list[0] = np.float128(gammainc_mem(z, a))
    for i in np.arange(1, m_max + 1):
        gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.power(a, i - 1. + z) * np.exp(-a)/factorial_ary[i]
    return gi_list

def gammainc_lo_fct_ary(m_max, z, b):
    #calculate naively an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
    #subtraction of nearly equal numbers causes loss of precision around m ~ 100
    gi_list = np.zeros(m_max + 1, dtype=np.float128)
    gi_list[0] = np.float128(gammainc_mem(z, 0, b))
    for i in np.arange(1, m_max + 1):
        gi_list[i] = (i - 1. + z)*gi_list[i-1]/i - np.power(b, i - 1. + z) * np.exp(-b)/factorial_ary[i]
    return gi_list
    
def gammainc_up_fct_ary_log(np.int m_max, np.float z, np.float a):
    #calculate in log-space an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
    cdef np.ndarray gi_list = np.zeros(m_max + 1, dtype=np.float128)
#    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)    
    gi_list[0] = np.float128(gammainc_mem(z, a))

    #cdef np.float i_t = 0.0
    cdef np.ndarray i_array = np.arange(1., m_max + 1., dtype=np.float128)
    # for i in np.arange(1., m_max + 1., dtype=np.float128):
    cdef Py_ssize_t i 
    for i in np.arange(1, m_max + 1):
        gi_list[i] = (i_array[i-1] - 1. + z)*gi_list[i-1]/i + np.exp((i_array[i-1] - 1. + z)*np.log(a) - a - log_factorial_ary[i])
        # gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.exp((i - 1. + z)*np.log(a) - a - log_factorial_ary[i])
    return np.asarray(gi_list,dtype=np.float64)


def gammainc_up_fct_ary_log_old(m_max, z, a):
    #calculate in log-space an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
    gi_list = np.zeros(m_max + 1, dtype=np.float128)
#    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)    
    gi_list[0] = np.float128(gammainc_mem(z, a))
    for i in np.arange(1., m_max + 1., dtype=np.float128):
        gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.exp((i - 1. + z)*np.log(a) - a - log_factorial_ary[i])
    return gi_list

# def gammainc_up_fct_ary_mid(m_max, z, a, b):
#     #calculate in log-space an array of middle incomplete gamma function divided by factorial (gamma(z + m, a, b)/m!) for m in [0, m_max]
#     gi_list = np.zeros(m_max + 1, dtype=np.float128)
# #    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)    
#     gi_list[0] = np.float128(gammainc_mem(z, a, b))
#     for i in np.arange(1., m_max + 1., dtype=np.float128):
#         gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.exp((i - 1. + z)*np.log(a) - a - log_factorial_ary[i])
#     return gi_list
    
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


def gammainc_lo_fct_ary_back_old(m_max,z,  b):
    #calculate in log-space an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
    #use backwards recursion to avoid loss of precision from subtraction
    gi_list = np.zeros(m_max + 1)
    gi_list[m_max] = np.float128(gammainc_mem(z + m_max, 0, b) / mp.factorial(m_max))
    b = np.float128(b)
    for i in np.arange(m_max, 0, -1, dtype=np.float128):
        gi_list[i-1] = (gi_list[i] + np.exp((i - 1. + z)*np.log(b) - b - log_factorial_ary[i]))*i/(i - 1. + z)
    return gi_list


# import numpy as np
# cimport numpy as np
# from mpmath import mp

# sp_max = 5000 #was 200

# def log_factorial(k):
#     #return log(k!)]
#     return np.sum(np.log(np.arange(1., k + 1., dtype=np.float)))

# #precalculate m! and log(m!) arrays for m in [0, sp_max - 1]
# #need to use mp.factorial, since np.math.factorial overflows for m ~ 170
# factorial_ary = np.array([np.float(mp.factorial(m)) for m in np.arange(sp_max)])
# log_factorial_ary = np.vectorize(log_factorial)(np.arange(sp_max))

# #memoize gamma functions from mpmath to increase speed
# gamma_mem = mp.memoize(mp.gamma)
# gammainc_mem = mp.memoize(mp.gammainc)
   
# def gammainc_up_fct_ary(m_max, z, a):
#     #calculate naively an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
#     #for typical PL parameters, overflows for m ~ 250
#     gi_list = np.zeros(m_max + 1, dtype=np.float)
# #    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)
#     gi_list[0] = np.float(gammainc_mem(z, a))
#     for i in np.arange(1, m_max + 1):
#         gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.power(a, i - 1. + z) * np.exp(-a)/factorial_ary[i]
#     return gi_list

# def gammainc_lo_fct_ary(m_max, z, b):
#     #calculate naively an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
#     #subtraction of nearly equal numbers causes loss of precision around m ~ 100
#     gi_list = np.zeros(m_max + 1, dtype=np.float)
#     gi_list[0] = np.float(gammainc_mem(z, 0, b))
#     for i in np.arange(1, m_max + 1):
#         gi_list[i] = (i - 1. + z)*gi_list[i-1]/i - np.power(b, i - 1. + z) * np.exp(-b)/factorial_ary[i]
#     return gi_list
    
# def gammainc_up_fct_ary_log(np.int m_max, np.float z, np.float a):
#     #calculate in log-space an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
#     cdef np.ndarray gi_list = np.zeros(m_max + 1, dtype=np.float)
# #    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)    
#     gi_list[0] = np.float(gammainc_mem(z, a))

#     #cdef np.float i_t = 0.0
#     cdef np.ndarray i_array = np.arange(1., m_max + 1., dtype=np.float)
#     # for i in np.arange(1., m_max + 1., dtype=np.float128):
#     cdef Py_ssize_t i 
#     for i in np.arange(1, m_max + 1):
#         gi_list[i] = (i_array[i-1] - 1. + z)*gi_list[i-1]/i + np.exp((i_array[i-1] - 1. + z)*np.log(a) - a - log_factorial_ary[i])
#         # gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.exp((i - 1. + z)*np.log(a) - a - log_factorial_ary[i])
#     return gi_list


# def gammainc_up_fct_ary_log_old(m_max, z, a):
#     #calculate in log-space an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
#     gi_list = np.zeros(m_max + 1, dtype=np.float128)
# #    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)    
#     gi_list[0] = np.float128(gammainc_mem(z, a))
#     for i in np.arange(1., m_max + 1., dtype=np.float128):
#         gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.exp((i - 1. + z)*np.log(a) - a - log_factorial_ary[i])
#     return gi_list

# # def gammainc_up_fct_ary_mid(m_max, z, a, b):
# #     #calculate in log-space an array of middle incomplete gamma function divided by factorial (gamma(z + m, a, b)/m!) for m in [0, m_max]
# #     gi_list = np.zeros(m_max + 1, dtype=np.float128)
# # #    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)    
# #     gi_list[0] = np.float128(gammainc_mem(z, a, b))
# #     for i in np.arange(1., m_max + 1., dtype=np.float128):
# #         gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.exp((i - 1. + z)*np.log(a) - a - log_factorial_ary[i])
# #     return gi_list
    
# def gammainc_lo_fct_ary_back(np.int m_max, np.float z, np.float  b):
#     #calculate in log-space an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
#     #use backwards recursion to avoid loss of precision from subtraction
#     cdef np.ndarray gi_list = np.zeros(m_max + 1,dtype=np.float)
#     gi_list[m_max] = np.float(gammainc_mem(z + m_max, 0, b) / mp.factorial(m_max))

#     cdef np.ndarray i_array = np.arange(m_max, 0, -1, dtype=np.float)

#     cdef Py_ssize_t i 
#     for i in np.arange(m_max, 0, -1):
#         gi_list[i-1] = (gi_list[i] + np.exp((i_array[-i] - 1. + np.float(z))*np.log(np.float(b)) - np.float(b) - log_factorial_ary[i]))*i_array[-i]/(i_array[-i] - 1. + np.float(z))
#     return gi_list


# def gammainc_lo_fct_ary_back_old(m_max,z,  b):
#     #calculate in log-space an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
#     #use backwards recursion to avoid loss of precision from subtraction
#     gi_list = np.zeros(m_max + 1)
#     gi_list[m_max] = np.float128(gammainc_mem(z + m_max, 0, b) / mp.factorial(m_max))
#     b = np.float128(b)
#     for i in np.arange(m_max, 0, -1, dtype=np.float128):
#         gi_list[i-1] = (gi_list[i] + np.exp((i - 1. + z)*np.log(b) - b - log_factorial_ary[i]))*i/(i - 1. + z)
#     return gi_list



















###################
###This is faster but not enough precision


# import numpy as np
# cimport numpy as np
# cimport cython
# from mpmath import mp

# cdef extern from "math.h":
#     long double log(long double x) nogil
#     long double exp(long double x) nogil

# sp_max = 5000 

# def log_factorial(k):
#     return np.sum(np.log(np.arange(1., k + 1., dtype=np.float128)))

# factorial_ary = np.array([np.float128(mp.factorial(m)) for m in np.arange(sp_max)])
# log_factorial_ary = np.vectorize(log_factorial)(np.arange(sp_max))

# gamma_mem = mp.gamma
# gammainc_mem = mp.gammainc

# def gammainc_up_fct_ary_log(m_max,double z, double a):
#     return gammainc_up_fct_ary_log_impl(m_max, z, a)

# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# cdef np.ndarray gammainc_up_fct_ary_log_impl(int m_max, long double z_t, long double a_t):
#     cdef long double z = <long double>z_t
#     cdef long double a = <long double>a_t
#     cdef long double[::1] gi_list = np.zeros(m_max + 1, dtype=np.float128)
#     gi_list[0] = gammainc_mem(z, a)
#     cdef long double t0
#     cdef long double t1
#     cdef Py_ssize_t i
#     for i in range(1, m_max + 1):
#         t0 = (i - 1. + z)
#         t1 = (i - 1. + z)*log(a) - a
#         gi_list[i] = t0*gi_list[i-1]/i + exp(t1 - log_factorial_ary[i])
#     return np.asarray(gi_list,np.float64)


# def gammainc_lo_fct_ary_back(m_max, double z, double b):
#     return gammainc_lo_fct_ary_back_impl(m_max, z, b)

# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# cdef np.ndarray gammainc_lo_fct_ary_back_impl(int m_max, double z_t, double b_t):
#     #calculate in log-space an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
#     #use backwards recursion to avoid loss of precision from subtraction
#     cdef long double z = <long double>z_t
#     cdef long double b = <long double>b_t
#     cdef long double[::1] gi_list = np.zeros(m_max + 1,dtype=np.float128)
#     cdef long double t0
#     cdef long double t1
#     gi_list[m_max] = gammainc_mem(z + m_max, 0, b) / mp.factorial(m_max)

#     #cdef np.ndarray i_array = np.arange(m_max, 0, -1, dtype=np.float)

#     cdef Py_ssize_t i 
#     for i in range(m_max, 0, -1):
#         t0 = i/(i - 1. + z) 
#         t1 = exp(i-1+z)*log(b) -b
#         gi_list[i-1] = (gi_list[i] + t1 - log_factorial_ary[i])*t0
#     return np.asarray(gi_list,np.float64)


# import numpy as np
# cimport numpy as np
# from mpmath import mp

# sp_max = 5000 #was 200

# def log_factorial(k):
#     #return log(k!)]
#     return np.sum(np.log(np.arange(1., k + 1., dtype=np.float128)))

# #precalculate m! and log(m!) arrays for m in [0, sp_max - 1]
# #need to use mp.factorial, since np.math.factorial overflows for m ~ 170
# factorial_ary = np.array([np.float128(mp.factorial(m)) for m in np.arange(sp_max)])
# log_factorial_ary = np.vectorize(log_factorial)(np.arange(sp_max))

# #memoize gamma functions from mpmath to increase speed
# # gamma_mem = mp.memoize(mp.gamma)
# # gammainc_mem = mp.memoize(mp.gammainc)
# gamma_mem = mp.gamma
# gammainc_mem = mp.gammainc
   
# def gammainc_up_fct_ary(m_max, z, a):
#     #calculate naively an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
#     #for typical PL parameters, overflows for m ~ 250
#     gi_list = np.zeros(m_max + 1, dtype=np.float128)
# #    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)
#     gi_list[0] = np.float128(gammainc_mem(z, a))
#     for i in np.arange(1, m_max + 1):
#         gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.power(a, i - 1. + z) * np.exp(-a)/factorial_ary[i]
#     return gi_list

# def gammainc_lo_fct_ary(m_max, z, b):
#     #calculate naively an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
#     #subtraction of nearly equal numbers causes loss of precision around m ~ 100
#     gi_list = np.zeros(m_max + 1, dtype=np.float128)
#     gi_list[0] = np.float128(gammainc_mem(z, 0, b))
#     for i in np.arange(1, m_max + 1):
#         gi_list[i] = (i - 1. + z)*gi_list[i-1]/i - np.power(b, i - 1. + z) * np.exp(-b)/factorial_ary[i]
#     return gi_list
    
# def gammainc_up_fct_ary_log(np.int m_max, np.float z, np.float a):
#     #calculate in log-space an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
#     cdef np.ndarray gi_list = np.zeros(m_max + 1, dtype=np.float128)
# #    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)    
#     gi_list[0] = np.float128(gammainc_mem(z, a))

#     #cdef np.float i_t = 0.0
#     cdef np.ndarray i_array = np.arange(1., m_max + 1., dtype=np.float128)
#     # for i in np.arange(1., m_max + 1., dtype=np.float128):
#     cdef Py_ssize_t i 
#     for i in np.arange(1, m_max + 1):
#         gi_list[i] = (i_array[i-1] - 1. + z)*gi_list[i-1]/i + np.exp((i_array[i-1] - 1. + z)*np.log(a) - a - log_factorial_ary[i])
#         # gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.exp((i - 1. + z)*np.log(a) - a - log_factorial_ary[i])
#     return np.asarray(gi_list,dtype=np.float64)


# def gammainc_up_fct_ary_log_old(m_max, z, a):
#     #calculate in log-space an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
#     gi_list = np.zeros(m_max + 1, dtype=np.float128)
# #    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)    
#     gi_list[0] = np.float128(gammainc_mem(z, a))
#     for i in np.arange(1., m_max + 1., dtype=np.float128):
#         gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.exp((i - 1. + z)*np.log(a) - a - log_factorial_ary[i])
#     return gi_list

# # def gammainc_up_fct_ary_mid(m_max, z, a, b):
# #     #calculate in log-space an array of middle incomplete gamma function divided by factorial (gamma(z + m, a, b)/m!) for m in [0, m_max]
# #     gi_list = np.zeros(m_max + 1, dtype=np.float128)
# # #    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)    
# #     gi_list[0] = np.float128(gammainc_mem(z, a, b))
# #     for i in np.arange(1., m_max + 1., dtype=np.float128):
# #         gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.exp((i - 1. + z)*np.log(a) - a - log_factorial_ary[i])
# #     return gi_list
    
# def gammainc_lo_fct_ary_back(np.int m_max, np.float z, np.float  b):
#     #calculate in log-space an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
#     #use backwards recursion to avoid loss of precision from subtraction
#     cdef np.ndarray gi_list = np.zeros(m_max + 1,dtype=np.float128)
#     gi_list[m_max] = np.float128(gammainc_mem(z + m_max, 0, b) / mp.factorial(m_max))

#     cdef np.ndarray i_array = np.arange(m_max, 0, -1, dtype=np.float128)

#     cdef Py_ssize_t i 
#     for i in np.arange(m_max, 0, -1):
#         gi_list[i-1] = (gi_list[i] + np.exp((i_array[-i] - 1. + np.float128(z))*np.log(np.float128(b)) - np.float128(b) - log_factorial_ary[i]))*i_array[-i]/(i_array[-i] - 1. + np.float128(z))
#     return np.asarray(gi_list,dtype=np.float64)


# def gammainc_lo_fct_ary_back_old(m_max,z,  b):
#     #calculate in log-space an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
#     #use backwards recursion to avoid loss of precision from subtraction
#     gi_list = np.zeros(m_max + 1)
#     gi_list[m_max] = np.float128(gammainc_mem(z + m_max, 0, b) / mp.factorial(m_max))
#     b = np.float128(b)
#     for i in np.arange(m_max, 0, -1, dtype=np.float128):
#         gi_list[i-1] = (gi_list[i] + np.exp((i - 1. + z)*np.log(b) - b - log_factorial_ary[i]))*i/(i - 1. + z)
#     return gi_list


# import numpy as np
# import math
# cimport numpy as np
# cimport cython
# from mpmath import mp

# cdef extern from "math.h":
#     double log(double x) nogil
#     double exp(double x) nogil

# sp_max = 5000 

# def log_factorial(k):
#     return np.sum(np.log(np.arange(1., k + 1., dtype=np.float)))

# cdef double[::1] factorial_ary = np.array([np.float(mp.factorial(m)) for m in np.arange(sp_max)])
# cdef double[::1] log_factorial_ary = np.vectorize(log_factorial)(np.arange(sp_max))

# gamma_mem = mp.gamma
# gammainc_mem = mp.gammainc

# def gammainc_up_fct_ary_log(m_max, z, a):
#     return gammainc_up_fct_ary_log_impl(m_max, z, a)

# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# cdef double[::1] gammainc_up_fct_ary_log_impl(int m_max, double z, double a):
#     cdef double[::1] gi_list = np.zeros(m_max + 1, dtype=np.float)
#     gi_list[0] = gammainc_mem(z, a)
#     cdef Py_ssize_t i
#     for i in range(1, m_max + 1):
#         t0 = (i - 1. + z)
#         t1 = (i - 1. + z)*log(a) - a
#         gi_list[i] = t0*gi_list[i-1]/i + exp(t1 - log_factorial_ary[i])
#         if math.isnan(gi_list[i]) or math.isinf(gi_list[i]) or gi_list[i]==0.0:
#             gi_list[i] = 10**-300.
#     if math.isnan(gi_list[0]) or math.isinf(gi_list[0]) or gi_list[0]==0.0:
#             gi_list[0] = 10**-300.
#     return gi_list


# def gammainc_up_fct_ary_log_old(m_max, z, a):
#     #calculate in log-space an array of upper incomplete gamma function divided by factorial (gamma(z + m, a, inf)/m!) for m in [0, m_max]
#     gi_list = np.zeros(m_max + 1, dtype=np.float128)
# #    gi_list[0] = np.float128(gamma_mem(z) - gammainc_mem(z, 0, a)) #faster than gammainc_mem(z, a)    
#     gi_list[0] = np.float128(gammainc_mem(z, a))
#     for i in range(1, m_max + 1):
#         gi_list[i] = (i - 1. + z)*gi_list[i-1]/i + np.exp((i - 1. + z)*np.log(a) - a - log_factorial_ary[i])
#     return gi_list


# def gammainc_lo_fct_ary_back(m_max, z, b):
#     return gammainc_lo_fct_ary_back_impl(m_max, z, b)

# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# cdef double[::1] gammainc_lo_fct_ary_back_impl(int m_max, double z, double b):
#     #calculate in log-space an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
#     #use backwards recursion to avoid loss of precision from subtraction
#     cdef double[::1] gi_list = np.zeros(m_max + 1,dtype=np.float)
#     gi_list[m_max] = gammainc_mem(z + m_max, 0, b) / mp.factorial(m_max)

#     #cdef np.ndarray i_array = np.arange(m_max, 0, -1, dtype=np.float)

#     cdef Py_ssize_t i 
#     for i in range(m_max, 0, -1):
#         t0 = i/(i - 1. + z) 
#         t1 = exp(i-1+z)*log(b) -b
#         gi_list[i-1] = (gi_list[i] + t1 - log_factorial_ary[i])*t0
#         if math.isnan(gi_list[i-1]) or math.isinf(gi_list[i-1])or gi_list[i-1]==0.0:
#             gi_list[i-1] = 10**-300.
#     if math.isnan(gi_list[m_max]) or math.isinf(gi_list[m_max]) or gi_list[m_max]==0.0:
#             gi_list[m_max] = 10**-300.
#     return gi_list




















# ###################
# ###This is faster but not enough precision


# # import numpy as np
# # cimport numpy as np
# # cimport cython
# # from mpmath import mp

# # cdef extern from "math.h":
# #     long double log(long double x) nogil
# #     long double exp(long double x) nogil

# # sp_max = 5000 

# # def log_factorial(k):
# #     return np.sum(np.log(np.arange(1., k + 1., dtype=np.float128)))

# # factorial_ary = np.array([np.float128(mp.factorial(m)) for m in np.arange(sp_max)])
# # log_factorial_ary = np.vectorize(log_factorial)(np.arange(sp_max))

# # gamma_mem = mp.gamma
# # gammainc_mem = mp.gammainc

# # def gammainc_up_fct_ary_log(m_max,double z, double a):
# #     return gammainc_up_fct_ary_log_impl(m_max, z, a)

# # @cython.boundscheck(False)
# # @cython.wraparound(False)
# # @cython.cdivision(True)
# # cdef np.ndarray gammainc_up_fct_ary_log_impl(int m_max, long double z_t, long double a_t):
# #     cdef long double z = <long double>z_t
# #     cdef long double a = <long double>a_t
# #     cdef long double[::1] gi_list = np.zeros(m_max + 1, dtype=np.float128)
# #     gi_list[0] = gammainc_mem(z, a)
# #     cdef long double t0
# #     cdef long double t1
# #     cdef Py_ssize_t i
# #     for i in range(1, m_max + 1):
# #         t0 = (i - 1. + z)
# #         t1 = (i - 1. + z)*log(a) - a
# #         gi_list[i] = t0*gi_list[i-1]/i + exp(t1 - log_factorial_ary[i])
# #     return np.asarray(gi_list,np.float64)


# # def gammainc_lo_fct_ary_back(m_max, double z, double b):
# #     return gammainc_lo_fct_ary_back_impl(m_max, z, b)

# # @cython.boundscheck(False)
# # @cython.wraparound(False)
# # @cython.cdivision(True)
# # cdef np.ndarray gammainc_lo_fct_ary_back_impl(int m_max, double z_t, double b_t):
# #     #calculate in log-space an array of lower incomplete gamma function divided by factorial (gamma(z + m, 0, b)/m!) for m in [0, m_max]
# #     #use backwards recursion to avoid loss of precision from subtraction
# #     cdef long double z = <long double>z_t
# #     cdef long double b = <long double>b_t
# #     cdef long double[::1] gi_list = np.zeros(m_max + 1,dtype=np.float128)
# #     cdef long double t0
# #     cdef long double t1
# #     gi_list[m_max] = gammainc_mem(z + m_max, 0, b) / mp.factorial(m_max)

# #     #cdef np.ndarray i_array = np.arange(m_max, 0, -1, dtype=np.float)

# #     cdef Py_ssize_t i 
# #     for i in range(m_max, 0, -1):
# #         t0 = i/(i - 1. + z) 
# #         t1 = exp(i-1+z)*log(b) -b
# #         gi_list[i-1] = (gi_list[i] + t1 - log_factorial_ary[i])*t0
# #     return np.asarray(gi_list,np.float64)