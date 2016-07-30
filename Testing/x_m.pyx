#------------------------------------------------------------------------------
# x_m and x_m_sum arrays for arbitary number of breaks
#------------------------------------------------------------------------------
##
## Calculates x_m and x_m_sum arrays for arbitary number of breaks.
## Uses dedicated functions for 1, 2 and 3 breaks as well as n > 3 breaks.
##
#------------------------------------------------------------------------------

import numpy as np
cimport numpy as np
from mpmath import mp
cimport cython

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef extern from "math.h":
    double log(double x) nogil
    double exp(double x) nogil
    double pow(double x, double y) nogil

import incgamma_fct as ig

def return_xs(double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data):
    """
    Returns arrays of x_m and x_m_sum
    """
    cdef int n_break = int((len(theta) - 2)/2) # TODO: handle exception in length of theta etc

    if n_break == 1:
        return return_xs_1break(theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data)
    if n_break == 2:
        return return_xs_2break(theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data)
    if n_break == 3:
        return return_xs_nbreak(theta, n_break, f_ary, df_rho_div_f_ary, PS_dist_compressed, data)
    else:
        return return_xs_nbreak(theta, n_break, f_ary, df_rho_div_f_ary, PS_dist_compressed, data)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
def return_xs_1break(double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data):

    cdef int k_max = int(max(data) + 1)

    cdef double A = float(theta[0])
    cdef double n1 = float(theta[1])
    cdef double n2 = float(theta[2])
    cdef double Sb = float(theta[3])
    
    cdef int npixROI = len(PS_dist_compressed)

    cdef double[:,::1] x_m_ary = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)

    cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)

    cdef double f2, df_rho_div_f2

    cdef double pref1_x_m_ary, pref2_x_m_ary

    cdef double x_m_sum_f, x_m_ary_f

    cdef Py_ssize_t f_index, p, k, n

    for f_index in range(len(f_ary)):
        
        f2 = float(f_ary[f_index])
        df_rho_div_f2 = df_rho_div_f_ary[f_index]
        
        g1_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n1, Sb * f2) 
        g2_ary_f = ig.incgamma_lo_fct_ary(k_max, 1. - n2, Sb * f2)
        
        pref1_x_m_ary =  pow(Sb * f2, n1)
        pref2_x_m_ary = pow(Sb * f2, n2)

        for p in range(npixROI):

            x_m_sum_f = A*Sb*f2*(1/(n1-1)+1/(1-n2)) * PS_dist_compressed[p]
            x_m_sum[p] += df_rho_div_f2*x_m_sum_f

            for k in range(data[p]+1):

                x_m_ary_f = A * (pref1_x_m_ary  * g1_ary_f[k] + pref2_x_m_ary * g2_ary_f[k] ) * PS_dist_compressed[p]
                x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f

    x_m_sum = np.asarray(x_m_sum) - np.asarray(x_m_ary)[:,0] # x_m_sum -= x_m_ary[:,0]

    return np.asarray(x_m_ary), np.asarray(x_m_sum)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
def return_xs_2break(double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data):

    cdef int k_max = int(np.max(data) + 1)

    cdef float A = float(theta[0])
    cdef float n1 = float(theta[1])
    cdef float n2 = float(theta[2])
    cdef float n3 = float(theta[3])
    cdef float Sb1 = float(theta[4])
    cdef float Sb2 = float(theta[5])

    cdef int npixROI = len(PS_dist_compressed)

    cdef double[:,::1] x_m_ary = np.zeros(shape=(npixROI,k_max+1))
    cdef double[::1] x_m_sum = np.zeros(npixROI)

    cdef double[::1] g0_ary_f = np.zeros(k_max + 1)
    cdef double[::1] g1_ary_f = np.zeros(k_max + 1)
    cdef double[::1] g2_ary_f = np.zeros(k_max + 1)

    cdef double f2, df_rho_div_f2

    cdef double pref0_x_m_ary, pref1_x_m_ary, pref2_x_m_ary

    cdef double first0_x_m_sum_ary, first1_x_m_sum_ary, first2_x_m_sum_ary
    cdef double second0_x_m_sum_ary, second1_x_m_sum_ary, second2_x_m_sum_ary

    cdef double x_m_sum_f, x_m_ary_f

    cdef Py_ssize_t f_index, p, k, n

    first0_x_m_sum_ary = 1/(1-n1)
    first1_x_m_sum_ary = 1/(1-n2)
    first2_x_m_sum_ary = 1/(1-n3)

    second0_x_m_sum_ary = -1.0
    second1_x_m_sum_ary = (1 - pow(Sb2/Sb1, 1-n2))
    second2_x_m_sum_ary = pow(Sb2/Sb1, 1-n2)


    for f_index in range(len(f_ary)):

        f2 = float(f_ary[f_index])
        df_rho_div_f2 = df_rho_div_f_ary[f_index]
        
        g0_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n1, Sb1 * f2)
        g1_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n2, Sb2 * f2) - ig.incgamma_up_fct_ary(k_max, 1. - n2, Sb1 * f2)
        g2_ary_f = ig.incgamma_lo_fct_ary(k_max, 1. - n3, Sb2 * f2)

        pref0_x_m_ary = pow(Sb1 * f2, n1) 
        pref1_x_m_ary =  pref0_x_m_ary * pow(Sb1 * f2, n2 - n1)
        pref2_x_m_ary = pref1_x_m_ary * pow(Sb2 * f2, n3 - n2)

        for p in range(npixROI):
            x_m_sum_f = A * (Sb1 *f2) * (first0_x_m_sum_ary*second0_x_m_sum_ary+ first1_x_m_sum_ary*second1_x_m_sum_ary + first2_x_m_sum_ary*second2_x_m_sum_ary) *PS_dist_compressed[p]
            x_m_sum[p] += df_rho_div_f2*x_m_sum_f

            for k in range(data[p]+1):
                x_m_ary_f = A * (pref0_x_m_ary * g0_ary_f[k] + pref1_x_m_ary  * g1_ary_f[k] + pref2_x_m_ary * g2_ary_f[k]) * PS_dist_compressed[p]
                x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f

    x_m_sum = np.asarray(x_m_sum) - np.asarray(x_m_ary)[:,0] # x_m_sum -= x_m_ary[:,0]

    return np.asarray(x_m_ary), np.asarray(x_m_sum)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
def return_xs_3break(double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data):

    cdef int k_max = int(np.max(data) + 1)

    cdef float A = float(theta[0])
    cdef float n1 = float(theta[1])
    cdef float n2 = float(theta[2])
    cdef float n3 = float(theta[3])
    cdef float n4 = float(theta[4])
    cdef float Sb1 = float(theta[5])
    cdef float Sb2 = float(theta[6])
    cdef float Sb3 = float(theta[7])

    cdef int npixROI = len(PS_dist_compressed)

    cdef double[:,::1] x_m_ary = np.zeros((npixROI,k_max + 1))
    cdef double[::1] x_m_sum = np.zeros(npixROI)

    cdef double[::1] g0_ary_f = np.zeros(k_max + 1)
    cdef double[::1] g1_ary_f = np.zeros(k_max + 1)
    cdef double[::1] g2_ary_f = np.zeros(k_max + 1)
    cdef double[::1] g3_ary_f = np.zeros(k_max + 1)

    cdef double f2, df_rho_div_f2

    cdef double pref0_x_m_ary, pref1_x_m_ary, pref2_x_m_ary, pref3_x_m_ary

    cdef double first0_x_m_sum_ary, first1_x_m_sum_ary, first2_x_m_sum_ary, first3_x_m_sum_ary
    cdef double second0_x_m_sum_ary, second1_x_m_sum_ary, second2_x_m_sum_ary, second3_x_m_sum_ary

    cdef double x_m_sum_f, x_m_ary_f

    cdef Py_ssize_t f_index, p, k, n

    first0_x_m_sum_ary = 1/(1-n1)
    first1_x_m_sum_ary = 1/(1-n2)
    first2_x_m_sum_ary = 1/(1-n3)
    first3_x_m_sum_ary = 1/(1-n4)

    second0_x_m_sum_ary = -1.0
    second1_x_m_sum_ary = (1 - pow(Sb2/Sb1, 1-n2))
    second2_x_m_sum_ary = pow(Sb2/Sb1, 1-n2) * (1 - pow(Sb3 / Sb2, 1- n3))
    second3_x_m_sum_ary = pow(Sb2/Sb1, 1-n2) * pow(Sb3 / Sb2, 1 - n3)

    for f_index in range(len(f_ary)):
        f2 = float(f_ary[f_index])
        df_rho_div_f2 = df_rho_div_f_ary[f_index]

        g0_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n1, Sb1 * f2) 
        g1_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n2, Sb2 * f2) - ig.incgamma_up_fct_ary(k_max, 1. - n2, Sb1 * f2)
        g2_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n3, Sb3 * f2) - ig.incgamma_up_fct_ary(k_max, 1. - n3, Sb2 * f2)
        g3_ary_f = ig.incgamma_lo_fct_ary(k_max, 1. - n4, Sb3 * f2)

        pref0_x_m_ary = pow(Sb1 * f2, n1)
        pref1_x_m_ary = pref0_x_m_ary * pow(Sb1 * f2, n2 - n1)
        pref2_x_m_ary = pref1_x_m_ary * pow(Sb2 * f2, n3 - n2)
        pref3_x_m_ary = pref2_x_m_ary * pow(Sb3 * f2,  n4 - n3)

        for p in range(npixROI):
            x_m_sum_f = (A * Sb1 * f2) * (first0_x_m_sum_ary*second0_x_m_sum_ary+ first1_x_m_sum_ary*second1_x_m_sum_ary + first2_x_m_sum_ary*second2_x_m_sum_ary + first3_x_m_sum_ary*second3_x_m_sum_ary)  * PS_dist_compressed[p]
            x_m_sum[p] += df_rho_div_f2*x_m_sum_f

            for k in range(data[p]+1):
                x_m_ary_f = A * (pref0_x_m_ary * g0_ary_f[k] + pref1_x_m_ary  * g1_ary_f[k] + pref2_x_m_ary * g2_ary_f[k]+ pref3_x_m_ary * g3_ary_f[k]) * PS_dist_compressed[p]
                x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f

    x_m_sum = np.asarray(x_m_sum) - np.asarray(x_m_ary)[:,0] # x_m_sum -= x_m_ary[:,0]

    return np.asarray(x_m_ary), np.asarray(x_m_sum)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
def return_xs_nbreak(double[::1] theta, int n_break, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data):

    cdef int k_max = int(np.max(data) + 1)

    cdef float A = float(theta[0])
    cdef double[::1] n_ary = theta[1:n_break + 2]
    cdef double[::1] Sb_ary = theta[n_break + 2: 2*n_break + 2]
    
    cdef int npixROI = len(PS_dist_compressed)

    cdef double[:,::1] x_m_ary = np.zeros((npixROI,k_max + 1))
    cdef double[::1] x_m_sum = np.zeros(npixROI)

    cdef double[:,::1] g_ary_f_ary = np.zeros((n_break + 1, k_max + 1))
    cdef double[::1] pref_x_m_ary = np.zeros(n_break + 1)

    cdef double[::1] first_x_m_sum_ary = np.zeros(n_break + 1)
    cdef double[::1] second_x_m_sum_ary = np.zeros(n_break + 1)

    cdef double[::1] gamma_ary = np.zeros(k_max + 1)

    cdef double f2, df_rho_div_f2

    cdef double x_m_sum_f, x_m_ary_f

    cdef Py_ssize_t f_index, p, k, n, m, i

    for f_index in range(len(f_ary)):

        f2 = float(f_ary[f_index])
        df_rho_div_f2 = df_rho_div_f_ary[f_index]

        gamma_ary = ig.incgamma_up_fct_ary(k_max, 1. - n_ary[0], Sb_ary[0] * f2)
        
        for m in range(k_max + 1):
            g_ary_f_ary[0][m] = gamma_ary[m]

        
        for i in range(1, n_break):
            gamma_ary = ig.incgamma_up_fct_ary(k_max, 1. - n_ary[i], Sb_ary[i] * f2) - ig.incgamma_up_fct_ary(k_max, 1. - n_ary[i], Sb_ary[i-1] * f2)
            for m in range(k_max + 1):
                g_ary_f_ary[i][m] = gamma_ary[m]

        gamma_ary = ig.incgamma_lo_fct_ary(k_max, 1. - n_ary[n_break], Sb_ary[n_break-1] * f2) 
        for m in range(k_max + 1):
            g_ary_f_ary[n_break][m] = gamma_ary[m]

        pref_x_m_ary[0] = pow(Sb_ary[0] * f2, n_ary[0]) * 1.
        for k in range(1, n_break):
            pref_x_m_ary[k] = pow(Sb_ary[0] * f2, n_ary[0]) * np.prod([pow(Sb_ary[j-1] * f2 , n_ary[j] - n_ary[j-1]) for j in range(1, k + 1)]) # make sure np.prod is okay precision-wise!
        pref_x_m_ary[n_break] = pow(Sb_ary[0] * f2, n_ary[0]) * np.prod([pow(Sb_ary[j-1] * f2, n_ary[j] - n_ary[j-1]) for j in range(1, n_break + 1)])

        second_x_m_sum_ary = np.asarray([1/(1 - n_ary[i]) for i in range(n_break+1)])

        first_x_m_sum_ary[0] = -1. 
        for i in range(1, n_break):
            first_x_m_sum_ary[i] = -1*np.prod([pow(Sb_ary[j] / Sb_ary[j - 1], - n_ary[j]) for j in range(1, i)]) * pow(Sb_ary[i-1], n_ary[i]) * (-pow(Sb_ary[i-1], 1-n_ary[i])+pow(Sb_ary[i], 1-n_ary[i])) / Sb_ary[0] 
        first_x_m_sum_ary[n_break] = np.prod([pow(Sb_ary[j] / Sb_ary[j - 1], - n_ary[j]) for j in range(1, n_break)]) * Sb_ary[n_break-1] / Sb_ary[0]

        for p in range(npixROI):
            
            x_m_sum_f = A * Sb_ary[0] * f2 * np.dot(first_x_m_sum_ary, second_x_m_sum_ary) * PS_dist_compressed[p]
            x_m_sum[p] += df_rho_div_f2 * x_m_sum_f
            
            for k in range(data[p]+1):

                x_m_ary_f = A * np.dot(pref_x_m_ary, np.transpose(g_ary_f_ary)[k]) * PS_dist_compressed[p]
                x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f

    x_m_sum = np.asarray(x_m_sum) - np.asarray(x_m_ary)[:,0] # x_m_sum -= x_m_ary[:,0]

    return np.asarray(x_m_ary), np.asarray(x_m_sum)
