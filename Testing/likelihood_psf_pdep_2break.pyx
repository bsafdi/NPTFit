import numpy as np
cimport numpy as np
cimport cython

import pulsars.special as spc

import logging
logger = logging.getLogger(__name__)

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef extern from "math.h":
    double log(double x) nogil
    double exp(double x) nogil
    double pow(double x, double y) nogil

import time

def log_nu_k_ary_PSF_exact_1_PS_2break(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data):
    return log_nu_k_ary_PSF_exact_1_PS_2break_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_nu_k_ary_PSF_exact_1_PS_2break_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data):

    cdef int k_max = np.max(data) + 1

    cdef double A = np.float(theta[0])
    cdef double n1 = np.float(theta[1])
    cdef double n2 = np.float(theta[2])
    cdef double n3 = np.float(theta[3])
    cdef double Sb1 = np.float(theta[4])
    cdef double Sb2 = np.float(theta[5])

    cdef int npixROI = len(xbg_PSF_compressed)

    cdef double f2 = 0.0
    cdef double df_rho_div_f2 = 0.0


    cdef double[:,::1] x_m_ary = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
    cdef double x_m_ary_f
    cdef double x_m_sum_f

    cdef double[::1] g0_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)

    cdef Py_ssize_t f_index, p, k, n


    #calculations for PS

    cdef int do_half = 0

    cdef double term1 = 0.0
    cdef double term2 = 0.0
    cdef double second_2_a = 0.0
    cdef double second_2_b = 0.0
    cdef double second_2_c = 0.0
    cdef double second_2_d = 0.0
    cdef double second_1_alpha = 0.0
    cdef double second_1_a = 0.0
    cdef double second_1_b = 0.0
    cdef double second_1_beta = 0.0
    cdef double second_1_c = 0.0
    cdef double second_1_d = 0.0



    for f_index in range(len(f_ary)):
        f2 = f_ary[f_index]
        df_rho_div_f2 = df_rho_div_f_ary[f_index]
        g0_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb2 * f2) 
        g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n2, Sb1 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n2, Sb2 * f2)
        g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n3, Sb1 * f2)
        term1 = (A * Sb1 * f2) \
                             * (1./(n2-1.) + 1./(1.-n3) - pow(Sb1 / Sb2, n2-1.)/(n2-1.) \
                                - (pow(Sb1 * f2, n2-1.) * g1_ary_f[0] + pow(Sb1 * f2, n3-1.) * g2_ary_f[0])   + pow( Sb1 / Sb2, n2 - 1.)*( 1./(n1-1.) - pow(Sb2*f2, n1 - 1.)*g0_ary_f[0] )  )
        second_1_alpha = A * pow(Sb2 * f2, n1)* pow(Sb1 / Sb2 , n2)
        second_1_a =  A  * pow(Sb1 * f2, n2)
        second_1_b = A * pow(Sb1 * f2, n3)

        for p in range(npixROI):
            x_m_sum_f = term1 * PS_dist_compressed[p]
            x_m_sum[p] += df_rho_div_f2*x_m_sum_f

            second_1_beta = second_1_alpha * PS_dist_compressed[p]
            second_1_c = second_1_a * PS_dist_compressed[p]
            second_1_d = second_1_b * PS_dist_compressed[p]
            for k in range(data[p]+1):
                x_m_ary_f = second_1_beta * g0_ary_f[k] + second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k] 
                x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f


    cdef double[::1] nu_ary = np.zeros(k_max + 1, dtype=DTYPE)

    cdef double f0_ary
    cdef double f1_ary

    cdef double[:] nu_mat = np.zeros((k_max+1), dtype=DTYPE)


    cdef double ll = 0.

    for p in range(npixROI):
        f0_ary = -(xbg_PSF_compressed[p] + x_m_sum[p])
        f1_ary = (xbg_PSF_compressed[p] + x_m_ary[p,1])
        nu_mat[0] = exp(f0_ary)
        nu_mat[1] = nu_mat[0] * f1_ary

        for k in range(2,data[p]+1):
            nu_mat[k] = 0.0
            for n in range(0, k - 1):
                nu_mat[k] += (k-n)/ float(k) * x_m_ary[p,k-n] * nu_mat[n]
            nu_mat[k] += f1_ary * nu_mat[k-1] / float(k)

     #   print 'At data[p] = ', data[p], 'nu_mat[data[p]] = ', nu_mat[data[p]]

        if nu_mat[data[p]] > 0:
            ll+=log( nu_mat[data[p]])
        else:
            ll+= -10.1**10.

    # if math.isnan(ll) or math.isinf(ll):
    #     ll = -10.1**10.

    return ll

