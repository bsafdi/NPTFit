import numpy as np
cimport numpy as np
# import healpy as hp
# import mpmath as mp
# import math
# import matplotlib.pyplot as plt
cimport cython
#from libcpp cimport bool

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

def sum_arrays(double[::1] array_1, double[::1] array_2):
    return sum_arrays_int(array_1, array_2)

def sub_arrays(double[::1] array_1, double[::1] array_2):
    return sub_arrays_int(array_1, array_2)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double[::1] sum_arrays_int(double[::1] array_1, double[::1] array_2):
    cdef int n = len(array_1)
    cdef double[::1] sum_array =  np.zeros(n,dtype=DTYPE)
    cdef Py_ssize_t i
    for i in range(n):
        sum_array[i] = array_1[i] + array_2[i]

    return sum_array

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double[::1] sub_arrays_int(double[::1] array_1, double[::1] array_2):
    cdef int n = len(array_1)
    cdef double[::1] sum_array =  np.zeros(n,dtype=DTYPE)
    cdef Py_ssize_t i
    for i in range(n):
        sum_array[i] = array_1[i] - array_2[i]

    return sum_array

def log_nu_k_ary_PSF_exact_3_PS(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, double[::1] PS_dist_compressed2,double[::1] PS_dist_compressed3, int[::1] data, double Sc = 100000.0, double[::1] x_m_sum2_t = np.array([-10. ,-10.],dtype=DTYPE), double[:,::1] x_m_ary2_t = np.zeros((2,2), dtype=DTYPE),double[::1] x_m_sum3_t = np.array([-10. ,-10.],dtype=DTYPE), double[:,::1] x_m_ary3_t = np.zeros((2,2), dtype=DTYPE)):
    return log_nu_k_ary_PSF_exact_3_PS_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, PS_dist_compressed2,PS_dist_compressed3, data, Sc, x_m_sum2_t , x_m_ary2_t,x_m_sum3_t , x_m_ary3_t)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_nu_k_ary_PSF_exact_3_PS_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, double[::1] PS_dist_compressed2,double[::1] PS_dist_compressed3, int[::1] data, double Sc , double[::1] x_m_sum2_t , double[:,::1] x_m_ary2_t,double[::1] x_m_sum3_t , double[:,::1] x_m_ary3_t ):
    #t_0 = time.time()
    # print 'k_max is :', k_max
    # print 'Number of pixels:', len(xbg_PSF_compressed)
    # print 'theta: ', theta

    #cdef np.float t_i = time.time()

    #print 'here theta: ', np.asarray(theta)

    cdef int k_max = np.max(data) + 1

    cdef int len_theta = len(theta) 

    cdef double A = np.float(theta[0])
    cdef double n1 = np.float(theta[1])
    cdef double n2 = np.float(theta[2])
    cdef double Sb = np.float(theta[3])
    
    cdef double A2, n12, n22, Sb2, A3, n13, n23, Sb3
    if len_theta > 4:
        A2 = np.float(theta[4])
        n12 = np.float(theta[5])
        n22 = np.float(theta[6])
        Sb2 = np.float(theta[7])

    if len_theta > 8:
        A3 = np.float(theta[8])
        n13 = np.float(theta[9])
        n23 = np.float(theta[10])
        Sb3 = np.float(theta[11])

    cdef int npixROI = len(xbg_PSF_compressed)

    cdef double[:,::1] x_m_ary2 = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)
    cdef double[::1] g1_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g2_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double x_m_ary_f2 = 0.0
    cdef double x_m_sum_f2 = 0.0
    cdef double f2 = 0.0
    cdef double df_rho_div_f2 = 0.0
    
    
    cdef double[:,::1] x_m_ary = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
    cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double x_m_ary_f = 0.0
    cdef double x_m_sum_f = 0.0

    cdef double[:,::1] x_m_ary3 = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum3 = np.zeros(npixROI, dtype=DTYPE)
    cdef double[::1] g1_ary_f3 = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g2_ary_f3 = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double x_m_ary_f3 = 0.0
    cdef double x_m_sum_f3 = 0.0


    cdef Py_ssize_t f_index, p, k, n


    #calculations for PS
    
    cdef int do_2 = 0
    cdef int do_3 = 0

    cdef double term1 = 0.0
    cdef double term2 = 0.0
    cdef double term3 = 0.0
    cdef double second_3_a = 0.0
    cdef double second_3_b = 0.0
    cdef double second_3_c = 0.0
    cdef double second_3_d = 0.0
    cdef double second_2_a = 0.0
    cdef double second_2_b = 0.0
    cdef double second_2_c = 0.0
    cdef double second_2_d = 0.0
    cdef double second_1_a = 0.0
    cdef double second_1_b = 0.0
    cdef double second_1_c = 0.0
    cdef double second_1_d = 0.0

    if len_theta <= 4:
        x_m_ary2 = x_m_ary2_t
        x_m_sum2 = x_m_sum2_t
        do_2 = 1
    else:
        x_m_ary2 = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
        x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)

    if len_theta <= 8:
        x_m_ary3 = x_m_ary3_t
        x_m_sum3 = x_m_sum3_t
        do_3 = 1
    else:
        x_m_ary3 = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
        x_m_sum3 = np.zeros(npixROI, dtype=DTYPE)
    
    if do_2 and do_3:
        for f_index in range(len(f_ary)):
            f2 = f_ary[f_index]
            df_rho_div_f2 = df_rho_div_f_ary[f_index]
            g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
            g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
            term1 = (A * Sb * f2) \
                                 * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
                                    - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
            second_1_a =  A  * pow(Sb * f2, n1)
            second_1_b = A * pow(Sb * f2, n2)
            for p in range(npixROI):
                x_m_sum_f = term1 * PS_dist_compressed[p]
                x_m_sum[p] += df_rho_div_f2*x_m_sum_f

                second_1_c = second_1_a * PS_dist_compressed[p]
                second_1_d = second_1_b * PS_dist_compressed[p]
                for k in range(data[p]+1):   #####take over here!!!
                    x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
                    x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f
    elif do_3:
        for f_index in range(len(f_ary)):
            f2 = f_ary[f_index]
            df_rho_div_f2 = df_rho_div_f_ary[f_index]
            g1_ary_f2 = spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sb2 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sc * f2)
            g2_ary_f2 = spc.gammainc_lo_fct_ary_back(k_max, 1. - n22, Sb2 * f2)
            g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
            g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
            term2 = (A2  * Sb2 * f2) \
                                 * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
                                    - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))

            term1 = (A * Sb * f2) \
                                 * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
                                    - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
            second_2_a =  A2  * pow(Sb2 * f2, n12)
            second_2_b = A2 * pow(Sb2 * f2, n22)

            second_1_a =  A  * pow(Sb * f2, n1)
            second_1_b = A * pow(Sb * f2, n2)

            for p in range(npixROI):
                x_m_sum_f2 = term2 * PS_dist_compressed2[p]
                x_m_sum2[p] += df_rho_div_f2*x_m_sum_f2
                x_m_sum_f = term1 * PS_dist_compressed[p]
                x_m_sum[p] += df_rho_div_f2*x_m_sum_f

                second_2_c = second_2_a * PS_dist_compressed2[p]
                second_2_d = second_2_b * PS_dist_compressed2[p]

                second_1_c = second_1_a * PS_dist_compressed[p]
                second_1_d = second_1_b * PS_dist_compressed[p]
                for k in range(data[p]+1):   #####take over here!!!
                    x_m_ary_f2 = second_2_c  * g1_ary_f2[k] + second_2_d * g2_ary_f2[k]
                    x_m_ary2[p,k] += df_rho_div_f2*x_m_ary_f2
                    x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
                    x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f
    else:
        for f_index in range(len(f_ary)):
            f2 = f_ary[f_index]
            df_rho_div_f2 = df_rho_div_f_ary[f_index]
            g1_ary_f3 = spc.gammainc_up_fct_ary_log(k_max, 1. - n13, Sb3 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n13, Sc * f2) 
            g2_ary_f3 = spc.gammainc_lo_fct_ary_back(k_max, 1. - n23, Sb3 * f2)
            g1_ary_f2 = spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sb2 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sc * f2) 
            g2_ary_f2 = spc.gammainc_lo_fct_ary_back(k_max, 1. - n22, Sb2 * f2)
            g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2) 
            g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
            term3 = (A3  * Sb3 * f2) \
                                 * (1./(n13-1.) + 1./(1.-n23) - pow(Sb3 / Sc, n13-1.)/(n13-1.) \
                                    - (pow(Sb3 * f2, n13-1.) * g1_ary_f3[0] + pow(Sb3 * f2, n23-1.) * g2_ary_f3[0]))
            term2 = (A2  * Sb2 * f2) \
                                 * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
                                    - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))

            term1 = (A * Sb * f2) \
                                 * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
                                    - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
            second_3_a =  A3  * pow(Sb3 * f2, n13)
            second_3_b = A3 * pow(Sb3 * f2, n23)

            second_2_a =  A2  * pow(Sb2 * f2, n12)
            second_2_b = A2 * pow(Sb2 * f2, n22)

            second_1_a =  A  * pow(Sb * f2, n1)
            second_1_b = A * pow(Sb * f2, n2)

            for p in range(npixROI):
                x_m_sum_f3 = term3 * PS_dist_compressed3[p]
                x_m_sum3[p] += df_rho_div_f2*x_m_sum_f3
                x_m_sum_f2 = term2 * PS_dist_compressed2[p]
                x_m_sum2[p] += df_rho_div_f2*x_m_sum_f2
                x_m_sum_f = term1 * PS_dist_compressed[p]
                x_m_sum[p] += df_rho_div_f2*x_m_sum_f

                second_3_c = second_3_a * PS_dist_compressed3[p]
                second_3_d = second_3_b * PS_dist_compressed3[p]

                second_2_c = second_2_a * PS_dist_compressed2[p]
                second_2_d = second_2_b * PS_dist_compressed2[p]

                second_1_c = second_1_a * PS_dist_compressed[p]
                second_1_d = second_1_b * PS_dist_compressed[p]
                for k in range(data[p]+1):   #####take over here!!!
                    x_m_ary_f3 = second_3_c  * g1_ary_f3[k] + second_3_d * g2_ary_f3[k]
                    x_m_ary3[p,k] += df_rho_div_f2*x_m_ary_f3
                    x_m_ary_f2 = second_2_c  * g1_ary_f2[k] + second_2_d * g2_ary_f2[k]
                    x_m_ary2[p,k] += df_rho_div_f2*x_m_ary_f2
                    x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
                    x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f

    
    cdef double[::1] x_m_ary_total = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double x_m_sum_total = 0.0


    cdef double[::1] nu_ary = np.zeros(k_max + 1, dtype=DTYPE)
    
    cdef double f0_ary = 0.0 # -(xbg_PSF_compressed + x_m_sum_total)
    cdef double f1_ary = 0.0 # (xbg_PSF_compressed + x_m_ary_total[1])

    cdef double[::1] nu_mat = np.zeros(k_max+1, dtype=DTYPE)

    #t0 = time.time()
    cdef double ll = 0.
    #cdef double temp = 0.
    for p in range(npixROI):
        x_m_sum_total = x_m_sum3[p]+x_m_sum2[p]+ x_m_sum[p]
        x_m_ary_total[0] = x_m_ary3[p,0]+x_m_ary2[p,0] + x_m_ary[p,0]
        x_m_ary_total[1] = x_m_ary3[p,1]+x_m_ary2[p,1] + x_m_ary[p,1]
        f0_ary = -(xbg_PSF_compressed[p] + x_m_sum_total)
        f1_ary = (xbg_PSF_compressed[p] + x_m_ary_total[1])
        nu_mat[0] = exp(f0_ary)
        nu_mat[1] = nu_mat[0] * f1_ary
        #print nu_mat[1,p]
        for k in range(2,data[p]+1):
            x_m_ary_total[k] = x_m_ary3[p,k]+x_m_ary2[p,k] + x_m_ary[p,k]
            nu_mat[k] = 0.0
            for n in range(0, k - 1):
                nu_mat[k] += (k-n)/ float(k) * x_m_ary_total[k-n] * nu_mat[n]
            nu_mat[k] += f1_ary * nu_mat[k-1] / float(k)

        if nu_mat[data[p]] > 0:
            ll+=log( nu_mat[data[p]])
        else:
            ll += -10.1**10.

    # if math.isnan(ll) ==True or math.isinf(ll) ==True:
    #     ll = -10.1**10.

    #print 'll = ', ll
    return ll

def log_nu_k_ary_PSF_exact_2_PS(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, double[::1] PS_dist_compressed2, int[::1] data, double Sc = 100000.0, double[::1] x_m_sum2_t = np.array([-10. ,-10.],dtype=DTYPE), double[:,::1] x_m_ary2_t = np.zeros((2,2), dtype=DTYPE)):
    return log_nu_k_ary_PSF_exact_2_PS_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, PS_dist_compressed2, data, Sc, x_m_sum2_t , x_m_ary2_t)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_nu_k_ary_PSF_exact_2_PS_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, double[::1] PS_dist_compressed2, int[::1] data, double Sc , double[::1] x_m_sum2_t , double[:,::1] x_m_ary2_t ):
    #t_0 = time.time()
    # print 'k_max is :', k_max
    # print 'Number of pixels:', len(xbg_PSF_compressed)
    # print 'theta: ', theta

    #cdef np.float t_i = time.time()

    cdef int len_theta = len(theta) 

    cdef int k_max = np.max(data) + 1

    cdef double A = np.float(theta[0])
    cdef double n1 = np.float(theta[1])
    cdef double n2 = np.float(theta[2])
    cdef double Sb = np.float(theta[3])
    
    cdef double A2, n12, n22, Sb2
    if len_theta > 4:
        A2 = np.float(theta[4])
        n12 = np.float(theta[5])
        n22 = np.float(theta[6])
        Sb2 = np.float(theta[7])

    cdef int npixROI = len(xbg_PSF_compressed)

    cdef double[:,::1] x_m_ary2 = np.zeros((npixROI, k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)
    cdef double[::1] g1_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g2_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double x_m_ary_f2 = 0.0
    cdef double x_m_sum_f2 = 0.0
    cdef double f2 = 0.0
    cdef double df_rho_div_f2 = 0.0
    
    
    cdef double[:,::1] x_m_ary = np.zeros((npixROI, k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
    cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double x_m_ary_f = 0.0
    cdef double x_m_sum_f = 0.0
    # cdef double f = 0.0
    # cdef double df_rho_div_f = 0.0
    

    #k_max = int(np.max(data))

    cdef Py_ssize_t f_index, p, k, n

    #cdef np.float t_0 = time.time()

    
    
    #calculations for PS
    
    cdef int do_half = 0

    cdef double term1 = 0.0
    cdef double term2 = 0.0
    cdef double second_2_a = 0.0
    cdef double second_2_b = 0.0
    cdef double second_2_c = 0.0
    cdef double second_2_d = 0.0
    cdef double second_1_a = 0.0
    cdef double second_1_b = 0.0
    cdef double second_1_c = 0.0
    cdef double second_1_d = 0.0


    if len_theta <= 4:
        x_m_ary2 = x_m_ary2_t
        x_m_sum2 = x_m_sum2_t
        do_half = 1
    else:
        x_m_ary2 = np.zeros((npixROI, k_max + 1), dtype=DTYPE)
        x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)
    
    if do_half:
        for f_index in range(len(f_ary)):
            f2 = f_ary[f_index]
            df_rho_div_f2 = df_rho_div_f_ary[f_index]
            g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
            g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
            term1 = (A * Sb * f2) \
                                 * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
                                    - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
            second_1_a =  A  * pow(Sb * f2, n1)
            second_1_b = A * pow(Sb * f2, n2)
            for p in range(npixROI):
                x_m_sum_f = term1 * PS_dist_compressed[p]
                x_m_sum[p] += df_rho_div_f2*x_m_sum_f

                second_1_c = second_1_a * PS_dist_compressed[p]
                second_1_d = second_1_b * PS_dist_compressed[p]
                for k in range(data[p]+1):   #####take over here!!!
                    x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
                    x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f
    else:
        #t0 = time.time()
        for f_index in range(len(f_ary)):
            f2 = f_ary[f_index]
            df_rho_div_f2 = df_rho_div_f_ary[f_index]
            #t0 = time.time()
            g1_ary_f2 = spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sb2 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sc * f2)
            g2_ary_f2 = spc.gammainc_lo_fct_ary_back(k_max, 1. - n22, Sb2 * f2)
            g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
            g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
           # ta += time.time() - t0
            #print 'g1, g2 time: ', time.time()-t0
            term2 = (A2  * Sb2 * f2) \
                                 * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
                                    - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))

            term1 = (A * Sb * f2) \
                                 * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
                                    - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
            second_2_a =  A2  * pow(Sb2 * f2, n12)
            second_2_b = A2 * pow(Sb2 * f2, n22)

            second_1_a =  A  * pow(Sb * f2, n1)
            second_1_b = A * pow(Sb * f2, n2)

            for p in range(npixROI):
                #t0 = time.time()
                # x_m_sum_f2[p] = (A2 * PS_dist_compressed2[p] * Sb2 * f2) \
                #                  * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
                #                     - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))
                x_m_sum_f2 = term2 * PS_dist_compressed2[p]
                x_m_sum2[p] += df_rho_div_f2*x_m_sum_f2
                x_m_sum_f = term1 * PS_dist_compressed[p]
                x_m_sum[p] += df_rho_div_f2*x_m_sum_f
                #tb += time.time() - t0

                second_2_c = second_2_a * PS_dist_compressed2[p]
                second_2_d = second_2_b * PS_dist_compressed2[p]

                second_1_c = second_1_a * PS_dist_compressed[p]
                second_1_d = second_1_b * PS_dist_compressed[p]
                for k in range(data[p]+1):   #####take over here!!!
                    #t0 = time.time()
                    x_m_ary_f2 = second_2_c  * g1_ary_f2[k] + second_2_d * g2_ary_f2[k]
                    #print 'At k,p = ', k, p, 'xm_ary_f2 = ', x_m_ary_f2[k,p]
                    x_m_ary2[p,k] += df_rho_div_f2*x_m_ary_f2
                    x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
                    #print 'At k,p = ', k, p, 'xm_ary_f = ', x_m_ary_f[k,p]
                    x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f
    
    cdef double[::1] x_m_ary_total = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double x_m_sum_total = 0.0


    cdef double[::1] nu_ary = np.zeros(k_max + 1, dtype=DTYPE)
    
    cdef double f0_ary = 0.0 # -(xbg_PSF_compressed + x_m_sum_total)
    cdef double f1_ary = 0.0 # (xbg_PSF_compressed + x_m_ary_total[1])

    cdef double[::1] nu_mat = np.zeros(k_max+1, dtype=DTYPE)

    #t0 = time.time()
    cdef double ll = 0.
    #cdef double temp = 0.
    for p in range(npixROI):
        x_m_sum_total = x_m_sum2[p]+ x_m_sum[p]
        x_m_ary_total[0] = x_m_ary2[p,0] + x_m_ary[p,0]
        x_m_ary_total[1] = x_m_ary2[p,1] + x_m_ary[p,1]
        f0_ary = -(xbg_PSF_compressed[p] + x_m_sum_total)
        f1_ary = (xbg_PSF_compressed[p] + x_m_ary_total[1])
        nu_mat[0] = exp(f0_ary)
        nu_mat[1] = nu_mat[0] * f1_ary
        #print nu_mat[1,p]
        for k in range(2,data[p]+1):
            x_m_ary_total[k] = x_m_ary2[p,k] + x_m_ary[p,k]
            nu_mat[k] = 0.0
            for n in range(0, k - 1):
                nu_mat[k] += (k-n)/ float(k) * x_m_ary_total[k-n] * nu_mat[n]
            nu_mat[k] += f1_ary * nu_mat[k-1] / float(k)

        if nu_mat[data[p]] > 0:
            ll+=log( nu_mat[data[p]])
        else:
            ll += -10.1**10.



    # if math.isnan(ll) ==True or math.isinf(ll) ==True:
    #     ll = -10.1**10.


    return ll

# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# @cython.initializedcheck(False)
# cdef double log_nu_k_ary_PSF_exact_2_PS_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, double[::1] PS_dist_compressed2, int[::1] data, double Sc , double[::1] x_m_sum2_t , double[:,::1] x_m_ary2_t ):
#     #t_0 = time.time()
#     # print 'k_max is :', k_max
#     # print 'Number of pixels:', len(xbg_PSF_compressed)
#     # print 'theta: ', theta

#     #cdef np.float t_i = time.time()

#     cdef int k_max = np.max(data) + 1

#     cdef double A = np.float(theta[0])
#     cdef double n1 = np.float(theta[1])
#     cdef double n2 = np.float(theta[2])
#     cdef double Sb = np.float(theta[3])
    
#     cdef double A2 = np.float(theta[4])
#     cdef double n12 = np.float(theta[5])
#     cdef double n22 = np.float(theta[6])
#     cdef double Sb2 = np.float(theta[7])

#     cdef int npixROI = len(xbg_PSF_compressed)

#     cdef double[:,::1] x_m_ary2 = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)
#     cdef double[::1] g1_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[::1] g2_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[:,::1] x_m_ary_f2 = np.zeros((k_max + 1, npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum_f2 = np.zeros(npixROI, dtype=DTYPE)
#     cdef double f2 = 0.0
#     cdef double df_rho_div_f2 = 0.0
    
    
#     cdef double[:,::1] x_m_ary = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
#     cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[:,::1] x_m_ary_f = np.zeros((k_max + 1, npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum_f = np.zeros(npixROI, dtype=DTYPE)
#     # cdef double f = 0.0
#     # cdef double df_rho_div_f = 0.0
    

#     #k_max = int(np.max(data))

#     cdef Py_ssize_t f_index, p, k, n

#     #cdef np.float t_0 = time.time()

    
    
#     #calculations for PS
    
#     cdef int do_half = 0

#     cdef double term1 = 0.0
#     cdef double term2 = 0.0
#     cdef double second_2_a = 0.0
#     cdef double second_2_b = 0.0
#     cdef double second_2_c = 0.0
#     cdef double second_2_d = 0.0
#     cdef double second_1_a = 0.0
#     cdef double second_1_b = 0.0
#     cdef double second_1_c = 0.0
#     cdef double second_1_d = 0.0


#     if x_m_sum2_t[0] != -10.0:
#         x_m_ary2 = x_m_ary2_t
#         x_m_sum2 = x_m_sum2_t
#         do_half = 1
#     else:
#         x_m_ary2 = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#         x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)
    
#     if do_half:
#         for f_index in range(len(f_ary)):
#             f2 = f_ary[f_index]
#             df_rho_div_f2 = df_rho_div_f_ary[f_index]
#             g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
#             g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
#             for p in range(npixROI):
#                 x_m_sum_f[p] = (A * PS_dist_compressed[p] * Sb * f2) \
#                                  * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
#                                     - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
#                 x_m_sum[p] += df_rho_div_f2*x_m_sum_f[p]
#                 for k in range(data[p]+1):   #####take over here!!!
#                     x_m_ary_f[k,p] = A  * (pow(Sb * f2, n1)  * g1_ary_f[k] * PS_dist_compressed[p] + pow(Sb * f2, n2) * g2_ary_f[k]*PS_dist_compressed[p])
#                     x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f[k,p]
#     else:
#         #t0 = time.time()
#         for f_index in range(len(f_ary)):
#             f2 = f_ary[f_index]
#             df_rho_div_f2 = df_rho_div_f_ary[f_index]
#             #t0 = time.time()
#             g1_ary_f2 = spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sb2 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sc * f2)
#             g2_ary_f2 = spc.gammainc_lo_fct_ary_back(k_max, 1. - n22, Sb2 * f2)
#             g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
#             g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
#            # ta += time.time() - t0
#             #print 'g1, g2 time: ', time.time()-t0
#             term2 = (A2  * Sb2 * f2) \
#                                  * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
#                                     - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))

#             term1 = (A * Sb * f2) \
#                                  * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
#                                     - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
#             second_2_a =  A2  * pow(Sb2 * f2, n12)
#             second_2_b = A2 * pow(Sb2 * f2, n22)

#             second_1_a =  A  * pow(Sb * f2, n1)
#             second_1_b = A * pow(Sb * f2, n2)

#             for p in range(npixROI):
#                 #t0 = time.time()
#                 # x_m_sum_f2[p] = (A2 * PS_dist_compressed2[p] * Sb2 * f2) \
#                 #                  * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
#                 #                     - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))
#                 x_m_sum_f2[p] = term2 * PS_dist_compressed2[p]
#                 x_m_sum2[p] += df_rho_div_f2*x_m_sum_f2[p]
#                 x_m_sum_f[p] = term1 * PS_dist_compressed[p]
#                 x_m_sum[p] += df_rho_div_f2*x_m_sum_f[p]
#                 #tb += time.time() - t0

#                 second_2_c = second_2_a * PS_dist_compressed2[p]
#                 second_2_d = second_2_b * PS_dist_compressed2[p]

#                 second_1_c = second_1_a * PS_dist_compressed[p]
#                 second_1_d = second_1_b * PS_dist_compressed[p]
#                 for k in range(data[p]+1):   #####take over here!!!
#                     #t0 = time.time()
#                     x_m_ary_f2[k,p] = second_2_c  * g1_ary_f2[k] + second_2_d * g2_ary_f2[k]
#                     #print 'At k,p = ', k, p, 'xm_ary_f2 = ', x_m_ary_f2[k,p]
#                     x_m_ary2[k,p] += df_rho_div_f2*x_m_ary_f2[k,p]
#                     x_m_ary_f[k,p] = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
#                     #print 'At k,p = ', k, p, 'xm_ary_f = ', x_m_ary_f[k,p]
#                     x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f[k,p]
    
#     cdef double[:,::1] x_m_ary_total = np.zeros((k_max + 1, npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum_total = np.zeros(npixROI, dtype=DTYPE)


#     cdef double[::1] nu_ary = np.zeros(k_max + 1, dtype=DTYPE)
    
#     cdef double[::1] f0_ary = np.zeros(npixROI, dtype=DTYPE) # -(xbg_PSF_compressed + x_m_sum_total)
#     cdef double[::1] f1_ary = np.zeros(npixROI, dtype=DTYPE) # (xbg_PSF_compressed + x_m_ary_total[1])

#     cdef double[:,::1] nu_mat = np.zeros((k_max+1, npixROI), dtype=DTYPE)

#     #t0 = time.time()
#     cdef double ll = 0.
#     #cdef double temp = 0.
#     for p in range(npixROI):
#         x_m_sum_total[p] = x_m_sum2[p]+ x_m_sum[p]
#         x_m_ary_total[0,p] = x_m_ary2[0,p] + x_m_ary[0,p]
#         x_m_ary_total[1,p] = x_m_ary2[1,p] + x_m_ary[1,p]
#         f0_ary[p] = -(xbg_PSF_compressed[p] + x_m_sum_total[p])
#         f1_ary[p] = (xbg_PSF_compressed[p] + x_m_ary_total[1,p])
#         nu_mat[0,p] = exp(f0_ary[p])
#         nu_mat[1,p] = nu_mat[0,p] * f1_ary[p]
#         #print nu_mat[1,p]
#         for k in range(2,data[p]+1):
#             x_m_ary_total[k,p] = x_m_ary2[k,p] + x_m_ary[k,p]
#             for n in range(0, k - 1):
#                 nu_mat[k,p] += (k-n)/ float(k) * x_m_ary_total[k-n,p] * nu_mat[n,p]
#             nu_mat[k,p] += f1_ary[p] * nu_mat[k-1,p] / float(k)

#         ll+=log( nu_mat[data[p],p])


#     if math.isnan(ll) ==True or math.isinf(ll) ==True:
#         ll = -10.1**10.


#     return ll


# def log_nu_k_ary_PSF_exact_1_PS_edep(double[:,::1] xbg_PSF_compressed, double[:,::1] theta, double[:,::1] f_ary, double[:,::1] df_rho_div_f_ary, double[:,::1] PS_dist_compressed, int[:,::1] data, double Sc = 1000.0):
#     return log_nu_k_ary_PSF_exact_1_PS_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data, Sc)


def log_nu_k_ary_PSF_exact_1_PS(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data, double Sc = 100000.0):
    return log_nu_k_ary_PSF_exact_1_PS_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data, Sc)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_nu_k_ary_PSF_exact_1_PS_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data, double Sc ):

    cdef int k_max = np.max(data) + 1

    cdef double A = np.float(theta[0])
    cdef double n1 = np.float(theta[1])
    cdef double n2 = np.float(theta[2])
    cdef double Sb = np.float(theta[3])

    cdef int npixROI = len(xbg_PSF_compressed)

    cdef double f2 = 0.0
    cdef double df_rho_div_f2 = 0.0


    cdef double[:,::1] x_m_ary = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
    cdef double x_m_ary_f
    cdef double x_m_sum_f

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
    cdef double second_1_a = 0.0
    cdef double second_1_b = 0.0
    cdef double second_1_c = 0.0
    cdef double second_1_d = 0.0



    for f_index in range(len(f_ary)):
        f2 = f_ary[f_index]
        df_rho_div_f2 = df_rho_div_f_ary[f_index]
        g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
        g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
        term1 = (A * Sb * f2) \
                             * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
                                - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
        second_1_a =  A  * pow(Sb * f2, n1)
        second_1_b = A * pow(Sb * f2, n2)

        for p in range(npixROI):
            x_m_sum_f = term1 * PS_dist_compressed[p]
            x_m_sum[p] += df_rho_div_f2*x_m_sum_f

            second_1_c = second_1_a * PS_dist_compressed[p]
            second_1_d = second_1_b * PS_dist_compressed[p]
            for k in range(data[p]+1):
                x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k] 
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

#=====================================
#Just return x's 
def return_xs(double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data, double Sc = 100000.0):

    cdef int k_max = np.max(data) + 1

    cdef double A = np.float(theta[0])
    cdef double n1 = np.float(theta[1])
    cdef double n2 = np.float(theta[2])
    cdef double Sb = np.float(theta[3])

    cdef int npixROI = len(PS_dist_compressed)

    cdef double f2 = 0.0
    cdef double df_rho_div_f2 = 0.0


    cdef double[:,::1] x_m_ary = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
    cdef double x_m_ary_f
    cdef double x_m_sum_f

    cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)

    cdef Py_ssize_t f_index, p, k, n


    #calculations for PS

    cdef double term1 = 0.0
    cdef double second_1_a = 0.0
    cdef double second_1_b = 0.0
    cdef double second_1_c = 0.0
    cdef double second_1_d = 0.0



    for f_index in range(len(f_ary)):
        f2 = f_ary[f_index]
        df_rho_div_f2 = df_rho_div_f_ary[f_index]
        g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
        g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
        term1 = (A * Sb * f2) \
                             * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
                                - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
        second_1_a =  A  * pow(Sb * f2, n1)
        second_1_b = A * pow(Sb * f2, n2)

        for p in range(npixROI):
            x_m_sum_f = term1 * PS_dist_compressed[p]
            x_m_sum[p] += df_rho_div_f2*x_m_sum_f

            second_1_c = second_1_a * PS_dist_compressed[p]
            second_1_d = second_1_b * PS_dist_compressed[p]
            for k in range(data[p]+1):
                x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k] 
                x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f




    # cdef int k_max = int(np.max(data) + 1)

    # cdef np.float A = np.float(theta[0])
    # cdef np.float n1 = np.float(theta[1])
    # cdef np.float n2 = np.float(theta[2])
    # cdef np.float Sb = np.float(theta[3])

    # cdef np.int npixROI = len(PS_dist_compressed)
    
    # cdef np.ndarray x_m_ary = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    # cdef np.ndarray x_m_sum = np.zeros(npixROI, dtype=DTYPE)
    
    # cdef np.ndarray g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    # cdef np.ndarray g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    # cdef np.ndarray x_m_ary_f = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    # cdef np.ndarray x_m_sum_f = np.zeros(npixROI, dtype=DTYPE)
    # cdef np.float f = 0.0
    # cdef np.float df_rho_div_f = 0.0
    
    # cdef Py_ssize_t f_index2
    # for f_index2 in range(len(f_ary)):
    #     f = f_ary[f_index2]
    #     df_rho_div_f = df_rho_div_f_ary[f_index2]
    #     g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f)
    #     g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f)
    #     x_m_ary_f = A  * (np.power(Sb * f, n1) * np.outer(g1_ary_f, PS_dist_compressed) + np.power(Sb * f, n2) * np.outer(g2_ary_f, PS_dist_compressed))
    #     x_m_sum_f = (A * PS_dist_compressed * Sb * f) \
    #         * (1./(n1-1.) + 1./(1.-n2) - np.power(Sb / Sc, n1-1.)/(n1-1.)\
    #            - (np.power(Sb * f, n1-1.) * g1_ary_f[0] + np.power(Sb * f, n2-1.) * g2_ary_f[0]))


    #     x_m_ary += df_rho_div_f*x_m_ary_f
    #     x_m_sum += df_rho_div_f*x_m_sum_f

    return np.asarray(x_m_ary), np.asarray(x_m_sum)

# def log_nu_k_ary_PSF_exact_1_PS(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data, double Sc = 1000.0):
#     return log_nu_k_ary_PSF_exact_1_PS_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data, Sc)

# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# @cython.initializedcheck(False)
# cdef double log_nu_k_ary_PSF_exact_1_PS_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data, double Sc ):

#     cdef int k_max = np.max(data) + 1

#     cdef double A = np.float(theta[0])
#     cdef double n1 = np.float(theta[1])
#     cdef double n2 = np.float(theta[2])
#     cdef double Sb = np.float(theta[3])

#     cdef int npixROI = len(xbg_PSF_compressed)

#     cdef double f2 = 0.0
#     cdef double df_rho_div_f2 = 0.0
    
    
#     cdef double[:,::1] x_m_ary = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
#     cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[:,::1] x_m_ary_f = np.zeros((k_max + 1, npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum_f = np.zeros(npixROI, dtype=DTYPE)

#     cdef Py_ssize_t f_index, p, k, n
    
    
#     #calculations for PS
    
#     cdef int do_half = 0

#     cdef double term1 = 0.0
#     cdef double term2 = 0.0
#     cdef double second_2_a = 0.0
#     cdef double second_2_b = 0.0
#     cdef double second_2_c = 0.0
#     cdef double second_2_d = 0.0
#     cdef double second_1_a = 0.0
#     cdef double second_1_b = 0.0
#     cdef double second_1_c = 0.0
#     cdef double second_1_d = 0.0

    

#     for f_index in range(len(f_ary)):
#         f2 = f_ary[f_index]
#         df_rho_div_f2 = df_rho_div_f_ary[f_index]
#         g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
#         g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
#         term1 = (A * Sb * f2) \
#                              * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
#                                 - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
#         second_1_a =  A  * pow(Sb * f2, n1)
#         second_1_b = A * pow(Sb * f2, n2)

#         for p in range(npixROI):
#             x_m_sum_f[p] = term1 * PS_dist_compressed[p]
#             x_m_sum[p] += df_rho_div_f2*x_m_sum_f[p]

#             second_1_c = second_1_a * PS_dist_compressed[p]
#             second_1_d = second_1_b * PS_dist_compressed[p]
#             for k in range(data[p]+1):            
#                 x_m_ary_f[k,p] = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
#                 x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f[k,p]

                
#     cdef double[::1] nu_ary = np.zeros(k_max + 1, dtype=DTYPE)
    
#     cdef double[::1] f0_ary = np.zeros(npixROI, dtype=DTYPE) 
#     cdef double[::1] f1_ary = np.zeros(npixROI, dtype=DTYPE) 

#     cdef double[:,::1] nu_mat = np.zeros((k_max+1, npixROI), dtype=DTYPE)


#     cdef double ll = 0.

#     for p in range(npixROI):
#         f0_ary[p] = -(xbg_PSF_compressed[p] + x_m_sum[p])
#         f1_ary[p] = (xbg_PSF_compressed[p] + x_m_ary[1,p])
#         nu_mat[0,p] = exp(f0_ary[p])
#         nu_mat[1,p] = nu_mat[0,p] * f1_ary[p]
        
#         for k in range(2,data[p]+1):
#             for n in range(0, k - 1):
#                 nu_mat[k,p] += (k-n)/ float(k) * x_m_ary[k-n,p] * nu_mat[n,p]
#             nu_mat[k,p] += f1_ary[p] * nu_mat[k-1,p] / float(k)

#         #print nu_mat[data[p],p]
#     #    print 'At data[p] = ', data[p], 'nu_mat[data[p]] = ', nu_mat[data[p],p]

#         ll+=log( nu_mat[data[p],p])

#     if math.isnan(ll) ==True or math.isinf(ll) ==True:
#         ll = -10.1**10.

#     return ll



# import numpy as np
# cimport numpy as np
# import healpy as hp
# import mpmath as mp
# import math
# import matplotlib.pyplot as plt
# cimport cython
# #from libcpp cimport bool

# import pulsars.special as spc

# import logging
# logger = logging.getLogger(__name__)

# DTYPE = np.float
# ctypedef np.float_t DTYPE_t

# cdef extern from "math.h":
#     double log(double x) nogil
#     double exp(double x) nogil
#     double pow(double x, double y) nogil

# import time

# def sum_arrays(double[::1] array_1, double[::1] array_2):
#     return sum_arrays_int(array_1, array_2)

# def sub_arrays(double[::1] array_1, double[::1] array_2):
#     return sub_arrays_int(array_1, array_2)

# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# @cython.initializedcheck(False)
# cdef double[::1] sum_arrays_int(double[::1] array_1, double[::1] array_2):
#     cdef int n = len(array_1)
#     cdef double[::1] sum_array =  np.zeros(n,dtype=DTYPE)
#     cdef Py_ssize_t i
#     for i in range(n):
#         sum_array[i] = array_1[i] + array_2[i]

#     return sum_array

# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# @cython.initializedcheck(False)
# cdef double[::1] sub_arrays_int(double[::1] array_1, double[::1] array_2):
#     cdef int n = len(array_1)
#     cdef double[::1] sum_array =  np.zeros(n,dtype=DTYPE)
#     cdef Py_ssize_t i
#     for i in range(n):
#         sum_array[i] = array_1[i] - array_2[i]

#     return sum_array

# def log_nu_k_ary_PSF_exact_3_PS(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, double[::1] PS_dist_compressed2,double[::1] PS_dist_compressed3, int[::1] data, double Sc = 1000.0, double[::1] x_m_sum2_t = np.array([-10. ,-10.],dtype=DTYPE), double[:,::1] x_m_ary2_t = np.zeros((2,2), dtype=DTYPE),double[::1] x_m_sum3_t = np.array([-10. ,-10.],dtype=DTYPE), double[:,::1] x_m_ary3_t = np.zeros((2,2), dtype=DTYPE)):
#     return log_nu_k_ary_PSF_exact_3_PS_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, PS_dist_compressed2,PS_dist_compressed3, data, Sc, x_m_sum2_t , x_m_ary2_t,x_m_sum3_t , x_m_ary3_t)


# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# @cython.initializedcheck(False)
# cdef double log_nu_k_ary_PSF_exact_3_PS_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, double[::1] PS_dist_compressed2,double[::1] PS_dist_compressed3, int[::1] data, double Sc , double[::1] x_m_sum2_t , double[:,::1] x_m_ary2_t,double[::1] x_m_sum3_t , double[:,::1] x_m_ary3_t ):
#     #t_0 = time.time()
#     # print 'k_max is :', k_max
#     # print 'Number of pixels:', len(xbg_PSF_compressed)
#     # print 'theta: ', theta

#     #cdef np.float t_i = time.time()

#     #print 'here theta: ', np.asarray(theta)

#     cdef int k_max = np.max(data) + 1
#     cdef int len_theta = len(theta)

#     cdef double A = float(theta[0])
#     cdef double n1 = float(theta[1])
#     cdef double n2 = float(theta[2])
#     cdef double Sb = float(theta[3])

#     cdef double A2, n12, n22, Sb2, A3, n13, n23, Sb3
    
#     if len_theta > 4:
#         A2 = float(theta[4])
#         n12 = float(theta[5])
#         n22 = float(theta[6])
#         Sb2 = float(theta[7])

#     if len_theta > 8:
#         A3 = float(theta[8])
#         n13 = float(theta[9])
#         n23 = float(theta[10])
#         Sb3 = float(theta[11])

#     cdef int npixROI = len(xbg_PSF_compressed)

#     cdef double[:,::1] x_m_ary2 = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)
#     cdef double[::1] g1_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[::1] g2_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double x_m_ary_f2 = 0.0
#     cdef double x_m_sum_f2 = 0.0
#     cdef double f2 = 0.0
#     cdef double df_rho_div_f2 = 0.0
    
    
#     cdef double[:,::1] x_m_ary = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
#     cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double x_m_ary_f = 0.0
#     cdef double x_m_sum_f = 0.0

#     cdef double[:,::1] x_m_ary3 = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum3 = np.zeros(npixROI, dtype=DTYPE)
#     cdef double[::1] g1_ary_f3 = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[::1] g2_ary_f3 = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double x_m_ary_f3 = 0.0
#     cdef double x_m_sum_f3 = 0.0


#     cdef Py_ssize_t f_index, p, k, n


#     #calculations for PS
    
#     cdef int do_2 = 0
#     cdef int do_3 = 0

#     cdef double term1 = 0.0
#     cdef double term2 = 0.0
#     cdef double term3 = 0.0
#     cdef double second_3_a = 0.0
#     cdef double second_3_b = 0.0
#     cdef double second_3_c = 0.0
#     cdef double second_3_d = 0.0
#     cdef double second_2_a = 0.0
#     cdef double second_2_b = 0.0
#     cdef double second_2_c = 0.0
#     cdef double second_2_d = 0.0
#     cdef double second_1_a = 0.0
#     cdef double second_1_b = 0.0
#     cdef double second_1_c = 0.0
#     cdef double second_1_d = 0.0


#     if len_theta > 8:
#         x_m_ary2 = x_m_ary2_t
#         x_m_sum2 = x_m_sum2_t
#         do_2 = 1
#     else:
#         x_m_ary2 = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#         x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)

#     if len_theta > 4:
#         x_m_ary3 = x_m_ary3_t
#         x_m_sum3 = x_m_sum3_t
#         do_3 = 1
#     else:
#         x_m_ary3 = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#         x_m_sum3 = np.zeros(npixROI, dtype=DTYPE)
    
#     if do_2 and do_3:
#         for f_index in range(len(f_ary)):
#             f2 = f_ary[f_index]
#             df_rho_div_f2 = df_rho_div_f_ary[f_index]
#             # g1_ary_f = sub_arrays_int( spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) , spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2) )
#             g1_ary_f =  spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2) 
#             g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
#             term1 = (A * Sb * f2) \
#                                  * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
#                                     - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
#             second_1_a =  A  * pow(Sb * f2, n1)
#             second_1_b = A * pow(Sb * f2, n2)
#             for p in range(npixROI):
#                 x_m_sum_f = term1 * PS_dist_compressed[p]
#                 x_m_sum[p] += df_rho_div_f2*x_m_sum_f

#                 second_1_c = second_1_a * PS_dist_compressed[p]
#                 second_1_d = second_1_b * PS_dist_compressed[p]
#                 for k in range(data[p]+1):   #####take over here!!!
#                     x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
#                     x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f
#     elif do_2:
#         for f_index in range(len(f_ary)):
#             f2 = f_ary[f_index]
#             df_rho_div_f2 = df_rho_div_f_ary[f_index]
#             # g1_ary_f2 = sub_arrays_int(spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sb2 * f2),spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sc * f2))
#             g1_ary_f2 = spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sb2 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sc * f2)
#             g2_ary_f2 = spc.gammainc_lo_fct_ary_back(k_max, 1. - n22, Sb2 * f2)
#             g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
#             #g1_ary_f = sub_arrays_int(spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) , spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2) )
#             g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
#             term2 = (A2  * Sb2 * f2) \
#                                  * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
#                                     - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))

#             term1 = (A * Sb * f2) \
#                                  * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
#                                     - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
#             second_2_a =  A2  * pow(Sb2 * f2, n12)
#             second_2_b = A2 * pow(Sb2 * f2, n22)

#             second_1_a =  A  * pow(Sb * f2, n1)
#             second_1_b = A * pow(Sb * f2, n2)

#             for p in range(npixROI):
#                 x_m_sum_f2 = term2 * PS_dist_compressed2[p]
#                 x_m_sum2[p] += df_rho_div_f2*x_m_sum_f2
#                 x_m_sum_f = term1 * PS_dist_compressed[p]
#                 x_m_sum[p] += df_rho_div_f2*x_m_sum_f

#                 second_2_c = second_2_a * PS_dist_compressed2[p]
#                 second_2_d = second_2_b * PS_dist_compressed2[p]

#                 second_1_c = second_1_a * PS_dist_compressed[p]
#                 second_1_d = second_1_b * PS_dist_compressed[p]
#                 for k in range(data[p]+1):   #####take over here!!!
#                     x_m_ary_f2 = second_2_c  * g1_ary_f2[k] + second_2_d * g2_ary_f2[k]
#                     x_m_ary2[k,p] += df_rho_div_f2*x_m_ary_f2
#                     x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
#                     x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f
#     else:
#         for f_index in range(len(f_ary)):
#             f2 = f_ary[f_index]
#             df_rho_div_f2 = df_rho_div_f_ary[f_index]
#             g1_ary_f3 =  spc.gammainc_up_fct_ary_log(k_max, 1. - n13, Sb3 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n13, Sc * f2) 
#             g2_ary_f3 = spc.gammainc_lo_fct_ary_back(k_max, 1. - n23, Sb3 * f2)
#             g1_ary_f2 = spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sb2 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sc * f2) 
#             g2_ary_f2 = spc.gammainc_lo_fct_ary_back(k_max, 1. - n22, Sb2 * f2)
#             g1_ary_f =  spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2) 
#             g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
#             term3 = (A3  * Sb3 * f2) \
#                                  * (1./(n13-1.) + 1./(1.-n23) - pow(Sb3 / Sc, n13-1.)/(n13-1.) \
#                                     - (pow(Sb3 * f2, n13-1.) * g1_ary_f3[0] + pow(Sb3 * f2, n23-1.) * g2_ary_f3[0]))
#             term2 = (A2  * Sb2 * f2) \
#                                  * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
#                                     - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))

#             term1 = (A * Sb * f2) \
#                                  * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
#                                     - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
#             second_3_a =  A3  * pow(Sb3 * f2, n13)
#             second_3_b = A3 * pow(Sb3 * f2, n23)

#             second_2_a =  A2  * pow(Sb2 * f2, n12)
#             second_2_b = A2 * pow(Sb2 * f2, n22)

#             second_1_a =  A  * pow(Sb * f2, n1)
#             second_1_b = A * pow(Sb * f2, n2)

#             for p in range(npixROI):
#                 x_m_sum_f3 = term3 * PS_dist_compressed3[p]
#                 x_m_sum3[p] += df_rho_div_f2*x_m_sum_f3
#                 x_m_sum_f2 = term2 * PS_dist_compressed2[p]
#                 x_m_sum2[p] += df_rho_div_f2*x_m_sum_f2
#                 x_m_sum_f = term1 * PS_dist_compressed[p]
#                 x_m_sum[p] += df_rho_div_f2*x_m_sum_f

#                 second_3_c = second_3_a * PS_dist_compressed3[p]
#                 second_3_d = second_3_b * PS_dist_compressed3[p]

#                 second_2_c = second_2_a * PS_dist_compressed2[p]
#                 second_2_d = second_2_b * PS_dist_compressed2[p]

#                 second_1_c = second_1_a * PS_dist_compressed[p]
#                 second_1_d = second_1_b * PS_dist_compressed[p]
#                 for k in range(data[p]+1):   #####take over here!!!
#                     x_m_ary_f3 = second_3_c  * g1_ary_f3[k] + second_3_d * g2_ary_f3[k]
#                     x_m_ary3[k,p] += df_rho_div_f2*x_m_ary_f3
#                     x_m_ary_f2 = second_2_c  * g1_ary_f2[k] + second_2_d * g2_ary_f2[k]
#                     x_m_ary2[k,p] += df_rho_div_f2*x_m_ary_f2
#                     x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
#                     x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f

    
#     cdef double[::1] x_m_ary_total = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double x_m_sum_total = 0.0


#     cdef double[::1] nu_ary = np.zeros(k_max + 1, dtype=DTYPE)
    
#     cdef double f0_ary = 0.0 # -(xbg_PSF_compressed + x_m_sum_total)
#     cdef double f1_ary = 0.0 # (xbg_PSF_compressed + x_m_ary_total[1])

#     cdef double[::1] nu_mat = np.zeros(k_max+1, dtype=DTYPE)

#     #t0 = time.time()
#     cdef double ll = 0.
#     #cdef double temp = 0.
#     for p in range(npixROI):
#         x_m_sum_total = x_m_sum3[p]+x_m_sum2[p]+ x_m_sum[p]
#         x_m_ary_total[0] = x_m_ary3[0,p]+x_m_ary2[0,p] + x_m_ary[0,p]
#         x_m_ary_total[1] = x_m_ary3[1,p]+x_m_ary2[1,p] + x_m_ary[1,p]
#         f0_ary = -(xbg_PSF_compressed[p] + x_m_sum_total)
#         f1_ary = (xbg_PSF_compressed[p] + x_m_ary_total[1])
#         nu_mat[0] = exp(f0_ary)
#         nu_mat[1] = nu_mat[0] * f1_ary
#         #print nu_mat[1,p]
#         for k in range(2,data[p]+1):
#             x_m_ary_total[k] = x_m_ary3[k,p]+x_m_ary2[k,p] + x_m_ary[k,p]
#             nu_mat[k] = 0.0
#             for n in range(0, k - 1):
#                 nu_mat[k] += (k-n)/ float(k) * x_m_ary_total[k-n] * nu_mat[n]
#             nu_mat[k] += f1_ary * nu_mat[k-1] / float(k)

#         ll+=log( nu_mat[data[p]])


#     if math.isnan(ll) ==True or math.isinf(ll) ==True:
#         ll = -10.1**10.

#     return ll

# def log_nu_k_ary_PSF_exact_2_PS(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, double[::1] PS_dist_compressed2, int[::1] data, double Sc = 1000.0, double[::1] x_m_sum2_t = np.array([-10. ,-10.],dtype=DTYPE), double[:,::1] x_m_ary2_t = np.zeros((2,2), dtype=DTYPE)):
#     return log_nu_k_ary_PSF_exact_2_PS_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, PS_dist_compressed2, data, Sc, x_m_sum2_t , x_m_ary2_t)


# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# @cython.initializedcheck(False)
# cdef double log_nu_k_ary_PSF_exact_2_PS_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, double[::1] PS_dist_compressed2, int[::1] data, double Sc , double[::1] x_m_sum2_t , double[:,::1] x_m_ary2_t ):
#     #t_0 = time.time()
#     # print 'k_max is :', k_max
#     # print 'Number of pixels:', len(xbg_PSF_compressed)
#     # print 'theta: ', theta

#     #cdef np.float t_i = time.time()

#     cdef int k_max = np.max(data) + 1

#     cdef int len_theta = len(theta)

#     cdef double A = float(theta[0])
#     cdef double n1 = float(theta[1])
#     cdef double n2 = float(theta[2])
#     cdef double Sb = float(theta[3])
    
#     cdef double A2, n12, n22, Sb2
    
#     if len_theta > 4:
#         A2 = float(theta[4])
#         n12 = float(theta[5])
#         n22 = float(theta[6])
#         Sb2 = float(theta[7])

#     cdef int npixROI = len(xbg_PSF_compressed)

#     cdef double[:,::1] x_m_ary2 = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)
#     cdef double[::1] g1_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[::1] g2_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double x_m_ary_f2 = 0.0
#     cdef double x_m_sum_f2 = 0.0
#     cdef double f2 = 0.0
#     cdef double df_rho_div_f2 = 0.0
    
    
#     cdef double[:,::1] x_m_ary = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
#     cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double x_m_ary_f = 0.0
#     cdef double x_m_sum_f = 0.0
#     # cdef double f = 0.0
#     # cdef double df_rho_div_f = 0.0
    

#     #k_max = int(np.max(data))

#     cdef Py_ssize_t f_index, p, k, n

#     #cdef np.float t_0 = time.time()

    
    
#     #calculations for PS
    
#     cdef int do_half = 0

#     cdef double term1 = 0.0
#     cdef double term2 = 0.0
#     cdef double second_2_a = 0.0
#     cdef double second_2_b = 0.0
#     cdef double second_2_c = 0.0
#     cdef double second_2_d = 0.0
#     cdef double second_1_a = 0.0
#     cdef double second_1_b = 0.0
#     cdef double second_1_c = 0.0
#     cdef double second_1_d = 0.0


#     if x_m_sum2_t[0] != -10.0:
#         x_m_ary2 = x_m_ary2_t
#         x_m_sum2 = x_m_sum2_t
#         do_half = 1
#     else:
#         x_m_ary2 = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#         x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)
    
#     if do_half:
#         for f_index in range(len(f_ary)):
#             f2 = f_ary[f_index]
#             df_rho_div_f2 = df_rho_div_f_ary[f_index]
#             g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
#             g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
#             term1 = (A * Sb * f2) \
#                                  * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
#                                     - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
#             second_1_a =  A  * pow(Sb * f2, n1)
#             second_1_b = A * pow(Sb * f2, n2)
#             for p in range(npixROI):
#                 x_m_sum_f = term1 * PS_dist_compressed[p]
#                 x_m_sum[p] += df_rho_div_f2*x_m_sum_f

#                 second_1_c = second_1_a * PS_dist_compressed[p]
#                 second_1_d = second_1_b * PS_dist_compressed[p]
#                 for k in range(data[p]+1):   #####take over here!!!
#                     x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
#                     x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f
#     else:
#         #t0 = time.time()
#         for f_index in range(len(f_ary)):
#             f2 = f_ary[f_index]
#             df_rho_div_f2 = df_rho_div_f_ary[f_index]
#             #t0 = time.time()
#             g1_ary_f2 = spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sb2 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sc * f2) 
#             g2_ary_f2 = spc.gammainc_lo_fct_ary_back(k_max, 1. - n22, Sb2 * f2)
#             g1_ary_f =  spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2) 
#             g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
#            # ta += time.time() - t0
#             #print 'g1, g2 time: ', time.time()-t0
#             term2 = (A2  * Sb2 * f2) \
#                                  * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
#                                     - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))

#             term1 = (A * Sb * f2) \
#                                  * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
#                                     - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
#             second_2_a =  A2  * pow(Sb2 * f2, n12)
#             second_2_b = A2 * pow(Sb2 * f2, n22)

#             second_1_a =  A  * pow(Sb * f2, n1)
#             second_1_b = A * pow(Sb * f2, n2)

#             for p in range(npixROI):
#                 #t0 = time.time()
#                 # x_m_sum_f2[p] = (A2 * PS_dist_compressed2[p] * Sb2 * f2) \
#                 #                  * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
#                 #                     - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))
#                 x_m_sum_f2 = term2 * PS_dist_compressed2[p]
#                 x_m_sum2[p] += df_rho_div_f2*x_m_sum_f2
#                 x_m_sum_f = term1 * PS_dist_compressed[p]
#                 x_m_sum[p] += df_rho_div_f2*x_m_sum_f
#                 #tb += time.time() - t0

#                 second_2_c = second_2_a * PS_dist_compressed2[p]
#                 second_2_d = second_2_b * PS_dist_compressed2[p]

#                 second_1_c = second_1_a * PS_dist_compressed[p]
#                 second_1_d = second_1_b * PS_dist_compressed[p]
#                 for k in range(data[p]+1):   #####take over here!!!
#                     #t0 = time.time()
#                     x_m_ary_f2 = second_2_c  * g1_ary_f2[k] + second_2_d * g2_ary_f2[k]
#                     #print 'At k,p = ', k, p, 'xm_ary_f2 = ', x_m_ary_f2[k,p]
#                     x_m_ary2[k,p] += df_rho_div_f2*x_m_ary_f2
#                     x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
#                     #print 'At k,p = ', k, p, 'xm_ary_f = ', x_m_ary_f[k,p]
#                     x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f
    
#     cdef double[::1] x_m_ary_total = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double x_m_sum_total = 0.0


#     cdef double[::1] nu_ary = np.zeros(k_max + 1, dtype=DTYPE)
    
#     cdef double f0_ary = 0.0 # -(xbg_PSF_compressed + x_m_sum_total)
#     cdef double f1_ary = 0.0 # (xbg_PSF_compressed + x_m_ary_total[1])

#     cdef double[::1] nu_mat = np.zeros(k_max+1, dtype=DTYPE)

#     #t0 = time.time()
#     cdef double ll = 0.
#     #cdef double temp = 0.
#     for p in range(npixROI):
#         x_m_sum_total = x_m_sum2[p]+ x_m_sum[p]
#         x_m_ary_total[0] = x_m_ary2[0,p] + x_m_ary[0,p]
#         x_m_ary_total[1] = x_m_ary2[1,p] + x_m_ary[1,p]
#         f0_ary = -(xbg_PSF_compressed[p] + x_m_sum_total)
#         f1_ary = (xbg_PSF_compressed[p] + x_m_ary_total[1])
#         nu_mat[0] = exp(f0_ary)
#         nu_mat[1] = nu_mat[0] * f1_ary
#         #print nu_mat[1,p]
#         for k in range(2,data[p]+1):
#             x_m_ary_total[k] = x_m_ary2[k,p] + x_m_ary[k,p]
#             nu_mat[k] = 0.0
#             for n in range(0, k - 1):
#                 nu_mat[k] += (k-n)/ float(k) * x_m_ary_total[k-n] * nu_mat[n]
#             nu_mat[k] += f1_ary * nu_mat[k-1] / float(k)

#         ll+=log( nu_mat[data[p]])


#     if math.isnan(ll) ==True or math.isinf(ll) ==True:
#         ll = -10.1**10.

#     return ll

# # @cython.boundscheck(False)
# # @cython.wraparound(False)
# # @cython.cdivision(True)
# # @cython.initializedcheck(False)
# # cdef double log_nu_k_ary_PSF_exact_2_PS_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, double[::1] PS_dist_compressed2, int[::1] data, double Sc , double[::1] x_m_sum2_t , double[:,::1] x_m_ary2_t ):
# #     #t_0 = time.time()
# #     # print 'k_max is :', k_max
# #     # print 'Number of pixels:', len(xbg_PSF_compressed)
# #     # print 'theta: ', theta

# #     #cdef np.float t_i = time.time()

# #     cdef int k_max = np.max(data) + 1

# #     cdef double A = np.float(theta[0])
# #     cdef double n1 = np.float(theta[1])
# #     cdef double n2 = np.float(theta[2])
# #     cdef double Sb = np.float(theta[3])
    
# #     cdef double A2 = np.float(theta[4])
# #     cdef double n12 = np.float(theta[5])
# #     cdef double n22 = np.float(theta[6])
# #     cdef double Sb2 = np.float(theta[7])

# #     cdef int npixROI = len(xbg_PSF_compressed)

# #     cdef double[:,::1] x_m_ary2 = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
# #     cdef double[::1] x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)
# #     cdef double[::1] g1_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
# #     cdef double[::1] g2_ary_f2 = np.zeros(k_max + 1, dtype=DTYPE)
# #     cdef double[:,::1] x_m_ary_f2 = np.zeros((k_max + 1, npixROI), dtype=DTYPE)
# #     cdef double[::1] x_m_sum_f2 = np.zeros(npixROI, dtype=DTYPE)
# #     cdef double f2 = 0.0
# #     cdef double df_rho_div_f2 = 0.0
    
    
# #     cdef double[:,::1] x_m_ary = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
# #     cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
# #     cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
# #     cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
# #     cdef double[:,::1] x_m_ary_f = np.zeros((k_max + 1, npixROI), dtype=DTYPE)
# #     cdef double[::1] x_m_sum_f = np.zeros(npixROI, dtype=DTYPE)
# #     # cdef double f = 0.0
# #     # cdef double df_rho_div_f = 0.0
    

# #     #k_max = int(np.max(data))

# #     cdef Py_ssize_t f_index, p, k, n

# #     #cdef np.float t_0 = time.time()

    
    
# #     #calculations for PS
    
# #     cdef int do_half = 0

# #     cdef double term1 = 0.0
# #     cdef double term2 = 0.0
# #     cdef double second_2_a = 0.0
# #     cdef double second_2_b = 0.0
# #     cdef double second_2_c = 0.0
# #     cdef double second_2_d = 0.0
# #     cdef double second_1_a = 0.0
# #     cdef double second_1_b = 0.0
# #     cdef double second_1_c = 0.0
# #     cdef double second_1_d = 0.0


# #     if x_m_sum2_t[0] != -10.0:
# #         x_m_ary2 = x_m_ary2_t
# #         x_m_sum2 = x_m_sum2_t
# #         do_half = 1
# #     else:
# #         x_m_ary2 = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
# #         x_m_sum2 = np.zeros(npixROI, dtype=DTYPE)
    
# #     if do_half:
# #         for f_index in range(len(f_ary)):
# #             f2 = f_ary[f_index]
# #             df_rho_div_f2 = df_rho_div_f_ary[f_index]
# #             g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
# #             g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
# #             for p in range(npixROI):
# #                 x_m_sum_f[p] = (A * PS_dist_compressed[p] * Sb * f2) \
# #                                  * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
# #                                     - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
# #                 x_m_sum[p] += df_rho_div_f2*x_m_sum_f[p]
# #                 for k in range(data[p]+1):   #####take over here!!!
# #                     x_m_ary_f[k,p] = A  * (pow(Sb * f2, n1)  * g1_ary_f[k] * PS_dist_compressed[p] + pow(Sb * f2, n2) * g2_ary_f[k]*PS_dist_compressed[p])
# #                     x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f[k,p]
# #     else:
# #         #t0 = time.time()
# #         for f_index in range(len(f_ary)):
# #             f2 = f_ary[f_index]
# #             df_rho_div_f2 = df_rho_div_f_ary[f_index]
# #             #t0 = time.time()
# #             g1_ary_f2 = spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sb2 * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n12, Sc * f2)
# #             g2_ary_f2 = spc.gammainc_lo_fct_ary_back(k_max, 1. - n22, Sb2 * f2)
# #             g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
# #             g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
# #            # ta += time.time() - t0
# #             #print 'g1, g2 time: ', time.time()-t0
# #             term2 = (A2  * Sb2 * f2) \
# #                                  * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
# #                                     - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))

# #             term1 = (A * Sb * f2) \
# #                                  * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
# #                                     - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
# #             second_2_a =  A2  * pow(Sb2 * f2, n12)
# #             second_2_b = A2 * pow(Sb2 * f2, n22)

# #             second_1_a =  A  * pow(Sb * f2, n1)
# #             second_1_b = A * pow(Sb * f2, n2)

# #             for p in range(npixROI):
# #                 #t0 = time.time()
# #                 # x_m_sum_f2[p] = (A2 * PS_dist_compressed2[p] * Sb2 * f2) \
# #                 #                  * (1./(n12-1.) + 1./(1.-n22) - pow(Sb2 / Sc, n12-1.)/(n12-1.) \
# #                 #                     - (pow(Sb2 * f2, n12-1.) * g1_ary_f2[0] + pow(Sb2 * f2, n22-1.) * g2_ary_f2[0]))
# #                 x_m_sum_f2[p] = term2 * PS_dist_compressed2[p]
# #                 x_m_sum2[p] += df_rho_div_f2*x_m_sum_f2[p]
# #                 x_m_sum_f[p] = term1 * PS_dist_compressed[p]
# #                 x_m_sum[p] += df_rho_div_f2*x_m_sum_f[p]
# #                 #tb += time.time() - t0

# #                 second_2_c = second_2_a * PS_dist_compressed2[p]
# #                 second_2_d = second_2_b * PS_dist_compressed2[p]

# #                 second_1_c = second_1_a * PS_dist_compressed[p]
# #                 second_1_d = second_1_b * PS_dist_compressed[p]
# #                 for k in range(data[p]+1):   #####take over here!!!
# #                     #t0 = time.time()
# #                     x_m_ary_f2[k,p] = second_2_c  * g1_ary_f2[k] + second_2_d * g2_ary_f2[k]
# #                     #print 'At k,p = ', k, p, 'xm_ary_f2 = ', x_m_ary_f2[k,p]
# #                     x_m_ary2[k,p] += df_rho_div_f2*x_m_ary_f2[k,p]
# #                     x_m_ary_f[k,p] = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
# #                     #print 'At k,p = ', k, p, 'xm_ary_f = ', x_m_ary_f[k,p]
# #                     x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f[k,p]
    
# #     cdef double[:,::1] x_m_ary_total = np.zeros((k_max + 1, npixROI), dtype=DTYPE)
# #     cdef double[::1] x_m_sum_total = np.zeros(npixROI, dtype=DTYPE)


# #     cdef double[::1] nu_ary = np.zeros(k_max + 1, dtype=DTYPE)
    
# #     cdef double[::1] f0_ary = np.zeros(npixROI, dtype=DTYPE) # -(xbg_PSF_compressed + x_m_sum_total)
# #     cdef double[::1] f1_ary = np.zeros(npixROI, dtype=DTYPE) # (xbg_PSF_compressed + x_m_ary_total[1])

# #     cdef double[:,::1] nu_mat = np.zeros((k_max+1, npixROI), dtype=DTYPE)

# #     #t0 = time.time()
# #     cdef double ll = 0.
# #     #cdef double temp = 0.
# #     for p in range(npixROI):
# #         x_m_sum_total[p] = x_m_sum2[p]+ x_m_sum[p]
# #         x_m_ary_total[0,p] = x_m_ary2[0,p] + x_m_ary[0,p]
# #         x_m_ary_total[1,p] = x_m_ary2[1,p] + x_m_ary[1,p]
# #         f0_ary[p] = -(xbg_PSF_compressed[p] + x_m_sum_total[p])
# #         f1_ary[p] = (xbg_PSF_compressed[p] + x_m_ary_total[1,p])
# #         nu_mat[0,p] = exp(f0_ary[p])
# #         nu_mat[1,p] = nu_mat[0,p] * f1_ary[p]
# #         #print nu_mat[1,p]
# #         for k in range(2,data[p]+1):
# #             x_m_ary_total[k,p] = x_m_ary2[k,p] + x_m_ary[k,p]
# #             for n in range(0, k - 1):
# #                 nu_mat[k,p] += (k-n)/ float(k) * x_m_ary_total[k-n,p] * nu_mat[n,p]
# #             nu_mat[k,p] += f1_ary[p] * nu_mat[k-1,p] / float(k)

# #         ll+=log( nu_mat[data[p],p])


# #     if math.isnan(ll) ==True or math.isinf(ll) ==True:
# #         ll = -10.1**10.


# #     return ll


# # def log_nu_k_ary_PSF_exact_1_PS_edep(double[:,::1] xbg_PSF_compressed, double[:,::1] theta, double[:,::1] f_ary, double[:,::1] df_rho_div_f_ary, double[:,::1] PS_dist_compressed, int[:,::1] data, double Sc = 1000.0):
# #     return log_nu_k_ary_PSF_exact_1_PS_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data, Sc)


# def log_nu_k_ary_PSF_exact_1_PS(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data, double Sc = 1000.0):
#     return log_nu_k_ary_PSF_exact_1_PS_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data, Sc)


# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)
# @cython.initializedcheck(False)
# cdef double log_nu_k_ary_PSF_exact_1_PS_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data, double Sc ):

#     cdef int k_max = np.max(data) + 1

#     cdef double A = np.float(theta[0])
#     cdef double n1 = np.float(theta[1])
#     cdef double n2 = np.float(theta[2])
#     cdef double Sb = np.float(theta[3])

#     cdef int npixROI = len(xbg_PSF_compressed)

#     cdef double f2 = 0.0
#     cdef double df_rho_div_f2 = 0.0


#     cdef double[:,::1] x_m_ary = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
#     cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
#     cdef double x_m_ary_f
#     cdef double x_m_sum_f

#     cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)

#     cdef Py_ssize_t f_index, p, k, n


#     #calculations for PS

#     cdef int do_half = 0

#     cdef double term1 = 0.0
#     cdef double term2 = 0.0
#     cdef double second_2_a = 0.0
#     cdef double second_2_b = 0.0
#     cdef double second_2_c = 0.0
#     cdef double second_2_d = 0.0
#     cdef double second_1_a = 0.0
#     cdef double second_1_b = 0.0
#     cdef double second_1_c = 0.0
#     cdef double second_1_d = 0.0



#     for f_index in range(len(f_ary)):
#         f2 = f_ary[f_index]
#         df_rho_div_f2 = df_rho_div_f_ary[f_index]
#         g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
#         g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
#         term1 = (A * Sb * f2) \
#                              * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
#                                 - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
#         second_1_a =  A  * pow(Sb * f2, n1)
#         second_1_b = A * pow(Sb * f2, n2)

#         for p in range(npixROI):
#             x_m_sum_f = term1 * PS_dist_compressed[p]
#             x_m_sum[p] += df_rho_div_f2*x_m_sum_f

#             second_1_c = second_1_a * PS_dist_compressed[p]
#             second_1_d = second_1_b * PS_dist_compressed[p]
#             for k in range(data[p]+1):
#                 x_m_ary_f = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k] 
#                 x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f


#     cdef double[::1] nu_ary = np.zeros(k_max + 1, dtype=DTYPE)

#     cdef double f0_ary
#     cdef double f1_ary

#     cdef double[:] nu_mat = np.zeros((k_max+1), dtype=DTYPE)


#     cdef double ll = 0.

#     for p in range(npixROI):
#         f0_ary = -(xbg_PSF_compressed[p] + x_m_sum[p])
#         f1_ary = (xbg_PSF_compressed[p] + x_m_ary[1,p])
#         nu_mat[0] = exp(f0_ary)
#         nu_mat[1] = nu_mat[0] * f1_ary

#         for k in range(2,data[p]+1):
#             nu_mat[k] = 0.0
#             for n in range(0, k - 1):
#                 nu_mat[k] += (k-n)/ float(k) * x_m_ary[k-n,p] * nu_mat[n]
#             nu_mat[k] += f1_ary * nu_mat[k-1] / float(k)

#      #   print 'At data[p] = ', data[p], 'nu_mat[data[p]] = ', nu_mat[data[p]]

#         ll+=log( nu_mat[data[p]])

#     if math.isnan(ll) or math.isinf(ll):
#         ll = -10.1**10.

#     return ll

# #=====================================
# #Just return x's 
# def return_xs(list theta, np.ndarray[DTYPE_t] f_ary, np.ndarray[DTYPE_t] df_rho_div_f_ary, np.ndarray PS_dist_compressed, np.ndarray data, np.float Sc = 1000.0):
#     cdef int k_max = int(np.max(data) + 1)

#     cdef np.float A = np.float(theta[0])
#     cdef np.float n1 = np.float(theta[1])
#     cdef np.float n2 = np.float(theta[2])
#     cdef np.float Sb = np.float(theta[3])

#     cdef np.int npixROI = len(PS_dist_compressed)
    
#     cdef np.ndarray x_m_ary = np.zeros((k_max + 1, npixROI), dtype=DTYPE)
#     cdef np.ndarray x_m_sum = np.zeros(npixROI, dtype=DTYPE)
    
#     cdef np.ndarray g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef np.ndarray g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
#     cdef np.ndarray x_m_ary_f = np.zeros((k_max + 1, npixROI), dtype=DTYPE)
#     cdef np.ndarray x_m_sum_f = np.zeros(npixROI, dtype=DTYPE)
#     cdef np.float f = 0.0
#     cdef np.float df_rho_div_f = 0.0
    
#     cdef Py_ssize_t f_index2
#     for f_index2 in range(len(f_ary)):
#         f = f_ary[f_index2]
#         df_rho_div_f = df_rho_div_f_ary[f_index2]
#         #g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f)
#         g1_ary_f =  spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f) 
#         g2_ary_f =  spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f) 
#         x_m_ary_f = A  * (np.power(Sb * f, n1) * np.outer(g1_ary_f, PS_dist_compressed) + np.power(Sb * f, n2) * np.outer(g2_ary_f, PS_dist_compressed))
#         x_m_sum_f = (A * PS_dist_compressed * Sb * f) \
#             * (1./(n1-1.) + 1./(1.-n2) - np.power(Sb / Sc, n1-1.)/(n1-1.)\
#                - (np.power(Sb * f, n1-1.) * g1_ary_f[0] + np.power(Sb * f, n2-1.) * g2_ary_f[0]))


#         x_m_ary += df_rho_div_f*x_m_ary_f
#         x_m_sum += df_rho_div_f*x_m_sum_f

#     return x_m_ary, x_m_sum

# # def log_nu_k_ary_PSF_exact_1_PS(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data, double Sc = 1000.0):
# #     return log_nu_k_ary_PSF_exact_1_PS_int(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data, Sc)

# # @cython.boundscheck(False)
# # @cython.wraparound(False)
# # @cython.cdivision(True)
# # @cython.initializedcheck(False)
# # cdef double log_nu_k_ary_PSF_exact_1_PS_int(double[::1] xbg_PSF_compressed, double[::1] theta, double[::1] f_ary, double[::1] df_rho_div_f_ary, double[::1] PS_dist_compressed, int[::1] data, double Sc ):

# #     cdef int k_max = np.max(data) + 1

# #     cdef double A = np.float(theta[0])
# #     cdef double n1 = np.float(theta[1])
# #     cdef double n2 = np.float(theta[2])
# #     cdef double Sb = np.float(theta[3])

# #     cdef int npixROI = len(xbg_PSF_compressed)

# #     cdef double f2 = 0.0
# #     cdef double df_rho_div_f2 = 0.0
    
    
# #     cdef double[:,::1] x_m_ary = np.zeros((k_max + 1,npixROI), dtype=DTYPE)
# #     cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
# #     cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
# #     cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
# #     cdef double[:,::1] x_m_ary_f = np.zeros((k_max + 1, npixROI), dtype=DTYPE)
# #     cdef double[::1] x_m_sum_f = np.zeros(npixROI, dtype=DTYPE)

# #     cdef Py_ssize_t f_index, p, k, n
    
    
# #     #calculations for PS
    
# #     cdef int do_half = 0

# #     cdef double term1 = 0.0
# #     cdef double term2 = 0.0
# #     cdef double second_2_a = 0.0
# #     cdef double second_2_b = 0.0
# #     cdef double second_2_c = 0.0
# #     cdef double second_2_d = 0.0
# #     cdef double second_1_a = 0.0
# #     cdef double second_1_b = 0.0
# #     cdef double second_1_c = 0.0
# #     cdef double second_1_d = 0.0

    

# #     for f_index in range(len(f_ary)):
# #         f2 = f_ary[f_index]
# #         df_rho_div_f2 = df_rho_div_f_ary[f_index]
# #         g1_ary_f = spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sb * f2) - spc.gammainc_up_fct_ary_log(k_max, 1. - n1, Sc * f2)
# #         g2_ary_f = spc.gammainc_lo_fct_ary_back(k_max, 1. - n2, Sb * f2)
# #         term1 = (A * Sb * f2) \
# #                              * (1./(n1-1.) + 1./(1.-n2) - pow(Sb / Sc, n1-1.)/(n1-1.) \
# #                                 - (pow(Sb * f2, n1-1.) * g1_ary_f[0] + pow(Sb * f2, n2-1.) * g2_ary_f[0]))
# #         second_1_a =  A  * pow(Sb * f2, n1)
# #         second_1_b = A * pow(Sb * f2, n2)

# #         for p in range(npixROI):
# #             x_m_sum_f[p] = term1 * PS_dist_compressed[p]
# #             x_m_sum[p] += df_rho_div_f2*x_m_sum_f[p]

# #             second_1_c = second_1_a * PS_dist_compressed[p]
# #             second_1_d = second_1_b * PS_dist_compressed[p]
# #             for k in range(data[p]+1):            
# #                 x_m_ary_f[k,p] = second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k]
# #                 x_m_ary[k,p] += df_rho_div_f2*x_m_ary_f[k,p]

                
# #     cdef double[::1] nu_ary = np.zeros(k_max + 1, dtype=DTYPE)
    
# #     cdef double[::1] f0_ary = np.zeros(npixROI, dtype=DTYPE) 
# #     cdef double[::1] f1_ary = np.zeros(npixROI, dtype=DTYPE) 

# #     cdef double[:,::1] nu_mat = np.zeros((k_max+1, npixROI), dtype=DTYPE)


# #     cdef double ll = 0.

# #     for p in range(npixROI):
# #         f0_ary[p] = -(xbg_PSF_compressed[p] + x_m_sum[p])
# #         f1_ary[p] = (xbg_PSF_compressed[p] + x_m_ary[1,p])
# #         nu_mat[0,p] = exp(f0_ary[p])
# #         nu_mat[1,p] = nu_mat[0,p] * f1_ary[p]
        
# #         for k in range(2,data[p]+1):
# #             for n in range(0, k - 1):
# #                 nu_mat[k,p] += (k-n)/ float(k) * x_m_ary[k-n,p] * nu_mat[n,p]
# #             nu_mat[k,p] += f1_ary[p] * nu_mat[k-1,p] / float(k)

# #         #print nu_mat[data[p],p]
# #     #    print 'At data[p] = ', data[p], 'nu_mat[data[p]] = ', nu_mat[data[p],p]

# #         ll+=log( nu_mat[data[p],p])

# #     if math.isnan(ll) ==True or math.isinf(ll) ==True:
# #         ll = -10.1**10.

# #     return ll