#------------------------------------------------
# NPTF Log Likelihood
#------------------------------------------------
##
## Calculation of Log Likelihood for a template
## fit involving non-poissonian templates
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
    return log_like_internal(xbg_PSF_compressed,theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data)


## Define cython function
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_like_internal(double[::1] xbg_PSF_compressed, 
                              double[::1] theta, 
                              double[::1] f_ary, 
                              double[::1] df_rho_div_f_ary, 
                              double[::1] PS_dist_compressed, 
                              int[::1] data):

    # Unpack theta into individual NPT parameters
    cdef double A = np.float(theta[0])
    cdef double n1 = np.float(theta[1])
    cdef double n2 = np.float(theta[2])
    cdef double n3 = np.float(theta[3])
    cdef double Sb1 = np.float(theta[4])
    cdef double Sb2 = np.float(theta[5])

    # Extract variables used to set array lengths
    cdef int k_max = np.max(data) + 1
    cdef int npixROI = len(xbg_PSF_compressed)

    # Setup variables used for indexing loops
    cdef Py_ssize_t f_index, p, k, n


    ## Calculate the required x_m values
    
    cdef double f2 = 0., df_rho_div_f2 = 0.
    cdef double term1 = 0., term2 = 0.
    cdef double second_1_alpha = 0., second_1_beta = 0.
    cdef double second_1_a = 0., second_1_b = 0.
    cdef double second_1_c = 0., second_1_d = 0.
    
    cdef double[:,::1] x_m_ary = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)
    cdef double x_m_ary_f, x_m_sum_f

    cdef double[::1] g0_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g1_ary_f = np.zeros(k_max + 1, dtype=DTYPE)
    cdef double[::1] g2_ary_f = np.zeros(k_max + 1, dtype=DTYPE)

    # Loop through f to perform the integral over df
    for f_index in range(len(f_ary)):
        f2 = f_ary[f_index]
        df_rho_div_f2 = df_rho_div_f_ary[f_index]
        g0_ary_f = igf.incgamma_up_fct_ary(k_max, 1. - n1, Sb2 * f2) 
        g1_ary_f = igf.incgamma_up_fct_ary(k_max, 1. - n2, Sb1 * f2) - igf.incgamma_up_fct_ary(k_max, 1. - n2, Sb2 * f2)
        g2_ary_f = igf.incgamma_lo_fct_ary(k_max, 1. - n3, Sb1 * f2)
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
            # Start at 1 as we don't need x_0
            for k in range(1,data[p]+1):
                x_m_ary_f = second_1_beta * g0_ary_f[k] + second_1_c  * g1_ary_f[k] + second_1_d * g2_ary_f[k] 
                x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f


    ## Determine log likelihood for model to see data
    # This amounts to calculating p_k in every pixel
    # and adding their log likelihoods

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
        # if so then penalise
        if pk[data[p]] > 0:
            ll += log(pk[data[p]])
        else:
            ll += -10.1**10.

    return ll
