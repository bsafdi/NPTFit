###############################################################################
# npll.pyx
###############################################################################
#
# Calculation of non-poissonian contribution to Log Likelihood, which is then
# combined with the poissonian contribution using the method of generating
# functions.
#
# Calculation broken into two parts:
#  1. Determination of x_m and x_m_sum, which depends on the number of 
#     non-poissonian templates & the number of breaks each one has; and
#  2. Determination of LL in terms of these using a recurrence relation
#
###############################################################################

import numpy as np
cimport numpy as np
cimport cython
import x_m

# Type used for all non-integer functions
DTYPE = np.float

# Setup cython functions
cdef extern from "math.h":
    double log(double x) nogil
    double exp(double x) nogil


def log_like(pt_sum_compressed, theta, f_ary, df_rho_div_f_ary, npt_compressed,
             data):
    """ Python wrapper for the log likelihood

    Organises the calculation of x_m values for multiple non-poissonian 
    templates

    Args:
        pt_sum_compressed: pixel-wise sum of Poissonian model templates
        theta: Array of non-poissonian parameters. Length of array is the
               number of NP templates, and each element in the array has
               the form: [A, n[1], .., n[j+1], Sb[1], .., Sb[j]]
        f_ary: Photon leakage probabilities characterizing PSF, sum(f_ary) = 1.0
        df_rho_div_f_ary: df*rho(f)/f for integrating over f as a sum
        npt_compressed: Array of non-poissonian templates. Length of array is
                        the number of PS templates, and each element in the
                        array is the pixel-wise normalization of the PS template
        data: The pixel-wise data

    Returns:
        double log likelihood

    """

    cdef int k_max = np.max(data) + 1
    cdef int npixROI = len(pt_sum_compressed)
    cdef double[:,::1] x_m_ary = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum = np.zeros(npixROI, dtype=DTYPE)

    cdef Py_ssize_t i
    cdef double[:,::1] x_m_ary_out = np.zeros((npixROI,k_max + 1), dtype=DTYPE)
    cdef double[::1] x_m_sum_out = np.zeros(npixROI, dtype=DTYPE)

    for i in range(len(theta)):
        # Check theta has the correct length
        assert(len(theta[i]) % 2 == 0), "theta has an invalid length!"

        x_m_ary_out, x_m_sum_out = x_m.return_xs(np.array(theta[i]), f_ary, df_rho_div_f_ary, 
                                       npt_compressed[i], data)
        x_m_ary += np.asarray(x_m_ary_out)
        x_m_sum += np.asarray(x_m_sum_out)

    return log_like_internal(pt_sum_compressed, data, x_m_ary, x_m_sum, k_max,
                             npixROI)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef double log_like_internal(double[::1] pt_sum_compressed, int[::1] data, 
                              double[:,::1] x_m_ary, double[::1] x_m_sum, 
                              int k_max, int npixROI):
    """ Calculation of log likelihood

    Take the x_m and x_m_sum values, which represent the information for the
    non-poissonian templates, combine these with sum of the poissonian
    templates (PT_sum) and calculate the full log likelihood.

    Calculation for the likelihood proceeds via a recursion relation for the pk
    values.

    Returns:
        double log likelihood 

    """

    cdef Py_ssize_t p, k, n
    cdef double ll = 0.
    cdef double f0_ary, f1_ary
    cdef double[:] pk = np.zeros((k_max+1), dtype=DTYPE)

    # Loop over pixels
    for p in range(npixROI):
        # Define p_0 (pk[0]) and p_1 (pk[1])
        # Then the remaining p_k are determined recursively up to the value of 
        # k = data in that pixel
        f0_ary = -(pt_sum_compressed[p] + x_m_sum[p])
        f1_ary = (pt_sum_compressed[p] + x_m_ary[p,1])
        pk[0] = exp(f0_ary)
        pk[1] = pk[0] * f1_ary

        for k in range(2, data[p]+1):
            pk[k] = 0.
            for n in range(0, k-1):
                pk[k] += (k-n)/float(k)*x_m_ary[p,k-n]*pk[n]
            pk[k] += f1_ary*pk[k-1]/float(k)

        # Need prob > 0 to define a LL
        # If a very bad fit can get prob = 0, if so then penalise to avoid this
        # region
        if pk[data[p]] > 0:
            ll += log(pk[data[p]])
        else:
            ll += -10.1**10.

    return ll
