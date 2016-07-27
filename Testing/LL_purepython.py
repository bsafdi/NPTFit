#------------------------------------------------
# NPTF Log Likelihood
#------------------------------------------------
## 
## Calculation of Log Likelihood for a template
## fit involving non-poissonian templates
##
#------------------------------------------------

import numpy as np
from mpmath import mp
gammainc_mem = mp.gammainc

sp_max = 5000
def log_factorial(k):
    return np.sum(np.log(np.arange(1., k + 1.)))
log_factorial_ary = np.vectorize(log_factorial)(np.arange(sp_max))

def gammainc_up_fct_ary_log(m_max, z, a):
    gi_list = np.zeros(m_max + 1)
    gi_list[0] = gammainc_mem(z, a)

    i_array = np.arange(1., m_max + 1.)
    for i in np.arange(1, m_max + 1):
        gi_list[i] = (i_array[i-1] - 1. + z)*gi_list[i-1]/i + np.exp((i_array[i-1] - 1. + z)*np.log(a) - a - log_factorial_ary[i])
    return np.asarray(gi_list)

def gammainc_lo_fct_ary_back(m_max, z, b):
    gi_list = np.zeros(m_max + 1)
    gi_list[m_max] = gammainc_mem(z + m_max, 0, b) / mp.factorial(m_max)

    i_array = np.arange(m_max, 0, -1)

    for i in np.arange(m_max, 0, -1):
        gi_list[i-1] = (gi_list[i] + np.exp((i_array[-i] - 1. + np.float128(z))*np.log(np.float128(b)) - np.float128(b) - log_factorial_ary[i]))*i_array[-i]/(i_array[-i] - 1. + np.float128(z))
    return np.asarray(gi_list)

def log_like(xbg_PSF_compressed, theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data):
    
    # Initialisation
    k_max = np.max(data)+1

    A = theta[0]
    n1 = theta[1]
    n2 = theta[2]
    n3 = theta[3]
    Sb1 = theta[4]
    Sb2 = theta[5]

    npixROI = len(xbg_PSF_compressed)

    f2 = 0.
    df_rho_div_f2 = 0.

    x_m_ary = np.zeros(shape=(npixROI,k_max+1))
    x_m_sum = np.zeros(npixROI)

    g0_ary_f = np.zeros(k_max + 1)
    g1_ary_f = np.zeros(k_max + 1)
    g2_ary_f = np.zeros(k_max + 1)

    # PS Calcs
    term1 = 0.
    term2 = 0.
    second_2_a = 0.
    second_2_b = 0.
    second_2_c = 0.
    second_2_d = 0.
    second_1_alpha = 0.
    second_1_a = 0.
    second_1_b = 0.
    second_1_beta = 0.
    second_1_c = 0.
    second_1_d = 0.

    for f_index in range(len(f_ary)):
        f2 = f_ary[f_index]
        df_rho_div_f2 = df_rho_div_f_ary[f_index]
        g0_ary_f = gammainc_up_fct_ary_log(k_max, 1. - n1, Sb2 * f2)
        g1_ary_f = gammainc_up_fct_ary_log(k_max, 1. - n2, Sb1 * f2) - gammainc_up_fct_ary_log(k_max, 1. - n2, Sb2 * f2)
        g2_ary_f = gammainc_lo_fct_ary_back(k_max, 1. - n3, Sb1 * f2)
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


    # Now calculate nu values
    nu_ary = np.zeros(k_max + 1)
    nu_mat = np.zeros(k_max+1)
    ll = 0.

    for p in range(npixROI):
        f0_ary = -(xbg_PSF_compressed[p] + x_m_sum[p])
        f1_ary = (xbg_PSF_compressed[p] + x_m_ary[p,1])
        nu_mat[0] = np.exp(f0_ary)
        nu_mat[1] = nu_mat[0] * f1_ary

        for k in range(2,data[p]+1):
            nu_mat[k] = 0.
            for n in range(0, k - 1):
                nu_mat[k] += (k-n)/ float(k) * x_m_ary[p,k-n] * nu_mat[n]
            nu_mat[k] += f1_ary * nu_mat[k-1] / float(k)

        if nu_mat[data[p]] > 0:
            ll+=np.log( nu_mat[data[p]])
        else:
            ll+= -10.1**10.

    return ll
