#------------------------------------------------------------------------------
# x_m and x_m_sum arrays for arbitary number of breaks
#------------------------------------------------------------------------------
##
## Calculates x_m and x_m_sum arrays for arbitary number of breaks.
## Uses dedicated functions for 1, 2 and 3 breaks as well as n > 3 breaks.
##
#------------------------------------------------------------------------------

import numpy as np
from mpmath import mp

import incgamma_fct as ig

def return_xs(theta, f_ary, df_rho_div_f_ary, PS_dist_compressed, data):
    """
    Returns arrays of x_m and x_m_sum
    """
    n_break = (len(theta) - 2)/2 # TODO: handle exception in length of theta
    print "The number of breaks is", n_break
    
    if n_break == 1:

        k_max = int(max(data) + 1)

        A = float(theta[0])
        n1 = float(theta[1])
        n2 = float(theta[2])
        Sb = float(theta[3])
        
        npixROI = len(PS_dist_compressed)

        x_m_ary = np.zeros((npixROI,k_max + 1))
        x_m_sum = np.zeros(npixROI)

        g1_ary_f = np.zeros(k_max + 1)
        g2_ary_f = np.zeros(k_max + 1)

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

                    x_m_ary_f = A  * (pref1_x_m_ary  * g1_ary_f[k] + pref2_x_m_ary * g2_ary_f[k] ) * PS_dist_compressed[p]
                    x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f

        x_m_sum -= x_m_ary[:,0]

    elif n_break == 2:

        k_max = int(np.max(data) + 1)

        A = float(theta[0])
        n1 = float(theta[1])
        n2 = float(theta[2])
        n3 = float(theta[3])
        Sb1 = float(theta[4])
        Sb2 = float(theta[5])

        Sb1, Sb2 = Sb2, Sb1 # Cheeky variable swap
        A = A * (Sb2/Sb1)**n2  # Cheeky rescaling

        npixROI = len(PS_dist_compressed)

        x_m_ary = np.zeros(shape=(npixROI,k_max+1))
        x_m_sum = np.zeros(npixROI)

        g0_ary_f = np.zeros(k_max + 1)
        g1_ary_f = np.zeros(k_max + 1)
        g2_ary_f = np.zeros(k_max + 1)

        for f_index in range(len(f_ary)):

            f2 = float(f_ary[f_index])
            df_rho_div_f2 = df_rho_div_f_ary[f_index]
            
            g0_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n1, Sb2 * f2)
            g1_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n2, Sb1 * f2) - ig.incgamma_up_fct_ary(k_max, 1. - n2, Sb2 * f2)
            g2_ary_f = ig.incgamma_lo_fct_ary(k_max, 1. - n3, Sb1 * f2)

            term1 = (A * Sb1 * f2) \
                                 * (1./(n2-1.) + 1./(1.-n3) - pow(Sb1 / Sb2, n2-1.)/(n2-1.) \
                                    - (pow(Sb1 * f2, n2-1.) * g1_ary_f[0] + pow(Sb1 * f2, n3-1.) * g2_ary_f[0])   + pow( Sb1 / Sb2, n2 - 1.)*( 1./(n1-1.) - pow(Sb2*f2, n1 - 1.)*g0_ary_f[0] )  )

            pref0_x_m_ary = pow(Sb2 * f2, n1)* pow(Sb1 / Sb2 , n2)
            pref1_x_m_ary =  pow(Sb1 * f2, n2)
            pref2_x_m_ary = pow(Sb1 * f2, n3)

            for p in range(npixROI):
                x_m_sum_f = term1 * PS_dist_compressed[p]
                x_m_sum[p] += df_rho_div_f2*x_m_sum_f

                for k in range(data[p]+1):
                    x_m_ary_f = A * (pref0_x_m_ary * g0_ary_f[k] + pref1_x_m_ary  * g1_ary_f[k] + pref2_x_m_ary * g2_ary_f[k]) * PS_dist_compressed[p]
                    x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f

    elif n_break == 3:

        k_max = int(np.max(data) + 1)

        A = float(theta[0])
        n1 = float(theta[1])
        n2 = float(theta[2])
        n3 = float(theta[3])
        n4 = float(theta[4])
        Sb1 = float(theta[5])
        Sb2 = float(theta[6])
        Sb3 = float(theta[7])

        Sb1, Sb3 = Sb3, Sb1 # Cheeky variable swap
        A = A * (Sb2/Sb1)**n3 * (Sb3/Sb2)**n2 # Cheeky rescaling

        npixROI = len(PS_dist_compressed)

        x_m_ary = np.zeros((npixROI,k_max + 1))
        x_m_sum = np.zeros(npixROI)

        g0_ary_f = np.zeros(k_max + 1)
        g1_ary_f = np.zeros(k_max + 1)
        g2_ary_f = np.zeros(k_max + 1)
        g3_ary_f = np.zeros(k_max + 1)

        for f_index in range(len(f_ary)):
            f2 = float(f_ary[f_index])
            df_rho_div_f2 = df_rho_div_f_ary[f_index]

            g0_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n1, Sb3 * f2) 
            g1_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n2, Sb2 * f2) - ig.incgamma_up_fct_ary(k_max, 1. - n2, Sb3 * f2)
            g2_ary_f = ig.incgamma_up_fct_ary(k_max, 1. - n3, Sb1 * f2) - ig.incgamma_up_fct_ary(k_max, 1. - n3, Sb2 * f2)
            g3_ary_f = ig.incgamma_lo_fct_ary(k_max, 1. - n4, Sb1 * f2)

            kappa = pow(Sb2 / Sb1,1-n3)
            kappa2 = kappa * pow(Sb3 / Sb2,1-n2) 
            term1 = (A * Sb1 * f2) \
                                 * (kappa/(n2-1.)+1./(n3-1.) + 1./(1.-n4) - pow(Sb1 / Sb2, n3-1.)/(n3-1.) - kappa * pow(Sb2 / Sb3, n2-1.)/(n2-1.) \
                                    - (kappa*pow(Sb2 * f2, n2-1.) * g1_ary_f[0]+pow(Sb1 * f2, n3-1.) * g2_ary_f[0] + pow(Sb1 * f2, n4-1.) * g3_ary_f[0])   + kappa2*( 1./(n1-1.) - pow(Sb3*f2, n1 - 1.)*g0_ary_f[0] )  )

            pref0_x_m_ary = pow(Sb3 * f2, n1)* pow(Sb1 / Sb2 , n3) * pow(Sb2 / Sb3 , n2)
            pref1_x_m_ary =  pow(Sb2 * f2, n2)* pow(Sb1 / Sb2 , n3)
            pref2_x_m_ary = pow(Sb1 * f2, n3)
            pref3_x_m_ary = pow(Sb1 * f2, n4)

            for p in range(npixROI):
                x_m_sum_f = term1 * PS_dist_compressed[p]
                x_m_sum[p] += df_rho_div_f2*x_m_sum_f

                for k in range(data[p]+1):
                    x_m_ary_f = A * (pref0_x_m_ary * g0_ary_f[k] + pref1_x_m_ary  * g1_ary_f[k] + pref2_x_m_ary * g2_ary_f[k]+ pref3_x_m_ary * g3_ary_f[k]) * PS_dist_compressed[p]
                    x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f

    else: # arbitrary number of breaks (application for > 3)

        k_max = int(np.max(data) + 1)

        A = float(theta[0])
        n_ary = [float(n_val) for n_val in theta[1:n_break + 2]]
        Sb_ary = [float(n_val) for n_val in theta[n_break + 2:2*n_break + 2]]
        npixROI = len(PS_dist_compressed)

        x_m_ary = np.zeros((npixROI,k_max + 1))
        x_m_sum = np.zeros(npixROI)

        g_ary_f_ary = np.zeros((n_break + 1, k_max + 1))
        pref_x_m_ary = np.zeros(n_break + 1)

        first_x_m_sum_ary = np.zeros(n_break + 1)
        second_x_m_sum_ary = np.zeros(n_break + 1)


        for f_index in range(len(f_ary)):

            f2 = float(f_ary[f_index])
            df_rho_div_f2 = df_rho_div_f_ary[f_index]

            g_ary_f_ary[0] = ig.incgamma_up_fct_ary(k_max, 1. - n_ary[0], Sb_ary[0] * f2) 
            g_ary_f_ary[1:n_break] = [ig.incgamma_up_fct_ary(k_max, 1. - n_ary[i], Sb_ary[i] * f2) - ig.incgamma_up_fct_ary(k_max, 1. - n_ary[i], Sb_ary[i-1] * f2) for i in range(1,n_break)] 
            g_ary_f_ary[n_break] = ig.incgamma_lo_fct_ary(k_max, 1. - n_ary[n_break], Sb_ary[n_break-1] * f2) 


            # # This is as in Nick's notes
            # pref_x_m_ary[0] = 1.
            # pref_x_m_ary[1:n_break] = [np.prod([pow(Sb_ary[j-1] * f2 , n_ary[j] - n_ary[j-1]) for j in range(1, k + 1)]) for k in range(1,n_break)] # make sure np.prod is okay precision-wise!
            # pref_x_m_ary[n_break] = np.prod([pow(Sb_ary[j-1] * f2, n_ary[j] - n_ary[j-1]) for j in range(1, n_break + 1)])

            # second_x_m_sum_ary = [1/(1 - n_ary[i]) for i in range(n_break+1)]

            first_x_m_sum_ary[0] = -1. 
            first_x_m_sum_ary[1:n_break] = [-np.prod([pow(Sb_ary[j] / Sb_ary[j - 1], - n_ary[j]) for j in range(1, i)]) * pow(Sb_ary[i-1], n_ary[i]) * (-pow(Sb_ary[i-1], 1-n_ary[i])+pow(Sb_ary[i], 1-n_ary[i])) / Sb_ary[0] for i in range(1,n_break)]
            first_x_m_sum_ary[n_break] = np.prod([pow(Sb_ary[j] / Sb_ary[j - 1], - n_ary[j]) for j in range(1, n_break)]) * Sb_ary[n_break-1] / Sb_ary[0]

            print first_x_m_sum_ary
            for p in range(npixROI):
                
                x_m_sum_f = A * Sb_ary[0] * f2 * np.dot(first_x_m_sum_ary, second_x_m_sum_ary) * PS_dist_compressed[p]
                x_m_sum[p] += df_rho_div_f2 * x_m_sum_f
                
                for k in range(data[p]+1):

                    x_m_ary_f = A * pow(Sb_ary[0] * f2, n_ary[0]) * np.dot(pref_x_m_ary, np.transpose(g_ary_f_ary)[k]) * PS_dist_compressed[p]
                    x_m_ary[p,k] += df_rho_div_f2*x_m_ary_f

        x_m_sum -= x_m_ary[:,0]

    return np.asarray(x_m_ary), np.asarray(x_m_sum)
