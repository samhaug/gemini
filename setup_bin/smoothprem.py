#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : smoothprem.py
Purpose : Make prem with smooth TZ ala Mike Thorne
Creation Date : 09-07-2017
Last Modified : Tue 11 Jul 2017 01:34:40 PM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt

def main():
    param = smooth_profile()
    param = np.array(smooth_profile())
    p = np.genfromtxt('/home/samhaug/work1/ScS_reverb_setup/models/prem.tvel',skip_header=2)

    plt.plot(6371-p[:,0],p[:,1])
    plt.plot(6371-p[:,0],p[:,2])
    plt.plot(6371-p[:,0],p[:,3])

    plt.plot(param[:,0],param[:,1],color='k',label='rho')
    plt.plot(param[:,0],param[:,2],color='r',label='vp')
    plt.plot(param[:,0],param[:,3],color='b',label='vs')
    plt.legend()
    plt.show()

def smooth_profile():
    re = 6371.
    param_list = []
    for dist in range(6372):
        x = dist/re
        if dist >= 6346.6:
            mu_p = 4.49100712
            rho_p = 3.38074821
            lam_p = 8.11061727
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 6346.6 and dist >= 6181.0:
            rho_p=2.691+.6924*x
            lam_p=4.1875+3.9382*x
            mu_p=2.1519+2.3481*x
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 6181.0 and dist >= 6121.0:
            R11 = 6181.0
            VS11 = 4.42997347
            VP11 = 8.00825250
            RHO11 = 3.36275081
            R12 = 6121.0
            VS12 = 4.66490000
            VP12 = 8.61666453
            RHO12 = 3.45368975
            m = (VS12-VS11)/(R12-R11)
            mu_p = m*(dist-R11)+VS11
            m = (VP12-VP11)/(R12-R11)
            lam_p = m*(dist-R11)+VP11
            m = (RHO12-RHO11)/(R12-R11)
            rho_p = m*(dist-R11)+RHO11
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 6121.0 and dist >= 6001.:
            rho_p=7.1089-3.8045*x
            lam_p=20.3926-12.2569*x
            mu_p=8.9496-4.4597*x
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 6001. and dist >= 5941.0:
            R11 = 6001.0
            VS11 = 4.74890000
            VP11 = 8.84752750
            RHO11 = 3.52534883
            R12 = 5941.0
            VS12 = 5.0200040
            VP12 = 9.28750292
            RHO12 = 3.76155793
            m = (VS12-VS11)/(R12-R11)
            mu_p = m*(dist-R11)+VS11
            m = (VP12-VP11)/(R12-R11)
            lam_p = m*(dist-R11)+VP11
            m = (RHO12-RHO11)/(R12-R11)
            rho_p = m*(dist-R11)+RHO11
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist >= 5771.0:
            rho_p=11.2494-8.0298*x
            lam_p=39.7027-32.6166*x
            mu_p=22.3512-18.5856*x
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 5771. and dist>=5731.0:
            rho_p=5.3197-1.4836*x
            lam_p=19.0957-9.8672*x
            mu_p=9.9839-4.9324*x
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 5731. and dist>=5671.0:
            R11 = 5731.0
            VS11 = 5.54698517
            VP11 = 10.21971143
            RHO11 = 3.98513532
            R12 = 5671.0
            VS12 = 6.03284432
            VP12 = 10.84473641
            RHO12 = 4.39943638
            m = (VS12-VS11)/(R12-R11)
            mu_p = m*(dist-R11)+VS11
            m = (VP12-VP11)/(R12-R11)
            lam_p = m*(dist-R11)+VP11
            m = (RHO12-RHO11)/(R12-R11)
            rho_p = m*(dist-R11)+RHO11
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 5671. and dist >= 5600:
            rho_p=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3
            lam_p=29.2766-23.6027*x+5.5242*x**2-2.5514*x**3
            mu_p=22.3459-17.2473*x-2.0834*x**2+0.9783*x**3
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 5600. and dist >= 3630:
            rho_p=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3
            lam_p=24.9520-40.4673*x+51.4832*x**2-26.6419*x**3
            mu_p=11.1671-13.7818*x+17.4575*x**2-9.2777*x**3
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 3630. and dist >= 3480:
            rho_p=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3
            lam_p=15.3891-5.3181*x+5.5242*x**2-2.5514*x**3
            mu_p=6.9254+1.4672*x-2.0834*x**2+0.9783*x**3
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 3480. and dist >= 1221.5:
            rho_p=12.5815-1.2638*x-3.6426*x**2-5.5281*x**3
            lam_p=11.0487-4.0362*x+4.8023*x**2-13.5732*x**3
            mu_p=0.
            param_list.append([dist,rho_p,lam_p,mu_p])
        elif dist < 1221.5:
            rho_p=13.0885-8.8381*x**2
            lam_p=11.2622-6.3640*x**2
            mu_p=3.6678-4.4475*x**2
            param_list.append([dist,rho_p,lam_p,mu_p])
    return param_list

main()
