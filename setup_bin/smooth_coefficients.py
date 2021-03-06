#!/home/samhaug/anaconda2/bin/python

import numpy as np
from matplotlib import pyplot as plt

def main():
    smooth_dict = make_smooth_dict()
    nlay,pco = get_header(smooth_dict)
    write_model('SMOOTHPREM',smooth_dict,nlay,pco)

def write_model(fname,sd,nlay,pco):
    with open(fname,'w') as f:
        f.write('#--------------------------------------------------------------------\n')
        f.write('#                      Earth model SMOOTHPREM                        \n')
        f.write('#       ( isotropic PREM with no ocean and smooth transition zone)   \n')
        f.write('#                                                                    \n')
        f.write('# Quality factors are taken at 1 Hz !!                               \n')
        f.write('# The ocean is removed by enlarging the thickness its underlying solid\n')
        f.write('# layer by the ocean thickness.                                      \n')
        f.write('#------------------------------------------------------------         \n')
        f.write('#                                                                     \n')
        f.write('# +++++++ Do not insert any comments below this line !! +++++++++++++++\n')
        f.write('%2d\n'%nlay)
        for ii in pco[::-1]:
            f.write('%2d'%ii)
        f.write('\n1.0\n0\n')
        f.write('-'*74+'\n')
        f.write('Radius Density  V-Pv  V-Ph  V-Sv  V-Sh  Qmu  Qk  Eta  \n')
        f.write('-'*74+'\n')
        for keys in reversed(sorted(sd.keys())):
            if type(sd[keys]['rho']) == float:
                f.write('%6.1f%9.4f%9.4f%9.4f%9.4f%9.4f%7.1f%9.1f%5.1f\n'%(
                         sd[keys]['h'],sd[keys]['rho'],sd[keys]['vp'],
                         sd[keys]['vp'],sd[keys]['vs'],sd[keys]['vs'],
                         sd[keys]['qmu'],sd[keys]['qk'],sd[keys]['eta']))
                f.write('     \n')
            elif type(sd[keys]['rho']) == tuple:
                for idx in range(len(sd[keys]['rho'])):
                    if idx == 0:
                        f.write('%6.1f%9.4f%9.4f%9.4f%9.4f%9.4f%7.1f%9.1f%5.1f\n'%(
                         sd[keys]['h'],sd[keys]['rho'][idx],sd[keys]['vp'][idx],
                         sd[keys]['vp'][idx],sd[keys]['vs'][idx],sd[keys]['vs'][idx],
                         sd[keys]['qmu'],sd[keys]['qk'],sd[keys]['eta'][idx]))
                    else:
                        f.write('%15.4f%9.4f%9.4f%9.4f%9.4f%21.1f\n'%(
                         sd[keys]['rho'][idx],sd[keys]['vp'][idx],
                         sd[keys]['vp'][idx],sd[keys]['vs'][idx],sd[keys]['vs'][idx],
                         sd[keys]['eta'][idx]))
                f.write('    \n')

        f.write('     \n')
        f.write('%6.1f\n'%6371.0)

def get_header(smooth_dict):
    nlay = len(smooth_dict.keys())
    pco = []
    for keys in smooth_dict:
        if type(smooth_dict[keys]['rho']) == float:
                pco.append(1)
        else:
            pco.append(len(smooth_dict[keys]['rho']))
    return nlay,pco

def make_smooth_dict():
    smooth_dict = {}

    smooth_dict[1] = {}
    smooth_dict[1]['h'] = (6346.6)
    smooth_dict[1]['vs'] = (4.49100712)
    smooth_dict[1]['rho'] = (3.38074821)
    smooth_dict[1]['vp'] = (8.11061727)
    smooth_dict[1]['qmu'] = (600.0)
    smooth_dict[1]['qk'] = (57823.0)
    smooth_dict[1]['eta'] = (1.0)

    smooth_dict[2] = {}
    smooth_dict[2]['h'] = 6181.0
    smooth_dict[2]['vs'] = (2.1519,2.3481)
    smooth_dict[2]['rho'] = (2.691,0.6924)
    smooth_dict[2]['vp'] = (4.1875,3.9382)
    smooth_dict[2]['qmu'] = (600.0)
    smooth_dict[2]['qk'] = (57823.0)
    smooth_dict[2]['eta'] = (1.0,0.0)

    smooth_dict[3] = {}
    R11 = 6181.0
    VS11 = 4.42997347
    VP11 = 8.00825250
    RHO11 = 3.36275081
    R12 = 6121.0
    VS12 = 4.66490000
    VP12 = 8.61666453
    RHO12 = 3.45368975
    smooth_dict[3]['h'] = R12

    m = ((VS12-VS11)/(R12/6371.-R11/6371.))
    VSB = -m*R12/6371.+VS12
    smooth_dict[3]['vs'] = (VSB,m)

    m = ((RHO12-RHO11)/(R12/6371.-R11/6371.))
    RHOB = -m*R12/6371.+RHO12
    smooth_dict[3]['rho'] = (RHOB,m)

    m = ((VP12-VP11)/(R12/6371.-R11/6371.))
    VPB = -m*R12/6371.+VP12
    smooth_dict[3]['vp'] = (VPB,m)

    smooth_dict[3]['qmu'] = (600.0)
    smooth_dict[3]['qk'] = (57823.0)
    smooth_dict[3]['eta'] = (1.0,0.0)

    smooth_dict[4] = {}
    smooth_dict[4]['h'] = 6001.
    smooth_dict[4]['rho'] = (7.1089,-3.8045)
    smooth_dict[4]['vp'] = (20.3926,-12.2569)
    smooth_dict[4]['vs'] = (8.9496,-4.4597)
    smooth_dict[4]['qmu'] = (80.0)
    smooth_dict[4]['qk'] = (57823.0)
    smooth_dict[4]['eta'] = (1.0,0.0)

    smooth_dict[5] = {}
    R11 = 6001.0
    VS11 = 4.74890000
    VP11 = 8.84752750
    RHO11 = 3.52534883
    R12 = 5941.0
    VS12 = 5.0200040
    VP12 = 9.28750292
    RHO12 = 3.76155793
    smooth_dict[5]['h'] = R12

    m = ((RHO12-RHO11)/(R12/6371.-R11/6371.))
    RHOB = -m*R12/6371.+RHO12
    smooth_dict[5]['rho'] = (RHOB,m)

    m = ((VP12-VP11)/(R12/6371.-R11/6371.))
    VPB = -m*R12/6371.+VP12
    smooth_dict[5]['vp'] = (VPB,m)

    m = ((VS12-VS11)/(R12/6371.-R11/6371.))
    VSB = -m*R12/6371.+VS12
    smooth_dict[5]['vs'] = (VSB,m)

    smooth_dict[5]['qmu'] = (80.0)
    smooth_dict[5]['qk'] = (57823.0)
    smooth_dict[5]['eta'] = (1.0,0.0)

    smooth_dict[6]={}
    smooth_dict[6]['h'] = 5771.
    smooth_dict[6]['rho'] = (11.2494,-8.0298)
    smooth_dict[6]['vp'] = (39.7027,-32.6166)
    smooth_dict[6]['vs'] = (22.3512,-18.5856)
    smooth_dict[6]['qmu'] = (80.0)
    smooth_dict[6]['qk'] = (57823.0)
    smooth_dict[6]['eta'] = (1.0,0.0)

    smooth_dict[7]={}
    smooth_dict[7]['h'] = 5731.
    smooth_dict[7]['rho'] = (5.3197,-1.4836)
    smooth_dict[7]['vp'] = (19.0957,-9.8672)
    smooth_dict[7]['vs'] = (9.9839,-4.9324)
    smooth_dict[7]['qmu'] = (143.0)
    smooth_dict[7]['qk'] = (57823.0)
    smooth_dict[7]['eta'] = (1.0,0.0)

    smooth_dict[8]={}
    R11 = 5731.0
    VS11 = 5.54698517
    VP11 = 10.21971143
    RHO11 = 3.98513532
    R12 = 5671.0
    VS12 = 6.03284432
    VP12 = 10.84473641
    RHO12 = 4.39943638
    smooth_dict[8]['h'] = R12

    m = ((RHO12-RHO11)/(R12/6371.-R11/6371.))
    RHOB = -m*R12/6371.+RHO12
    smooth_dict[8]['rho'] = (RHOB,m)

    m = ((VP12-VP11)/(R12/6371.-R11/6371.))
    VPB = -m*R12/6371.+VP12
    smooth_dict[8]['vp'] = (VPB,m)

    m = ((VS12-VS11)/(R12/6371.-R11/6371.))
    VSB = -m*R12/6371.+VS12
    smooth_dict[8]['vs'] = (VSB,m)

    smooth_dict[8]['qmu'] = (143.0)
    smooth_dict[8]['qk'] = (57823.0)
    smooth_dict[8]['eta'] = (1.0,0.0)

    smooth_dict[9]={}
    smooth_dict[9]['h'] = 5600.
    smooth_dict[9]['rho'] = (7.9565,-6.4761,5.5283,-3.0807)
    smooth_dict[9]['vp'] = (29.2766,-23.6027,5.5242,-2.5514)
    smooth_dict[9]['vs'] = (22.3459,-17.2473,-2.0834,0.9783)
    smooth_dict[9]['qmu'] = (312.0)
    smooth_dict[9]['qk'] = (57823.0)
    smooth_dict[9]['eta'] = (1.0,0.0,0.0,0.0)

    smooth_dict[10]={}
    smooth_dict[10]['h'] = 3630.
    smooth_dict[10]['rho'] = (7.9565,-6.4761,5.5283,-3.0807)
    smooth_dict[10]['vp'] = (24.9520,-40.4673,51.4832,-26.6419)
    smooth_dict[10]['vs'] = (11.1671,-13.7818,17.4575,-9.2777)
    smooth_dict[10]['qmu'] = (312.0)
    smooth_dict[10]['qk'] = (57823.0)
    smooth_dict[10]['eta'] = (1.0,0.0,0.0,0.0)

    smooth_dict[11]={}
    smooth_dict[11]['h'] = 3480.
    smooth_dict[11]['rho'] =(7.9565,-6.4761,5.5283,-3.0807)
    smooth_dict[11]['vp'] = (15.3891,-5.3181,5.5242,-2.5514)
    smooth_dict[11]['vs'] = (6.9254,1.4672,-2.0834,.9783)
    smooth_dict[11]['qmu'] = (312.0)
    smooth_dict[11]['qk'] = (57823.0)
    smooth_dict[11]['eta'] = (1.0,0.0,0.0,0.0)

    smooth_dict[12]={}
    smooth_dict[12]['h'] = 1221.5
    smooth_dict[12]['rho'] = (12.5815,-1.2638,-3.6426,-5.5281)
    smooth_dict[12]['vp'] = (11.0487,-4.0362,4.8023,-13.5732)
    smooth_dict[12]['vs'] = (0.0,0.0,0.0,0.0)
    smooth_dict[12]['qmu'] = (-1.0)
    smooth_dict[12]['qk'] = (57823.0)
    smooth_dict[12]['eta'] = (1.0,0.0,0.0,0.0)

    smooth_dict[13]={}
    smooth_dict[13]['h'] = 0.0
    smooth_dict[13]['rho'] = (13.0885,0.0,-8.8381)
    smooth_dict[13]['vp'] = (11.2622,0.0,-6.3640)
    smooth_dict[13]['vs'] = (3.6678,0.0,-4.4475)
    smooth_dict[13]['qmu'] = (84.6)
    smooth_dict[13]['qk'] = (1327.7)
    smooth_dict[13]['eta'] = (1.0,0.0,0.0)

    return smooth_dict

main()
