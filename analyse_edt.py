#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 20:06:39 2026

@author: diogo
"""

import numpy as np
from poretools import io
from scipy import ndimage
import matplotlib.pylab as plt 

geo = io.readMHD( "single_sphere.mhd" )

#edt_zabot = np.sqrt( io.readRAW(  "edt-squared-50c-int32.raw"  , (50,50,50) , data = np.int32 ) )
edt_scipy = ndimage.distance_transform_edt( geo  )**2
edt_lbpm =  io.readRAW(  "distance.raw", tuple( np.array( geo.shape ) + 2 ) , data = np.float64 )[1:-1,1:-1,1:-1]

fig, axes = plt.subplots(1, 3)  # 1 row, 2 columns

z = 4

# axes[0].imshow(edt_scipy[z] )
# axes[0].set_title("Scipy")

# diff = np.abs( edt_scipy - edt_lbpm[z] )
# axes[1].imshow( edt_scipy[z] - edt_lbpm[z] *geo[z])
# axes[1].set_title("Difference")

# axes[2].imshow( edt_lbpm[z] )
# axes[2].set_title("LBPM")
