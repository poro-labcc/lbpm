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

sl = slice(0,5),slice(45,None),slice(45,None)

geo = io.readMHD( "sphere-pack-50c-uint8.mhd" )

#edt_zabot = np.sqrt( io.readRAW(  "edt-squared-50c-int32.raw"  , (50,50,50) , data = np.int32 ) )
edt_scipy = ndimage.distance_transform_edt( geo  )**2
edt_lbpm =  io.readRAW(  "distance.raw", tuple( np.array( geo.shape ) + 2 ) , data = np.int32 )[1:-1,1:-1,1:-1]

fig, axes = plt.subplots(1, 3)  # 1 row, 2 columns

z =1

axes[0].imshow(edt_scipy[sl][z] )
axes[0].set_title("Scipy")

axes[1].imshow(  50*np.round( np.abs(edt_scipy[sl][z] - edt_lbpm[sl][z])) )
axes[1].set_title("Difference")

axes[2].imshow( edt_lbpm[sl][z] )
axes[2].set_title("LBPM")



geo = io.readMHD( "sphere-pack-50c-uint8.mhd" )
io.writeRAW( geo[sl], "slice.raw")
