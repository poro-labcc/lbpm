#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 18:15:51 2026

@author: diogo
"""

import numpy as np
from  poretools import io
n = 8
r = 2

x, y, z = np.indices((n, n, n))

center = n // 2

dist2 = (x-center)**2 + (y-center)**2 + (z-center)**2

mat = np.ones((n, n, n), dtype=int)
mat[dist2 <= r**2] = 0

io.writeRAW( mat.astype( np.uint8) , "single_sphere.raw" ) 