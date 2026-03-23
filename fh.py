#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 10:52:40 2026

@author: diogo
"""

import numpy as np
from poretools import io
from scipy import ndimage
import matplotlib.pylab as plt 

geo = io.readMHD( "sphere-pack-50c-uint8.mhd" )

def edt_1d(f):
    """
    Transformada 1D de Felzenszwalb-Huttenlocher.

    Entrada
    -------
    f : array 1D de floats
        f[q] = 0 nos sites de interesse
        f[q] = +inf nos demais

    Saída
    -----np.array( [ 0, 1e20, 1e20, 1e20, 0, 0 ] )
    g : array 1D de floats
        g[x] = min_q ((x - q)^2 + f[q])
    """
    n = len(f)
    g = np.empty(n, dtype=np.float64)

    v = np.empty(n, dtype=np.int64)      # índices das parábolas no envelope
    z = np.empty(n + 1, dtype=np.float64)  # fronteiras entre parábolas

    k = 0
    v[0] = 0
    z[0] = -np.inf
    z[1] = np.inf

    def intersection(q, vk):
        # posição s onde a parábola de q passa a vencer a de vk
        return ((f[q] + q * q) - (f[vk] + vk * vk)) / (2.0 * (q - vk))

    # construção do envelope inferior
    for q in range(1, n):
        
        # Encontra a interseção entre a paráb
        
        s = intersection(q, v[k])
        while s <= z[k]:
            k -= 1
            s = intersection(q, v[k])
        k += 1
        v[k] = q
        z[k] = s
        z[k + 1] = np.inf

    # avaliação do envelope
    k = 0
    for x in range(n):
        while z[k + 1] < x:
            k += 1
        dx = x - v[k]
        g[x] = dx * dx + f[v[k]]

    return g

def edt_3d(binary, target=1):
    """
    EDT exata 3D para um volume binário.

    Parâmetros
    ----------
    binary : ndarray 3D
        Volume binário.
    target : int
        Valor considerado como 'objeto' (distância até ele).

    Retorna
    -------
    dist2 : ndarray 3D
        Distância euclidiana ao quadrado.
    dist : ndarray 3D
        Distância euclidiana.
    """
    binary = np.asarray(binary)
    assert binary.ndim == 3

    inf = 1e20

    # f = 0 nos voxels-alvo, inf nos demais
    f = np.where(binary == target, 0.0, inf).astype(np.float64)

    nx, ny, nz = f.shape

    # Passo no eixo 0
    for y in range(ny):
        for z in range(nz):
            f[:, y, z] = edt_1d(f[:, y, z])

    # Passo no eixo 1
    for x in range(nx):
        for z in range(nz):
            f[x, :, z] = edt_1d(f[x, :, z])

    # Passo no eixo 2
    for x in range(nx):
        for y in range(ny):
            f[x, y, :] = edt_1d(f[x, y, :])

    dist2 = f
    dist = np.sqrt(dist2)
    return dist

edt_fh = edt_3d( geo , target = 0)
edt_scipy = ndimage.distance_transform_edt( geo )