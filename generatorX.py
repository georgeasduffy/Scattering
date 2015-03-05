# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 14:13:28 2015

@author: gaduffy2
"""

import matplotlib.pyplot as plt
from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, GammaPSD, BinnedPSD
from pytmatrix import orientation, radar, tmatrix_aux, refractive
import numpy as np
import scipy.io as sio 
from netCDF4 import Dataset
import pdb

def mix(a, b, d, ar):
    # a, d in cgs, use bruggeman equation
    m = a*d**b
    v = (1./6.)*ar*(d)**3*np.pi
    density = m/v
    n = 1.782
    kappa = complex(0, 7.302*10**-3)
    IceDens = .917
    e1 = n**2 - kappa**2
    e2 = 2*n*kappa
    B = e1 + e2
    f = density/IceDens
    if f > 1:
        return complex(n, 7.302*10**-3)
    else:
        mixed = (B*(1+2*f)-(2*f-2))/((2+f)+B*(1-f))
        e1out = np.real(mixed)
        e2out = np.imag(mixed)
        nout = np.sqrt((np.sqrt(e1out**2 +e2out**2)+e1out)/2)
        kout = np.sqrt((np.sqrt(e1out**2 +e2out**2)-e1out)/2)
        m =complex(nout, kout)
    return m
    
ar = np.linspace(0.2,1, num = 16)
b = np.linspace(1.8, 2.3, num = 10)
a = np.linspace(0.002, 0.012, num = 20)
d = np.logspace(-2, 2, num = 1000) #mm
wave = (8.5714, 21.4286)
backscatter = Dataset('scatsout_a.cdf', 'w')

AR = backscatter.createDimension('Area Ratio', len(ar.tolist()))
B = backscatter.createDimension('Mass Parameter b', len(b.tolist()))
A = backscatter.createDimension('Mass Parameter a', len(a.tolist()))
D = backscatter.createDimension('Maximum Diameter', len(d.tolist()))
WAVE = backscatter.createDimension('Wavelength', len(wave))
xsections = backscatter.createVariable('backscat', 'f4', ('Area Ratio', 'Mass Parameter b', 'Mass Parameter a', 'Maximum Diameter', 'Wavelength'))

for i in range(len(a)):
    print(a[i])
    for j in range(len(b)):
        print(b[j])
        for k in range(len(ar)):
            print(ar[k])
            for l in range(len(d)):
                for n in range(len(wave)):

                     if wave[n]/d[l] > 0.15:
                         dmod = d[l] * (ar[k]**(1./3.))
                         newm = mix(a[i], b[j], d[l]*0.1, ar[k])
                         oblate = Scatterer(radius = dmod/2., wavelength = wave[n], m = newm, axis_ratio = 1/ar[k])
                         xsection = 10*np.log10(radar.refl(oblate))
                         xsections[i,j,k,l,n] = 10*np.log10(radar.refl(oblate))
                     else:
                         xsections[i,j,k,l,n] = np.nan

backscatter.close()
