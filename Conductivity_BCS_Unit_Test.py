#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:48:03 2020

@author: wenninger
"""
import matplotlib.pylab as plt
import numpy as np
import os
import scipy.constants as const
import sys

sys.path.insert(1, '../Helper')

from Fundamental_BCS_Equations import conductivity_BCS
from plotxy import plotcomplex,pltsettings,newfig,lbl

directory = 'Conductivity_BCS_Unit_Test/' 
if not os.path.exists(directory):
        os.makedirs(directory)
        
title = newfig('Nb_NbTiN')
freqrange = np.arange(100,2000,1)
deltaNb =1.4e-3*const.e
deltaNbTiN =2.2e-3*const.e
#Reproduction of Kooi's thesis figure 4.21
nb = []
nbTiN = []
for f in freqrange:
    nb.append(conductivity_BCS(freq=f*1e9,te=4.2,delta=deltaNb))
    nbTiN.append(conductivity_BCS(freq=f*1e9,te=4.2,delta=deltaNbTiN))
nb = np.array(nb)
nbTiN = np.array(nbTiN)
plt.plot(freqrange,nb.real, label=('$\\sigma_{1,\\mathrm{Nb}}$'))
plt.plot(freqrange,nb.imag, label=('$\\sigma_{2,\\mathrm{Nb}}$'))
plt.plot(freqrange,nbTiN.real, label=('$\\sigma_{1,\\mathrm{NbTiN}}$'))
plt.plot(freqrange,nbTiN.imag, label=('$\\sigma_{2,\\mathrm{NbTiN}}$'))
plt.vlines( 2e-9*deltaNb/const.h,0,10,alpha=1,linestyles='dashed')#,label='$f_\\mathrm{gap,Nb}$')  
plt.vlines( 2e-9*deltaNbTiN/const.h,0,10,alpha=1,linestyles='dashed')#,label='$f_\\mathrm{gap,NbTiN}$')  
pltsettings(directory+title,xlabel=lbl['GHz'],xlim=[100,2000],ylim = [1e-3,10] ,ylabel=lbl['normalisedSigma'],yscale='log',close=True,legendColumns=1,fileformat='.pdf') 

title = newfig('Nb_Temperature_Dependence')
freqrange = np.arange(400,1000,1)
temperatures = [3,4,6]
for t in temperatures:
    nbtemp = []
    for f in freqrange:
        nbtemp.append(conductivity_BCS(freq=f*1e9,te=t,delta=deltaNb))
    nbtemp = np.array(nbtemp)
    plt.plot(freqrange,nbtemp.real, label=('$\\sigma_{1} %.1f \\mathrm{K}$'%t))
    plt.plot(freqrange,nbtemp.imag, label=('$\\sigma_{2} %.1f \\mathrm{K}$'%t))
pltsettings(directory+title,xlabel=lbl['GHz'],xlim=[400,1000],ylim = [1e-3,3] ,ylabel=lbl['normalisedSigma'],yscale='log',close=True,legendColumns=1,fileformat='.pdf') 


title = newfig('Nb_Real_Temperature_Dependence')
freqrange = np.arange(400,900,1)
temperatures = [1.5,2,3,4,6,8]
for t in temperatures:
    nbtemp = []
    for f in freqrange:
        nbtemp.append(conductivity_BCS(freq=f*1e9,te=t,delta=deltaNb))
    nbtemp = np.array(nbtemp)
    plt.plot(freqrange,nbtemp.real, label=('$\\sigma_{1} (%.1f \\mathrm{K}$)'%t))
pltsettings(directory+title,xlabel=lbl['GHz'],xlim=[400,900],ylim = [1e-6,3] ,ylabel=lbl['normalisedSigma'],yscale='log',close=True,legendColumns=1,fileformat='.pdf') 

