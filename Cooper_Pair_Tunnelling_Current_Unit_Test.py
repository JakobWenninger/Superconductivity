#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:33:44 2020

@author: wenninger
"""
from argparse import ArgumentParser
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import os
import sys

sys.path.insert(0, '../Helper')

from Cooper_Pair_Tunnelling_Current import cooperPairTunnellingCurrent,magneticFluxQuantum
from Fundamental_BCS_Equations import cooperPair_Binding_Energy_over_T,cooperPair_Binding_Energy_0K
from plotxy import newfig,pltsettings,lbl,plot

parser = ArgumentParser()
parser.add_argument('-f', '--folder', action='store',default = 'Default_Folder', help='The folder in which the result is stored in.')
args = parser.parse_args()

directory = 'Cooper_Pair_Tunnelling_Current_Unit_Test/'+args.folder+'/'
if not os.path.exists(directory):
        os.makedirs(directory)
        
        
title = newfig('Flux_Dependency')   
fluxes = np.arange(0,15,.001) * magneticFluxQuantum
te = 4 #K, The actual temperature
tC = 9.2 # K for Nb
delta = cooperPair_Binding_Energy_over_T(t_over_Tc =te/tC,delta0 = cooperPair_Binding_Energy_0K(tC), tDebye =  276)[0] 
rN =10
current = cooperPairTunnellingCurrent(fluxes,delta,te,rN)
plt.plot(fluxes,current)
pltsettings(save=directory+title,fileformat='.pdf',disp = True,close=True, xlabel=lbl['Wb'],ylabel=lbl['A'], 
                xlim=None,ylim=None,title=None,legendColumns=1,skip_legend=True)

title = newfig('Normal_Resistance_Dependency')   
fluxes = 0
te = 4 #K, The actual temperature
delta = cooperPair_Binding_Energy_over_T(t_over_Tc =te/tC,delta0 = cooperPair_Binding_Energy_0K(tC), tDebye =  276)[0] #9.2 K for Nb
rN =np.arange(5,40,.01)
current = cooperPairTunnellingCurrent(fluxes,delta,te,rN)
plt.plot(rN,current)
pltsettings(save=directory+title,fileformat='.pdf',disp = True,close=False, xlabel=lbl['Rn'],ylabel=lbl['A'], 
                xlim=None,ylim=None,title=None,legendColumns=1,skip_legend=True)
