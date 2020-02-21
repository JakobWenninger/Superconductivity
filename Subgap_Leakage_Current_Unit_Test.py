#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 13:48:07 2020

@author: wenninger
"""

from argparse import ArgumentParser
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import os
import sys

sys.path.insert(0, '../Helper')

from Subgap_Leakage_Current import subgapLeakageCurrent
from Fundamental_BCS_Equations import cooperPair_Binding_Energy_over_T,cooperPair_Binding_Energy_0K
from plotxy import newfig,pltsettings,lbl,plot

parser = ArgumentParser()
parser.add_argument('-f', '--folder', action='store',default = 'Default_Folder', help='The folder in which the result is stored in.')
args = parser.parse_args()

directory = 'Subgap_Leakage_Current_Unit_Test/'+args.folder+'/'
if not os.path.exists(directory):
        os.makedirs(directory)
        
tC = 9.2 # K for Nb
rN =13
bias=np.arange(0,3e-3,0.001e-3)

title = newfig('Temperature_Dependence')
temperatures=[3,4,5,6,7]
for te in temperatures:
    delta = cooperPair_Binding_Energy_over_T(t_over_Tc =te/tC,delta0 = cooperPair_Binding_Energy_0K(tC), tDebye =  276)[0] 
    current = subgapLeakageCurrent(v = bias,t = te,d = delta,rN = rN)
    plt.plot(bias*1e3,current*1e6, label = '%.1f K'%te )
pltsettings(save=directory+title,fileformat='.pdf',disp = True,close=False, xlabel=lbl['mV'],ylabel=lbl['uA'], 
                xlim=None,ylim=None,title=None,legendColumns=1,skip_legend=False)
    

title = newfig('Normal_Resistance_Dependence')
te=4
delta = cooperPair_Binding_Energy_over_T(t_over_Tc =te/tC,delta0 = cooperPair_Binding_Energy_0K(tC), tDebye =  276)[0] 
rNs = np.arange(10,15,1)
for rN in rNs:
    current = subgapLeakageCurrent(v = bias,t = te,d = delta,rN = rN)
    plt.plot(bias*1e3,current*1e6, label = '$%.1f\\,[\\Omega]$'%rN )
pltsettings(save=directory+title,fileformat='.pdf',disp = True,close=False, xlabel=lbl['mV'],ylabel=lbl['uA'], 
                xlim=None,ylim=None,title=None,legendColumns=1,skip_legend=False)
    
