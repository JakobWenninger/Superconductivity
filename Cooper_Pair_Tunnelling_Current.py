#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:14:52 2020

@author: wenninger
"""
import numpy as np
import scipy.constants as const
from scipy.special import jv

magneticFluxQuantum = np.divide(const.h,2*const.e)

def cooperPairTunnellingCurrent(magneticFlux,delta,te,rN):
    '''This function computes the Cooper pair tunnelling current through an SIS junction.
    Formula from Tinkham Superconductivity [2004] equation 6.10 with the addition 
    of the magnetic flux dependency from Garret [2018].

    Only tested/intended for one parameter to be an array.
    
    inputs
    ------
    magneticFlux: float or 1d array
        The magnetic flux applied through the SIS junction.
    delta: float
        The Cooper pair binding energy at the temperature te.
    te: float
        The physical temperature of the junction.
    rN: float or 1d array
        The normal resistance of the junction.
        
    returns
    -------
    float or array:
        The Cooper pair tunnelling current.
    '''
    def cooperPairTunnellingCurrent_without_Magnetic_Flux(delta,te,rN):
        '''This function computes the Cooper pair tunnelling current 
            through an SIS junction without a magnetic flux present.
        Formula from Tinkham Superconductivity [2004] equation 6.10.
        
        
        inputs
        ------
        delta: float
            The Cooper pair binding energy at the temperature te.
        te: float
            The physical temperature of the junction.
        rN: float 
            The normal resistance of the junction.
            
        returns
        -------
        float or array:
            The Cooper pair tunnelling current with no magnetic flux present.
        '''
        #np.reciprocal can't deal with integers
        return np.multiply(np.multiply(np.reciprocal(rN*1.),
                                           np.divide(np.multiply(np.pi,delta),2*const.e)),
                               np.tanh(np.divide(delta,np.multiply(2*const.k,te))))
    x = np.divide(magneticFlux,magneticFluxQuantum)
    with np.errstate(divide='ignore', invalid='ignore'):# This is necessary for evaluation of the raw data
       fluxcoefficent = 2* np.divide(jv(1,x),x)
    current = cooperPairTunnellingCurrent_without_Magnetic_Flux(delta,te,rN)
    if isinstance(x,float):
        if x !=0:
            current = np.multiply(current,fluxcoefficent)
    else:        
        current = np.full(x.shape,current)
        current[np.nonzero(x)] = np.multiply(current[np.nonzero(x)],fluxcoefficent[np.nonzero(x)])
    return current
    
    