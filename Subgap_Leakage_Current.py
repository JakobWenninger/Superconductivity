#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 13:34:00 2020

@author: wenninger
"""
import numpy as np
import scipy.constants as const
from scipy.special import kn


def subgapLeakageCurrent(v,t,d,rN):
    '''This function computes the subgap leakage current following the equation of Van Duzer and Turner [1981].
    The equation assumes the superconductors to be identical and that the temperature is well below the critical temperature.
    The equation is only valid below the gap voltage, so to say in the subgap region.
    
    inputs
    ------
    v: float or 1d array
        The bias voltage of the SIS junction. The result is only valid as the bias voltage is below the gap voltage.
        Note that the unit is V and not mV.
    t: float
        The temperature of the 
    d: flaot
        The Cooper pair binding energy delta of the superconducting material.
    rN: float
        The normal resistance.
    
    returns
    -------
    float or 1d array  
        The subgap current of the SIS junction.
    '''
    eV = lambda v: np.multiply(const.e,v)
    eV2kT = lambda v,t : np.divide(eV(v),2*const.k*t)
    temdecay = lambda t,d : np.exp(-np.divide(d,np.multiply(t,const.k))) # constant against voltage
    sqr = lambda v,d: np.sqrt(np.divide(2*d,np.add(eV(v),2*d)))
    sinhK0 = lambda v,t: np.multiply(np.sinh(eV2kT(v,t)),kn(0,eV2kT(v,t)))#approaches expectation value
    return np.multiply(np.multiply(np.multiply(np.multiply(
                                                np.divide(2,
                                                      np.multiply(const.e,rN)),
                                                temdecay(t,d)),
                                                sqr(v,d)),
                                                np.add(eV(v),d)),
                                                sinhK0(v,t))
                                                