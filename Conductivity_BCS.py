#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:51:30 2020

@author: wenninger
"""
import numpy as np
import scipy.constants as const
from scipy.integrate import quad
from scipy.special import expit

def conductivity_two_fluid(freq,sigmaN,tN,london):
    '''This function computes the complex conductivity according to the two-fluid model.
    The equations are from Kautz 1978.
    
    inputs
    ------
    freq: float
        The frequency of the signal in/on the material.
    sigmaN: float
        The normal state conductivity at the critical temperature.
    tN: float
        The temperature normalised to the critical temperature.
    london: float
        The london penetration depth of the material at 0 K.
        
    returns
    -------
    complex float
        The complex conductivity.
    '''
    sigma1 = np.multiply(sigmaN,tN)
    sigma2 = np.divide(np.subtract(1,np.power(tN,4)),np.multiply(2*np.pi*freq*const.mu_0,np.square(london)))
    return np.add(sigma1 , 1j*sigma2)

def conductivity_BCS(freq,sigmaN=1,te=4,delta= 1.5e-3*const.e):
    '''This function computes the complex conductivity according to the BCS weak-coupling theory.
    The equations are from Mattis and Bardeen 1958.
    Evaluated with Kooi's PhD thesis.
    
    inputs
    ------
    freq: float
        The frequency of the signal in/on the material.
    sigmaN: float
        The normal state conductivity at the critical temperature.
    te: float
        The physical temperature of the material.
    delta: float
        The Cooper pair binding energy at temperature T.
        
    returns
    -------
    complex float
        The complex conductivity.
    '''
    hbarOmega = np.multiply(const.hbar, np.multiply(2*np.pi,freq))
    #expit(x) = 1/(1+exp(-x))
    f = lambda epsilon: expit(-np.divide(epsilon,np.multiply(const.k,te)))          
    
    epsilon1 = lambda epsilon: np.sqrt(np.subtract(np.square(epsilon),np.square(delta)))
    epsilon2 = lambda epsilon: np.sqrt(np.subtract(np.square(np.add(epsilon,hbarOmega)),np.square(delta)))
    g = lambda epsilon: abs(np.divide(np.add(np.square(epsilon),np.add(np.square(delta),np.multiply(hbarOmega,epsilon))),
                                  np.multiply(epsilon1(epsilon),epsilon2(epsilon))))
    sigma1integral1 = lambda epsilon: np.multiply(np.divide(2,hbarOmega),
                                                  np.multiply(np.subtract(f(epsilon),f(np.add(epsilon,hbarOmega))),g(epsilon))
                                                  )
    sigma1integral2 = lambda epsilon: np.multiply(np.divide(1,hbarOmega),
                                                  np.multiply(np.subtract(1,2*f(np.add(epsilon,hbarOmega))),g(epsilon))
                                                  )
    sigma2integral = lambda epsilon: np.multiply(np.divide(1,hbarOmega),
                                                 np.divide(np.multiply(np.subtract(1,2*f(np.add(epsilon,hbarOmega))),
                                                                       np.add(np.square(epsilon),np.add(np.square(delta),np.multiply(hbarOmega,epsilon)))
                                                      ),
                                                  np.multiply(np.sqrt(np.subtract(np.square(delta),np.square(epsilon))),epsilon2(epsilon))))
                                                 
    sigma1 = quad(sigma1integral1,delta+1e-40,delta*15)[0]
    
    if hbarOmega>2*delta:
        sigma1 += quad(sigma1integral2,np.subtract(delta+1e-40,hbarOmega),-delta-1e-40)[0]
        sigma2 = quad(sigma2integral,-delta+1e-40,delta-1e-40)[0]
    else:
        sigma2 = quad(sigma2integral,delta-hbarOmega+1e-40,delta-1e-40)[0]
        
    return np.multiply(sigmaN,np.add(sigma1,1j*sigma2))
    
     
     
     