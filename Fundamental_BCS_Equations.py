#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:29:36 2019

@author: wenninger
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy.special import expit

def nS_over_nN0(x =np.arange(-.001,.001,0.00001)*const.e,te=4,teC=9.2):
    '''Returns the number of Cooper pairs over the number of electrons in normal state at 0 K, dependent on the energy.
    \frac{N_\text{S}(T)}{N_\text{N}(T=0\,\text{K})} Tinkham equ. 3.73
    
    Default values are used for testing
    
    inputs
    ------
    x: 1D array
        The energy relative to the fermi surface in J.
    te: float or int
        The actual temperature of the junction.
    teC: float or int
        The critical temperture of the superconductor.
        
    returns
    -------
    1d array: 
        The nS_over_nN0 value for each energy level relative to the fermi surface.
    '''
    x=np.array(x)
    dCP = cooperPair_Binding_Energy_over_T(te/teC,cooperPair_Binding_Energy_0K(teC))[0]
    ret = np.zeros_like(x)
    cond=(x>dCP) | (x<-dCP) # the condition of values not zero
    ret[cond] = np.abs(x[cond])/np.sqrt(np.subtract(np.multiply(x[cond],x[cond]),dCP*dCP))
    return ret

def _nS_over_nN0_Test(energyRange=np.arange(-.01*const.e,.01*const.e,.0001*const.e),te=4,teC=10):
    '''Test function of :fun: nS_over_nN0.
    ------
    inputs
    ------
    energyRange: 1D array
        energy values relative to the fermi surface in J tested
    te: float or int
        The actual temperature of the junction
    teC: float or int
        The criticl temperture of the superconductor
    '''
    plt.plot(energyRange,nS_over_nN0(energyRange,te,teC))
    plt.axvline(cooperPair_Binding_Energy_over_T(te/teC,cooperPair_Binding_Energy_0K(teC)))
    plt.axvline(cooperPair_Binding_Energy_over_T(te/teC,cooperPair_Binding_Energy_0K(teC)))
    plt.xlabel('Energy [J]')
    plt.ylabel('Cooper Pair Density [m$^{-3}$]')
    plt.show()

def normalised_CooperPair_Binding_Energy_Tc_Approximation(t_over_Tc):
    '''This function computes the approximation value for the Cooper pair binding energy close to the critical temperature Tc.
    Note that the input and output values are normalised to the critical temperature and the Cooper pair binding energy at 0 K.
    
    inputs
    ------
    t_over_Tc: float or array
        The actual temperature normalised to the critical temperature.
        T/T_c
    
    returns
    float or array
        The cooper pair binding energy normalised on the critical temperature at 0K.
        Delta/Delta_0
    '''
    return 1.74*np.sqrt(1 - t_over_Tc)

def cooperPair_Binding_Energy_0K(tC, bCSfactor=3.52/2.):
    '''This function computes teh Cooper pair binding energy at 0 K using the critical temperature and the BCS factor.

    inputs
    ------
    tC: float or array
        The criticl temperture of the superconductor
    bCSfactor: float
        The BCS factor which is usually 3.5/2.
        Note that the value differs for Niobium.
        
    returns
    -------
    float or array: 
        The Cooper Pair binding energy at 0 K for given Tc.
    '''
    return np.multiply(bCSfactor*const.k,tC)

def normalised_CooperPair_Binding_Energy_over_T(t_over_Tc =np.arange(0,1.001,.001),delta0 = cooperPair_Binding_Energy_0K(9.2), tDebye =  276):
    '''This function computes the normalised Cooper pair binding energy over temperature as given by the  implicite BCS equation (Tinkham equation 3.53)
    This function follows https://physics.stackexchange.com/questions/54200/superconducting-gap-temperature-dependence-how-to-calculate-this-integral

    inputs
    ------
    t_over_Tc: 1d np.array or float
        The temperatures, normalised on the critical temperature, which are evaluated.
    delta0: float
        The Cooper Pair binding energy at 0 K
    tDebye: float
        The Debye temperature of the electron gas in the superconductor.
    
    returns
    -------
    1d array
        The Cooper pair binding energy.
    '''
    omegaDebye = const.k/const.hbar*tDebye
    #Equation to integrate
    funcToIntegrate = lambda z , sig , tau: np.divide(np.tanh(np.multiply(.882 , np.multiply(np.divide(sig , tau),np.sqrt(1 + np.multiply(z , z))))),np.sqrt(1 + np.multiply(z , z)))
    #Equation to Solve
    funcToSolve = lambda sig , *tau: + quad(funcToIntegrate,0,np.divide(const.hbar*omegaDebye,np.multiply(delta0,sig)),args=(sig,tau))[0] -np.arcsinh(const.hbar*omegaDebye/delta0) 
    if isinstance(t_over_Tc*1.0,float):
        return fsolve(funcToSolve,.7, args=(t_over_Tc))
    else:
        deltaDelta0 = []
        for i in t_over_Tc: # probe for T/TC values
            deltaDelta0.append(fsolve(funcToSolve,.7, args=(i)))
        return np.hstack(deltaDelta0)
    #return np.vstack([t_over_Tc,deltaDelta0])

def cooperPair_Binding_Energy_over_T(t_over_Tc =np.arange(0,1.001,.001),delta0 = cooperPair_Binding_Energy_0K(9.2), tDebye =  276):
    '''This function computes the Cooper pair binding energy over temperature as given by the  implicite BCS equation (Tinkham equation 3.53)
    This function follows https://physics.stackexchange.com/questions/54200/superconducting-gap-temperature-dependence-how-to-calculate-this-integral

    inputs
    ------
    t_over_Tc: 1d np.array or float
        The temperatures, normalised on the critical temperature, which are evaluated.
    delta0: float
        The Cooper Pair binding energy at 0 K
    tDebye: float
        The Debye temperature of the electron gas in the superconductor.
    
    returns
    -------
    1d array 
        The Cooper pair binding energy.
    '''
    return delta0*normalised_CooperPair_Binding_Energy_over_T(t_over_Tc,delta0,tDebye)

def gap_Voltage_over_T(t_over_Tc =np.arange(0,1.001,.001),delta0 = 1.8*const.k*9.2, tDebye =  276):
    '''This function computes the gap voltage over temperature as given by the  implicite BCS equation (Tinkham equation 3.53)
    This function follows https://physics.stackexchange.com/questions/54200/superconducting-gap-temperature-dependence-how-to-calculate-this-integral

    inputs
    ------
    t_over_Tc: 1d np.array
        The temperatures, normalised on the critical temperature, which are evaluated.
    delta0: float
        The Cooper Pair binding energy at 0 K
    tDebye: float
        The Debye temperature of the electron gas in the superconductor.
    
    returns
    -------
    1d array
        The Gap Voltage associated to the given temperatures
    '''
    return cooperPair_Binding_Energy_over_T(t_over_Tc,delta0,tDebye)/const.e

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
    
     
     
     
     
     
     
     
     
     
