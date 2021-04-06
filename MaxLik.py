# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 13:57:11 2020

@author: Jeremy Raskop

Thesis Olivier Morrin: https://tel.archives-ouvertes.fr/tel-01066655v2
Papier Lvovsky: https://iopscience.iop.org/article/10.1088/1464-4266/6/6/014/meta

Maybe use qutip?
"""

sigma_0 = 1

import numpy as np
import scipy.special as sp
import math
import matplotlib.pyplot as plt


def hermite(n,x):
    return sp.eval_hermite(n,x)

def wavefunction_fock(n,x,theta,sigma_0 = 1):
    # Normalisation constant:
    N = pow(np.sqrt(2*np.pi) * sigma_0 * 2**n * math.factorial(n), -0.5 )
    # Scaling of x:
    x = x/sigma_0/np.sqrt(2)  
    # Wavefunction:
    return np.exp(1j*n*theta) * N * hermite(n,x) * np.exp(-x**2/2)
    
def projection_operator(m,n,theta,x):
    return wavefunction_fock(m,x,theta)*np.conj(wavefunction_fock(n,x,theta))




x = np.linspace(-6,6,100)
plt.plot(x,wavefunction_fock(2,x,0)*np.conj(wavefunction_fock(2,x,0)))






