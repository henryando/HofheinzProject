#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 00:31:45 2020

@author: sammarkoff
"""
#generates plot 1c

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from math import *

#intitialize states and operators
psi0 = tensor(fock(2,1), fock(10,0))
a  = tensor(qeye(2), destroy(10))
sm = tensor(destroy(2), qeye(10))

#fixed cavity quibit interaction term
omega = 2 * np.pi * 19e-3

#coherence times (Gdq = dephasing rate)
T1q = 650
T2q = 150
T1r = 3500
Gdq = ((1/T2q)-(1/(2*T1q)))/2 

#collapse operators
Ca = sqrt(1/T1r)*a
Cq1 = sqrt(1/T1q)*sm
Cq2 = sqrt(Gdq)*(2*sm*sm.dag() - 1)

#generate matrix of evolutions over a detuning frequency range of 40 to -40 MHz 
times = np.linspace(0.0, 500, 1000)
DRange = np.linspace(40.0e-3,-40.0e-3,100)
n = DRange.size
result = np.zeros(100)
resultMatrix1 = np.zeros((100,1000))
for i in range(n):
    H = 0.5 * omega * sm.dag() * a + 0.5 * omega * sm * a.dag()+2*pi*DRange[i]*sm.dag()*sm
    result = mesolve(H, psi0, times, [Ca,Cq1,Cq2], [sm.dag()*sm])
    resultMatrix[i] = result.expect[0]

#plot
plt.imshow(resultMatrix,aspect = 3.0, extent=[0,500,-40,40])
plt.xlabel('Time (ns)')
plt.ylabel('Frequency (MHz)')
plt.colorbar(label = 'Excited State Probability')
plt.show()




