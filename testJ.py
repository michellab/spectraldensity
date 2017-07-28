import numpy as np
import sys,os
import argparse
import scipy
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit
from scipy.fftpack import fft
from math import log,log10

def func(x, a, b):
    return a*np.exp(-b*x)

def func2(x, a1, k1, a2, k2):
    return a1*np.exp(-k1*x)+a2*np.exp(-k2*x)

xinf = 32768
timestep = 5.0
A1 = 1.0
TC1 = 1000.0
A2 = 100.0
TC2 = 100.0
ACF = np.array(np.zeros(xinf))
for x in range(0,xinf):
    t = x*timestep
    #ACF[x] = func(x,A1,1/TC1)
    ACF[x] = func2(x,A1,1/TC1,A2,1/TC2)
J = fft(ACF)
print (J)
scale = 1#2/5.
Jscl = []
for x in range(0,len(J)):
    Jx = J[x]
    Jsclx = scale*Jx
    Jscl.append(Jsclx)
Jscl = np.array(Jscl)
stream = open('test.dat','w')
N = len(Jscl)
T = 1/timestep
ws = np.linspace(0.0,1.0/(2.0*T),N/2)
stream.write("-5 %e\n" % (Jscl[0]))
for x in range(1,int(len(Jscl)/2)):
    # must check what w is exactly
    w = log10(ws[x])
    val = Jscl[x]
    #print (w,val)
    stream.write("%e %e\n" % (w,val))
stream.close()
