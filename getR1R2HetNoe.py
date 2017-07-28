#!/usr/bin/python
import sys
import numpy as np
from math import log10

#R1=d00/4(J(wI-wS) + 3J(wI)+6J(wI+wS))
#R2=d00/8(4J(0) + J(wI-wS)+3J(wI)+6J(wS) +6J(wI+wS))
#HetNOE= 1 + d00/4gammaH/gammaN (6J(wI+wS) - J(wI-wS))

# d00= { gammaH gammaN (h/(8pi)) / r_hn^3)^2
# wI is the frequency of Nitrogen
# wS is the frequency of the Proton
# gammaH gyromagnetic ratio of Proton
# gammaN gyromagnetic ratio of Nitrogen
# rNH = bond length (say 0.96 Angstrom)
# h Planck's constant 6.626E-34 J.s-1
# On a 600 MHz spectrometer
# wH = 600 MHz --> log10(600E6*1E-12) = -3.22
# wN = 60 MHz --> log10(600E-6*1E-12) = -4.30

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def interpolate(wvals,freq,Jvals):
    nearest_freq = find_nearest(wvals,freq)
    nearest_freq_minone = nearest_freq - 1
    nearest_freq_plusone = nearest_freq + 1
    if (   abs( wvals[nearest_freq_minone] - wI)
           < abs( wvals[nearest_freq_plusone] - wI) ):
        next_nearest_freq = nearest_freq_minone
    else:
        next_nearest_freq = nearest_freq_plusone
    # print (freq,wvals[nearest_freq],wvals[next_nearest_freq])
    if wvals[nearest_freq] < wvals[next_nearest_freq]:
        lower = nearest_freq
        upper = next_nearest_freq
    else:
        lower = next_nearest_freq
        upper = nearest_freq

    # print (lower,upper)
    d = (freq-wvals[lower])/(wvals[upper]-wvals[lower])
    # print (d)
    J_freq = d*Jvals[upper]+(1-d)*Jvals[lower]
    # print (J_freq)
    return J_freq

freq = 600E6 # 600 MHz
wS = log10(freq*1E-12) # -3.22 at 600 MHz, input in ps-1
wI = log10(freq/10.0*1E-12) # -4.30 at 600 Mhz
wISplus = log10((freq+freq/10.0)*1E-12)
wISminus = log10((freq-freq/10.0)*1E-12)

try:
    sd = sys.argv[1:]
except IndexError:
    print ("USAGE is script spectral-density.dat")
    sys.exit(-1)

sd.sort()

results = {}

for sdfile in sd:
    base,idx = sdfile.split("-")
    idx = idx.rstrip(".dat")
    idx = int(idx)
    stream = open(sdfile,'r')
    buffer = stream.readlines()
    stream.close()

    warray = []
    Jarray = []
    for line in buffer:
        w, J = line.split()
        w = float(w)
        J = float(J)
        warray.append(w)
        Jarray.append(J)

    wvals = np.array(warray)
    Jvals = np.array(Jarray)

    J_0 = Jvals[0]
    J_wI = interpolate(wvals,wI,Jvals)
    J_wS = interpolate(wvals,wS,Jvals)
    J_wISplus = interpolate(wvals,wISplus,Jvals)
    J_wISminus = interpolate(wvals,wISminus,Jvals)
    # print (0, wI, wISminus, wS, wISplus)
    # print (J_0, J_wI, J_wISminus, J_wS, J_wISplus)

    d00 = 1.0# Pain to work out constant...
    gammaH = 42.577#  MHz.T-1
    gammaN = -4.316#  MHz.T-1

    R1 = (d00/4.0) * ( J_wISminus + 3*J_wI + 6*J_wISplus )
    R2 = (d00/8.0) * ( 4 * J_0 + J_wISminus + 3*J_wI + 6*J_wS + 6*J_wISplus )
    HetNOE = 1 + (d00/4.)*(6*J_wISplus-J_wISminus)*(1/R1)*(gammaH/gammaN)
    #print (idx,R1,R2,HetNOE)
    results[idx] = (R1,R2,HetNOE)

keys = list(results.keys())
keys.sort()

for k in keys:
    print (k,results[k][0],results[k][1],results[k][2])
