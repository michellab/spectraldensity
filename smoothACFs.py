#!/usr/bin/python
#@author Julien Michel March 2017

import numpy as np
import sys,os
import argparse
import scipy
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit
from scipy.fftpack import fft
from math import log,log10

parser = argparse.ArgumentParser(description="Smooth MD derived ACFs for spectral density predictions.",
                                 epilog="smoothACFS is analyse_freenrg is distributed under the GPL. For more information please visit ",
                                 prog="smoothACFS", add_help=False)

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-h', '--help', action='store_true', dest='hi')

parser.add_argument('-i', '--input', nargs=1,
                    help="Supply the name of the xvg file containing raw ACFs.")


# Fitting A exp(-b x)
def func(x, a, b):
    return a*np.exp(-b*x)

def loadACFS(inputfile):
    stream = open(inputfile,'r')
    buffer = stream.readlines()
    stream.close()
    idx = 0
    acfs = {}
    acfs[0] = []
    for line in buffer:
        if ( line.startswith("#") or
             line.startswith("@") ):
            continue
        if line.startswith("&"):
            idx += 1
            acfs[idx] = []
            continue
        elems = line.split()
        val = float(elems[1])
        acfs[idx].append(val)
    del acfs[idx]

    return acfs

if __name__ == '__main__':
    args = parser.parse_args()
    if args.input:
        input_acfs = args.input[0]
    else:
        input_acfs = None
    if input_acfs is None:
        print ("Please supply an input xvg file!")
        sys.exit(-1)
    rawACFS = loadACFS(input_acfs)

    #For each ACF
    region2max = 0.3
    region2min = 0.05
    timestep = 5.0 # in ps
    tinf = 163840 # ps gives 32768 data points
    for ACFidx in rawACFS.keys():
        # 1) Find interval where C(t) value decays from 0.3 to 0.05
        ACF = np.array(rawACFS[ACFidx])
        xmax = 0
        xmin = 99999
        for x in range(0,len(ACF)):
            if ACF[x] < region2max:
                xmax = x
                break
        for x in range(0,len(ACF)):
            if ACF[x] < region2min:
                xmin = x
                break
        #print (xmax,xmin)
        region2y = ACF[xmax:xmin]
        #print (region2y)
        region2x = []
        for x in range(0,len(region2y)):
            region2x.append(x*timestep)
        #print (region2x)
        # 2) Fit region2 to single exponential
        popt, pcov = curve_fit(func, region2x, region2y, p0=(1, 1e-6))
        # FIXME: Estimate uncertainties from pcov
        # 3) Assemble smoothed ACF
        # Part A) if t < tmax: raw ACF
        # Part B) elif t < tmin: weighted average of ACF + fit over [tmax,tmin]
        # Part C) elif t < tinf: fit
        tmax = xmax*timestep
        tmin = xmin*timestep
        smoothACF = []
        x = -1
        for t in range(0,tinf,int(timestep)):
            x += 1
            if (t < tmax):
                smoothACF.append(ACF[x])
            elif (t < tmin):
                weight = (t-tmax)/(tmin-tmax)
                val = weight * func(t-tmax,popt[0],popt[1]) + (1-weight)*ACF[x]
                smoothACF.append(val)
            elif (t < tinf):
                val = func(t-tmax,popt[0],popt[1])
                smoothACF.append(val)
        #print (smoothACF)
        stream = open('smoothACF-%s.dat' % ACFidx,'w')
        for x in range(0,len(smoothACF)):
            t = x*timestep
            val = smoothACF[x]
            stream.write("%8.2f %e\n" % (t,val))
        stream.close()
        # FFT
        # see https://docs.scipy.org/doc/scipy-0.18.1/reference/tutorial/fftpack.html
        J = fft(np.array(smoothACF))
        # Scale result
        # But Robustelli et al. have http://pubs.acs.org/doi/full/10.1021/ct400654r
        # J(w) = (2/5) Re [ \int C(t) exp (-iwt) dt ]
        scale = 2/5.
        Jscl = []
        for x in range(0,len(J)):
            Jx = J[x]
            Jsclx = scale*Jx
            Jscl.append(Jsclx)
        Jscl = np.array(Jscl)
        stream = open('spectraldensity-%s.dat' % ACFidx,'w')
        N = len(Jscl)
        T = 1/timestep
        ws = np.linspace(0.0,1.0/(2.0*T),N/2)
        for x in range(1,int(len(Jscl)/2)):
            # must check what w is exactly
            w = log10(ws[x])
            val = Jscl[x]
            #print (w,val)
            stream.write("%e %e\n" % (w,val))
        stream.close()
        #sys.exit(-1)
