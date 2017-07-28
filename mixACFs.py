#!/usr/bin/python
#@author Julien Michel July 2017

import numpy as np
import sys,os
import argparse
import scipy
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit
from scipy.fftpack import fft
from math import log,log10

parser = argparse.ArgumentParser(description="Mix an internal ACF with a global ACF for J(w) predictions.",
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
    TC = 8200.0 # in ps Exp tumbling time found by Arun in the literature
    #TC = 5000.0 # in ps
    region2max = 0.3
    region2min = 0.05
    timestep = 5.0 # in ps
    xinf = 65536 # gives 327680 ps
    for ACFidx in rawACFS.keys():
        ACF = np.array(rawACFS[ACFidx])
        print (ACF)
        observedACF = np.array(np.zeros(xinf))
        #for x in range(0,len(ACF)):
        for x in range(0,xinf):
            time = timestep*x
            globalACF = func(time,1.0,1/TC)
            if (time < 100000):
                observedACF[x] = ACF[x]*globalACF
            else:
                observedACF[x] = ACF[-1]*globalACF
        print (observedACF)
        stream = open('observedACF-%s.dat' % ACFidx,'w')
        for x in range(0,len(observedACF)):
            t = x*timestep
            val = observedACF[x]
            stream.write("%8.2f %e\n" % (t,val))
        stream.close()
        # FFT
        # see https://docs.scipy.org/doc/scipy-0.18.1/reference/tutorial/fftpack.html
        J = fft(np.array(observedACF))
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
        ws = np.linspace(0.0,0.1/(2.0*T),N/2)
        stream.write("-6 %e\n" % (Jscl[0]))
        for x in range(1,int(len(Jscl)/2)):
            # must check what w is exactly
            w = log10(ws[x])
            val = Jscl[x]
            #print (w,val)
            stream.write("%e %e\n" % (w,val))
         #import pdb ; pdb.set_trace()
        stream.close()
         #sys.exit(-1)
