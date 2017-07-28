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

parser = argparse.ArgumentParser(description="extracts S2 order parameters from autocorrelation functions.",
                                 epilog="getS2 is distributed under the GPL. For more information please visit ",
                                 prog="getS2", add_help=False)

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-h', '--help', action='store_true', dest='hi')

parser.add_argument('-i', '--input', nargs=1,
                    help="Supply the name of the xvg file containing raw ACFs.")


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
    timestep = 5.0 # in ps
    tinf = 163840 # ps gives 32768 data points
    tmin = 50000.0
    tol = 0.10
    for ACFidx in rawACFS.keys():
        t = 0.0
        c = 0
        S2 = 0.0
        ACF = np.array(rawACFS[ACFidx])
        for x in range(0,len(ACF)):
            t += timestep
            #print t
            if (t > tmin):
                S2 += ACF[x]
                c += 1
            last = ACF[x]
        S2 /= c
        
        if ( abs(S2-last) > tol):
            print ("%s 0.1 #NOT CONVERGED!" % ACFidx)
        else:
            print ("%s %8.5f" % (ACFidx,S2))
        #sys.exit(-1)
