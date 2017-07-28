#!/usr/bin/python
import os,sys,math
import numpy as np

S2files = sys.argv[1:]

#print (weights)

avgS2 = {}
cumweight = 0.0
for S2file in S2files:
    cumweight += 1.0
    stream = open("%s" % S2file,"r")
    buffer = stream.readlines()
    stream.close()
    for line in buffer:
        #print (line)
        if line.find("NOT CONVERGED") > 0:
            id, S2, trash1, trash2 = line.split()
        else:
            id, S2 = line.split()
        S2 = float(S2)
        id = int(id)
        try:
            avgS2[id]
        except KeyError:
            avgS2[id] = []
        avgS2[id].append(S2)

#print (avgS2)

S2keys = list(avgS2.keys())
S2keys.sort()

for key in S2keys:
    array = np.array(avgS2[key])
    mean = array.mean()
    sterr = array.std()/math.sqrt(len(array))
    print ("%8d %8.5f %8.5f" % (key, mean, sterr))
