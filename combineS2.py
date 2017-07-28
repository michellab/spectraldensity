#!/usr/bin/python
import os,sys,math

stream = open('weights.dat','r')
buffer = stream.readlines()
stream.close()
weights = {}
for line in buffer:
    elems = line.split()
    weights[elems[0]] = float(elems[1])

#print (weights)

avgS2 = {}
cumweight = 0.0
for state in weights.keys():
    cumweight += weights[state]
    stream = open("%s/S2tol01.dat" % state,"r")
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
            avgS2[id] = 0.0
        avgS2[id] += weights[state]*S2

#print (avgS2)

S2keys = list(avgS2.keys())
S2keys.sort()

for key in S2keys:
    avgS2[key] = avgS2[key]/cumweight
    print ("%8d %8.5f " % (key, avgS2[key]))
