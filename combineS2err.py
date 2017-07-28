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

# We map back atomic indices to Cyp residues
# list below is N-terminal residue (was excluded)
# and all Pro residues (no NH)
map_indices = [1,4,16,30,58,95,105]

avgS2 = {}
cumweight = 0.0
for state in weights.keys():
    cumweight += weights[state]
    stream = open("%s/S2tol01-2blocks.dat" % state,"r")
    buffer = stream.readlines()
    stream.close()
    for line in buffer:
        #print (line)
        id, S2, err = line.split()
        S2 = float(S2)
        err = float(err)
        id = int(id)
        offset = 0
        for value in map_indices:
            if (id >= (value-1)):
                offset += 1
        id = id + offset
        try:
            avgS2[id]
        except KeyError:
            avgS2[id] = [0.0, 0.0]
        avgS2[id][0] += weights[state]*S2
        avgS2[id][1] += weights[state]*math.pow(err,2)

#print (avgS2)

S2keys = list(avgS2.keys())
S2keys.sort()

for key in S2keys:
    avg = avgS2[key][0]/cumweight
    err = math.sqrt(avgS2[key][1]/cumweight)
    print ("%8d %8.5f %8.5f" % (key, avg, err))
