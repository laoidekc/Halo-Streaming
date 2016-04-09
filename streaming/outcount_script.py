#!/usr/bin/env python
from numpy import *

data = loadtxt("outcount.dat")

def column(matrix, i):
    return [row[i] for row in matrix]

output = column(data,1)

for i in range(len(output)-1):
	data[i,1] = data[i+1,1] - data[i,1]
	data[i,0] = data[i+1,0]

data = delete(data,len(output)-1,axis=0)

savetxt('outcount_diffs.dat',data)
