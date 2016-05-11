#!/usr/bin/env python
from numpy import *

data = loadtxt("outcountdata.txt")

def column(matrix, i):
    return [row[i] for row in matrix]

output = column(data,1)
output = delete(output,len(output)-1,axis=0)
output = delete(output,len(output)-1,axis=0)
data = delete(data,0,axis=0)

for i in range(len(output)):
	data[i,1] = data[i+1,1] - data[i,1]
	data[i,0] = data[i+1,0]
	if data[i,1]<0.0036 and i>500 and i<10000:
		print i
data = delete(data,len(output),axis=0)

savetxt('outcount_diffs.dat',data)
