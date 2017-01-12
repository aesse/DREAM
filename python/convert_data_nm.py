# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 22:12:13 2016

@author: aesseb
"""

import math
import sys
import random

#base = "EGR1_sample_testing_data_7"
base = sys.argv[1]

data = open(base+".out")

valid = []
data = data.readlines()
header = ""
for d in data:
    d = d.strip().split("\t")
    if d[0].startswith("chrm"):
        header = "\t".join(d)
        continue    
    for pos in [6,7]:
        try:
            peak = int(d[pos])
            if peak == 0:
                d[pos] = "0"
            elif peak <= 50:
                d[pos] = "1"
            elif peak <= 200:
                d[pos] = "2"
            elif peak <= 500:
                d[pos] = "3"
            elif peak <= 1000:
                d[pos] = "4"
            elif peak <= 20000:
                d[pos] = "5"
            elif peak <= 100000:
                d[pos] = "6"
            else:
                d[pos] = "7"
        except ValueError:
            d[pos] = "Null"
    d[8] = str(math.log(float(d[8])+1))
    d[19] = str(math.log(float(d[19])+1))
    valid.append("\t".join(d))


print len(valid)

splitV = int(math.floor(len(valid)/10.0))

points = random.sample(xrange(0, len(valid)), len(valid))
for v in range(0, 10):
    start = v * splitV
    end = (v+1) * splitV
    out = open(base+"_convert_"+str(v)+".out", 'w')
    out.write(header+"\n")
    for val in points[start:end]:
        out.write(valid[val]+"\n")
    out.close()

    