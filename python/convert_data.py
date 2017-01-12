# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 13:56:03 2016

@author: aesseb
"""
import math
import sys
import random

#base = "EGR1_sample_testing_data_7"
base = sys.argv[1]

data = open(base+".out")

valid = []
invalid = []
data = data.readlines()
header = ""
for d in data:
    d = d.strip().split("\t")
    if d[0].startswith("chrm"):
        header = "\t".join(d)
        continue    
    for pos in [9,10,12]:
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
    d[11] = str(math.log(float(d[11])+1))
    d[23] = str(math.log(float(d[23])+1))
    if d[8] == "True":
        valid.append("\t".join(d))
    elif d[8] == "False":
        d[6] = "Null"
        d[7] = "Null"
        invalid.append("\t".join(d))

print len(valid)
print len(invalid)

validCount = 5
invalidCount = 5

if len(invalid) < 400000:
    if len(invalid)%100000 > 2:
        invalidCount = 2
        validCount = 8
    else:
        invalidCount = 1
        validCount = 9
elif len(invalid) == 0:
    invalidCount = 0
    validCount = 10

if len(valid) < 400000:
    if len(valid)%100000 > 2:
        invalidCount = 8
        validCount = 2
    else:
        invalidCount = 9
        validCount = 1
        
splitI = int(math.floor(len(invalid)/float(invalidCount)))
splitV = int(math.floor(len(valid)/float(validCount)))

points = random.sample(xrange(0, len(valid)), len(valid))
for v in range(0, validCount):
    start = v * splitV
    end = (v+1) * splitV
    out = open(base+"_valid_convert_"+str(v)+".out", 'w')
    out.write(header+"\n")
    for val in points[start:end]:
        out.write(valid[val]+"\n")
    out.close()
    
points = random.sample(xrange(0, len(invalid)), len(invalid))
for iv in range(0, invalidCount):
    start = iv * splitI
    end = (iv+1) * splitI
    out = open(base+"_invalid_convert_"+str(iv)+".out", 'w')
    out.write(header+"\n")
    for ival in points[start:end]:
        out.write(invalid[ival]+"\n")
    out.close()

    