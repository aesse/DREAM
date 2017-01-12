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

valid = open(base+"_valid.out", 'w')
invalid = open(base+"_invalid.out", 'w')

data = data.readlines()
header = ""
for d in data:
    d = d.strip().split("\t")
    if d[0].startswith("chrm"):
        header = "\t".join(d)
        valid.write(header)
        invalid.write(header)
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
        valid.write("\t".join(d)+"\n")
    elif d[8] == "False":
        d[6] = "Null"
        d[7] = "Null"
        invalid.write("\t".join(d)+"\n")

print len(valid)
print len(invalid)

valid.close()
invalid.close()