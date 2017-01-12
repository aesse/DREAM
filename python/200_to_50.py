# -*- coding: utf-8 -*-
"""
Created on Mon Aug 08 13:14:51 2016

@author: aesseb
"""

chipData = "ChIP_labels_mixed_s.bed"
chipData = open(chipData)
chipData = chipData.readlines()

out = open("ChIP_labels_mixed_50.bed", 'w')
disp = False
for i in range(len(chipData)-1):
    if i%100000 == 0:
        print i
    chip = chipData[i].strip().split("\t")
    c,s,e = chip
    chipNext = chipData[i+1].strip().split("\t")
    cn,sn,en = chipNext

    #Does the next entry overlap with the current?
    #Then write 50bps based on current start and move on
    if int(sn) < int(e):
        out.write("%s\t%s\t%d\n"%(c,s,int(s)+50))
        if disp:
            print "write simple"
    #if the next entry doesn't overlap, write full 200bps into 50bps
    elif int(sn) >= int(e):
        if disp:
            print "write extended"
        out.write("%s\t%s\t%d\n"%(c,s,int(s)+50))
        for j in range(3):
            s = int(s) + 50
            out.write("%s\t%d\t%d\n"%(c,s,s+50))
    #if we're transitioning between chromosomes, write full 200bps into 50bps
    #In this case, sn would be < e but with mismatched chromosome label
    elif not c == cn:
        if disp:
            print "write next chrmsm"
        out.write("%s\t%s\t%d\n"%(c,s,int(s)+50))
        for j in range(3):
            s = int(s) + 50
            out.write("%s\t%d\t%d\n"%(c,s,s+50))   
out.close()
    
    