# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 09:24:45 2016

@author: aesseb
"""

import seqdata as sq
import pysam

#windows = sq.BedFile("../DNASE/peaks/relaxed/DNASE.HepG2.relaxed.narrowPeak", 'Optional')
windows = sq.BedFile("../ChIP-seq/peaks/conservative/Unions/ChIPseq.full.union.bed", 'Optional')

label_list = ["../ChIP-seq/labels/SRF.train.labels.H1-hESC.bed.gz"]

pysam_files = []
for lab in label_list:
    pysam_files.append(pysam.TabixFile(lab))

print "Files loaded"  

valEr = 0
for py in range(0, len(pysam_files)):
    pyFile = pysam_files[py]
    name = label_list[py]
    out = open("ChIP_labels_restricted.out", 'w')
    for window in windows:
        c,s,e = window.loc()
        try:
            hits = pyFile.fetch(c, s, e)
            if hits is not None:
                for h in hits:
                    h = h.split("\t")
                    out.write("%s\t%s\t%s\t%s\n"%(h[0], h[1], h[2], h[3]))
            out.flush()
        except ValueError:
            valEr += 1
            continue
    out.close