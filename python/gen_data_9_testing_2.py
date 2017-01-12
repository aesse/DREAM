# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:10:31 2016

@author: aesseb
"""

"""
Data dependent on TF:
    Motif
Data dependent on Cell:
    Gene expression
    DHS peaks
    DHS general footprint
Data dependent on TF and Cell:
    DHS Motif Footprint (cell too)
    ChIP-seq labels - split into cell type!
Independent data:
    Distance to TSS
"""

import seqdata as sq
from sequence import *
import math
import string
import sys
import pysam

def getHits(bedFile, c,s,e):
    hits = bedFile.fetch(c,s,e)
    c = 0
    for h in hits:
        c += 1
    return c

def getClosest(bedFile, c,s,e):
    cont = True
    snew = s
    enew = e
    try:
        while cont:
            if s-snew > 500000:
                return "Null"
            hits = bedFile.fetch(c,snew,enew)
            if hits is not None:
                minDist = 10000000000
                result = False
                for hit in hits:
                    result = True
                    hit = hit.split("\t")
                    sd,ed = int(hit[1]),int(hit[2])
                    if (s<sd and e>ed) or (s<ed and ed<e) or (s<sd and sd<e) or (sd<s and ed>e):
                        return 0
                    else:
                        dist = min(abs(s-sd), abs(s-ed), abs(e-sd), abs(e-ed))
                    if dist < minDist:
                        minDist = dist
                if result:
                    return minDist
                else:
                    snew -= 1000
                    enew += 1000
    except ValueError:
        return "Null"

chromosomeSizes = {"chr1":249250621,
                    "chr2":243199373,"chr3":198022430,"chr4":191154276,"chr5":180915260,
                    "chr6":171115067,"chr7":159138663,"chrX":155270560,"chr8":146364022,
                    "chr9":141213431,"chr10":135534747,"chr11":135006516,"chr12":133851895,
                    "chr13":115169878,"chr14":107349540,"chr15":102531392,"chr16":90354753,
                    "chr17":81195210,"chr18":78077248,"chr20":63025520,"chrY":59373566,
                    "chr19":59128983,"chr22":51304566,"chr21":48129895}

#tf = "CTCF"
#tfid = "MA0139.1"
#motif_size = 19
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "99" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat+tf
#cell = ["PC-3"]
##cell = "induced_pluripotent_stem_cell"
#locFile = ""
#chrm = ""

#tf = "E2F1"
#tfid = "MA0024.3"
#motif_size = 12
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "11" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat+tf
#cell = "K562"

#tf = "EGR1"
#tfid = "MA0162.2"
#motif_size = 14
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "122" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat+tf
#cell = "liver"

#tf = "FOXA1"
#tfid = "MA0148.3"
#motif_size = 15
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "108" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat+tf
#cell = "liver"

#tf = "Foxa2"
#tfid = "MA0047.2"
#motif_size = 12
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "30" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat+tf
#cell = "liver"

#tf = "Gabpa"
#tfid = "MA0062.2"
#motif_size = 11
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "40" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat+tf
#cell = "liver"

#tf = "Hnf4a"
#tfid = "MA0114.3"
#motif_size = 16
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "84" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat+tf
#cell = "liver"

#tf = "JUND"
#tfid = "MA0491.1"
#motif_size = 11
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "157" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat+tf
#cell = "liver"

#tf = "MAX"
#tfid = "MA0058.3"
#motif_size = 10
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "37" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat+tf
#cell = "liver"

#tf = "REST"
#tfid = "MA0138.2"
#motif_size = 21
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "98" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat+tf
#cell = "liver"

#tf = "NANOG"
#tfid = "NANOG_HUMAN.H10MO.A"
#motif_size = 17
#tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
#tfnum = "522" #Found by looking inside one of the footprint result folders
#tf_lab = tfnum+"-"+tfid_flat
#cell = "induced_pluripotent_stem_cell"

#tf = "TAF1"
#cell = "liver"

tf = sys.argv[1]
tfid = sys.argv[2]
motif_size = int(sys.argv[3])
tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
tfnum = sys.argv[4] #Found by looking inside one of the footprint result folders
tf_lab = tfnum+"-"+tfid_flat+tf
cell = sys.argv[5]
#cell = "induced_pluripotent_stem_cell"
locFile = sys.argv[6]
chrm = sys.argv[7]

#Use fasta file to allow motif searching for each window
#locations = readFastaFile("../Locations/ChIP_label_locs_cons_50_b_ran_set_"+str(motif_size)+".fa")
locations = readFastaFile(locFile)

#Load PWM/motif to search through each 50bp window
if tf == "NANOG" or tf == "ATF2":
    motifPWM = open("Motifs/"+tfid+".pwm")
    tf_lab = tfnum+"-"+tfid_flat
else:
    motifPWM = open("Motifs/"+tfid+"_"+tf+".pwm")
motifPWM = motifPWM.readlines()
motif = [] # Require list of distributions
for m in motifPWM:
    m = m.strip().split("\t")
    d = {}
    d['A'] = m[0]
    d['C'] = m[1]
    d['G'] = m[2]
    d['T'] = m[3]
    d['N'] = m[4]
    motif.append(Distrib(DNA_Alphabet_wN, d))
motif = PWM(motif)
motifRev = motif.getRC()

motifMat = motif.getMatrix()
consensus = ""
lowConsensus = ""
for col in range(0, len(motifMat[0])):
    maxVal = -100
    maxRow = 0
    minVal = 100
    minRow = 0
    for row in range(0, len(motifMat)):
        val = motifMat[row][col]
        if val > maxVal:
            maxVal = val
            maxRow = row
        if val < minVal and val > -20:
            minVal = val
            minRow = row
    consensus += motif.alphabet[maxRow]
    lowConsensus += motif.alphabet[minRow]
pwmScore, pwmIndex, valid = motif.maxscore(consensus)
apwmScore, apwmIndex, avalid = motif.maxscore(lowConsensus)

maxPWMScore = round(pwmScore, 4)
minPWMScore = round(apwmScore, 4)

positives = ["0", "3", "4", "12", "11", "13"]

out = open(tf+"_"+cell+"_testing_data_chrm_"+chrm+".out", 'w')
header = "chrm\tstart\tstop\tTF\tCell\tGCWindow\tMotifScore\tNormScore\tValidity\tdhsPeakLabel\tdhsGFootLabel\tdhsGFootCount\tdhsSFootLabel\tGeneName\tDistance\tlogDistance\tExpressionCluster\tPromoterState\tLocation\tExpression\tlogExpression\tExpression2\tlogExpression2\tOverlaps\n"
out.write(header)

dhsPeakFile = "DNASE_peaks_testing/DNASE."+cell+".conservative_"+chrm+".narrowPeak.gz"
dhsGenFile = "Gen_footprints_testing/"+cell+"_footprints_"+chrm+".bed.gz"
dhsSpecFile= "Spec_footprints_testing/footprints_"+cell+"/beds/"+tf_lab+".top.full.bed.gz"
#To calculate distance
geneLocs = "Gene_loc_testing/labelled_clusters_"+cell+"_"+chrm+".out"

#bedfiles
geneLocs = sq.BedFile(geneLocs, 'Optional')

#tabixFiles
dhsPeaks = pysam.TabixFile(dhsPeakFile)
dhsGprints = pysam.TabixFile(dhsGenFile)
dhsSprints = pysam.TabixFile(dhsSpecFile)

bedFiles = [motif, dhsPeaks, dhsGprints, dhsSprints, geneLocs]
bedLabels = ["Motif", "dhsPeak", "dhsGFoot", "dhsSFoot", "geneLoc"]   

print "Data loaded"     

for locs in locations: #Searching through fasta file
    split = locs.name.split(':')
    chrom = split[0].strip('>')
    s = int(split[1].split('-')[0])+(25) #Each fasta entry is extended by 25bp 
    e = int(split[1].split('-')[1])-(25) #to allow extended motif scanning
    sequence = locs.sequence.tostring()[25-motif_size:75+motif_size]
    
    loc = sq.BedEntry(chrom, s, e)
    labels = []
    labels.append(tf)
    labels.append(cell)
    c,s,e = loc.loc()

    for b in range(len(bedFiles)):
        data = bedLabels[b]
        if data == "dhsPeak" or data == "dhsSFoot" or data == "dhsGFoot":
            closeDistance = getClosest(bedFiles[b], c,s,e)
            labels.append(str(closeDistance))
#                print "Closest = " + str(closeDistance)
            if data == "dhsGFoot":
                try:
                    labels.append(str(getHits(bedFiles[b], c, s-225, e+225)))
                except ValueError:
                    labels.append("Null")
        elif data == "Motif":
            pwmScore, pwmIndex, valid = motif.maxscore(sequence)
            pwmRevScore, pwmRevIndex, rValid = motifRev.maxscore(sequence)
            maxScore = max(pwmScore, pwmRevScore)
            maxScore = round(maxScore, 4)
            valid = valid or rValid
            GCWindow = (float(sequence.count('G')) + float(sequence.count('C'))) / len(sequence)
            if valid:
                normScore = (maxScore-minPWMScore)/(maxPWMScore-minPWMScore)
            else:
                normScore = "Null"
            labels = labels + [str(GCWindow), str(maxScore), str(normScore), str(valid)]
        elif data == "geneLoc":
            closest = bedFiles[b].closest(loc)
            #Use closest because at this stage we don't care about multiple overlaps
            if closest is not None:
                dist, entry, insp = closest
                chrm, start, end = entry.loc()
                strand = entry.name
                #Use strand because we want distance to TSS - depends on gene orientation
                if strand == "+":
                    distance = abs((s+(e-s)) - start) #distance between middle of window and start
                    promoter = sq.BedEntry(chrm, start-2000, start+2000)
                elif strand == "-":
                    distance = abs((s+(e-s)) - end) #distance between middle of window and end
                    promoter = sq.BedEntry(chrm, end-2000, end+2000)
                else:
                    print "Location: " + c+str(s)+str(e) + " had no strand information for ", entry.loc()
                    sys.exit()
                    
                #Count overlaps in window
                nStart = s-500000
                if nStart < 0:
                    nStart = 1
                nEnd = e+500000
                if nEnd > chromosomeSizes[c]:
                    nEnd = chromosomeSizes[c]
                window = sq.BedEntry(c, nStart, nEnd)
                overlaps, entries = bedFiles[b].countOverlap(window)                    
                    
                if sq.overlap(loc.loc(), promoter.loc()) > 0: #Overlaps with promoter
                    labels = labels + [entry.strand, str(distance),str(math.log(distance+1)), entry.thickStart, entry.summit, "Promoter", str(entry.thickEnd), str(math.log(float(entry.thickEnd)+0.0001)), "Null", "Null", str(overlaps)]
                elif dist == 0: #overlaps gene but not promoter
                    labels = labels + [entry.strand, str(distance),str(math.log(distance+1)), entry.thickStart, entry.summit,"Genic", str(entry.thickEnd), str(math.log(float(entry.thickEnd)+0.0001)), "Null", "Null", str(overlaps)]
                else:
                    if overlaps > 0 and len(entries) > 0:
                        maxExpr = -1000.0
                        minExpr = 100000.0
                        for ent in entries:
                            val = float(ent.thickEnd)
                            if val < minExpr:
                                minExpr = val
                            if val > maxExpr:
                                maxExpr = val
                        labels = labels + [entry.strand, str(distance),str(math.log(distance+1)), entry.thickStart, entry.summit,"Distal", str(maxExpr), str(math.log(float(maxExpr)+0.0001)), str(minExpr), str(math.log(float(minExpr)+0.0001)), str(overlaps)]
                    else:
                        labels = labels + [entry.strand, str(distance),str(math.log(distance+1)), entry.thickStart, entry.summit,"Distal", "0", "0", "0", "0", str(overlaps)]
            else:
                print "Location: " +c+str(s)+str(e) + " had no closest entry in " + bedLabels[b]
                labels = labels + ["Null"]*11           
        else:
            print "Invalid data label"
        
    line = "\t".join(labels)
    fin = "%s\t%s\t%s\t%s\n"%(c,s,e,line)
    out.write(fin)
    out.flush()
out.close()




