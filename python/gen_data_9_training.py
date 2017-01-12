# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 09:25:02 2016

@author: aesseb
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 10:28:48 2016

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

#tfDetails = {"Arid3a":["MA0151.1", 6,"111", ["HepG2"]],
#            "Atf3":["MA0605.1",8, "209",["H1-hESC","HepG2","K562"]],
#            "ATF7":["MA0834.1", 14, "438", ["K562","GM12878","HepG2"]],
#            "CEBPB":["MA0466.2",10,"133",["HepG2","K562","H1-hESC",""HeLa-S3"]],
#            "CREB1":["MA0018.2",8,"9",["GM12878","HepG2","K562","H1-hESC"]],
#            "CTCF":["MA0139.1", 19, "99",  ["H1-hESC", "HeLa-S3", "HepG2", "K562"]],
#            "E2F1":["MA0024.3",12,"11",["GM12878", "HeLa-S3"]],
#            "E2F6":["MA0471.1",11,"138",["H1-hESC","HeLa-S3"]],
#            "EGR1":["MA0162.2",14,"122",["GM12878","H1-hESC"]],
#            "FOXA1":["MA0148.3", 15, "108",["HepG2"]], 
#            "Foxa2":["MA0047.2", 12, "30", ["HepG2"]], 
#            "Gabpa":["MA0062.2",11,"40",["HepG2","H1-hESC","HeLa-S3","GM12878","MCF-7"]],
#            "GATA3":["MA0037.2",8,"22",["A549"]],
#            "Hnf4a":["MA0114.3",16, "84", ["HepG2"]], 
#            "JUND":["MA0491.1",11,"157",["HepG2","HCT116","HeLa-S3","MCF-7","K562"]],
#            "MAFK":["MA0496.1",15,"162",["HepG2","H1-hESC", "HeLa-S3","GM12878","IMR-90"]],
#            "MAX":["MA0058.3",10,"37",["HepG2","A549","HCT116","K562","H1-hESC","HeLa-S3","GM12878"]],
#            "Myc":["MA0147.2",10,"107",["A549","K562","MCF-7","HeLa-S3"]],
#            "NANOG":["NANOG_HUMAN.H10MO.A", 17,"522", ["H1-hESC"]], 
#            "REST":["MA0138.2",21,"98",["HepG2","H1-hESC","HeLa-S3","MCF-7","Panc1"]],
#            "RFX5":["MA0510.2",16,"176",["HeLa-S3","GM12878","MCF-7"]],
#            "SPI1":["MA0080.4", 14,"57",["GM12878"]],
#            "SRF":["MA0083.3",16,"59",["HepG2","HCT116","K562", "H1-hESC", "GM12878"]],
#            "STAT3":["MA0144.2", 11, "104",["HeLa-S3"]],
#            "Tcf12":["MA0521.1",11,"187",["MCF-7", "H1-hESC", "GM12878"]],
#            "TCF7L2":["MA0523.1",14,"189",["HCT116", "HeLa-S3", "Panc1"]],
#            "TEAD4":["MA0809.1",10,"413",["HepG2", "A549", "HCT116", "K562", "H1-hESC"]],
#            "YY1":["MA0095.2",12,"68",["HepG2", "HCT116", "H1-hESC", "GM12878"]],
#            "ZNF143":["MA0088.2",16,"62",["HepG2", "GM12878","H1-hESC", "HeLa-S3"]]}

#tfDetails = {"Arid3a":["MA0151.1", 6,"111", ["HepG2"]],
#            "FOXA1":["MA0148.3", 15, "108",["HepG2"]], 
#            "CTCF":["MA0139.1", 19, "99",  ["H1-hESC", "HeLa-S3", "HepG2", "K562"]],
#            "E2F6":["MA0471.1",11,"138",["H1-hESC","HeLa-S3"]],
#            "EGR1":["MA0162.2",14,"122",["GM12878","H1-hESC"]],
#            "ZNF143":["MA0088.2",16,"62",["HepG2", "GM12878","H1-hESC", "HeLa-S3"]]}

#Use fasta file to allow motif searching for each window
#locations = readFastaFile("../Locations/ChIP_label_locs_cons_50_b_ran_set_"+str(motif_size)+".fa")

tf = sys.argv[1]
tfid = sys.argv[2]
motif_size = int(sys.argv[3])
tfid_flat = tfid.translate(string.maketrans("",""), string.punctuation)
tfnum = sys.argv[4] #Found by looking inside one of the footprint result folders
tf_lab = tfnum+"-"+tfid_flat+tf
cells = sys.argv[5].split(",") #comma separated list
locFile = sys.argv[6]

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

locations = readFastaFile(locFile)

dataSize = len(locations)/len(cells)

out = open(tf+"_training_data.out", 'w')
header = "chrm\tstart\tstop\tTF\tCell\tGCWindow\tMotifScore\tNormScore\tValidity\tdhsPeakLabel\tdhsGFootLabel\tdhsGFootCount\tdhsSFootLabel\tGeneName\tDistance\tlogDistance\tExpressionCluster\tPromoterState\tLocation\tExpression\tlogExpression\tExpression2\tlogExpression2\tOverlaps\tChipStates\tChipLabel\tBoundState\n"
out.write(header)
count = 1
index = 0
for cell in cells:

    dhsPeakFile = "DHS_peaks_training/DNASE."+cell+".conservative_s.narrowPeak.gz"
    dhsGenFile = "Gen_footprint_training/"+cell+"_footprints_merged.bed.gz"
    dhsSpecFile= "Spec_footprint_training/footprints_"+cell+"/beds/"+tf_lab+".top.full.bed.gz"
    chipSeqLabs = "Labels/"+tf.upper()+".train.labels."+cell+".bed.gz" 
    #To calculate distance
    geneLocs = "Gene_training/labelled_clusters_"+cell+".out"
    
    #bedfiles
    geneLocs = sq.BedFile(geneLocs, 'Optional')
    
    #tabixFiles
    dhsPeaks = pysam.TabixFile(dhsPeakFile)
    dhsGprints = pysam.TabixFile(dhsGenFile)
    dhsSprints = pysam.TabixFile(dhsSpecFile)
    chipSeqLabs = pysam.TabixFile(chipSeqLabs)
    
    bedFiles = [motif, dhsPeaks, dhsGprints, dhsSprints, geneLocs, chipSeqLabs]
    bedLabels = ["Motif", "dhsPeak", "dhsGFoot", "dhsSFoot", "geneLoc", "chipSeqLabs"]        

    for locs in locations[index:dataSize*count]: #Searching through fasta file
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
            if data == "chipSeqLabs":
                hits = bedFiles[b].fetch(c, s, e)
                if hits is not None:
                    labs = []
                    for o in hits:
                        o = o.split("\t")
                        labs.append(o[3])
                    states = [labs.count("U"), labs.count("A"), labs.count("B")]
                    labels.append(",".join(map(str, states)))
                    if states[2] == 4:
                        state = "0"
                    elif states[1] == 4:
                        state = "1"
                    elif states[0] == 4:
                        state = "2"
                    elif states[0] == 0:
                        if states[1] == 1:
                            state = "3"
                        elif states[1] == 2:
                            state = "4"
                        elif states[1] == 3:
                            state = "5"
                    elif states[0] == 1:
                        if states[1] == 0:
                            state = "12"
                        elif states[1] == 1:
                            state = "11"
                        elif states[1] == 2:
                            state = "6"
                        elif states[1] == 3:
                            state = "7"
                    elif states[0] == 2:
                        if states[1] == 0:
                            state = "13"
                        elif states[1] == 2:
                            state = "8"
                        elif states[1] == 1:
                            state = "9"
                    elif states[0] == 3:
                        state = "10"
                    else:
                        state = "-1"
                    labels.append(state)
                    if (state in positives):
                        labels.append("True")
                    else:
                        labels.append("False")
                else:
                    labels = labels + ["Null"]*3
                
            elif data == "dhsPeak" or data == "dhsSFoot" or data == "dhsGFoot":
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
        
    index = dataSize*count
    count += 1
out.close()




