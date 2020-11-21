# -*- coding : utf-8 -*-

import sys,os,re

# Arguments
pyScriptName = sys.argv[0]
samFileName = sys.argv[1]
fileName, fileExtention = os.path.splitext(samFileName)
samOutputName = fileName.split('/')
samOutputName = samOutputName[-1]+".fasta"

# Check SAM format
with open(samFileName, "r") as sam_file :
    stop = 0
    for i in range (1):
        check = sam_file.readline()
        if not (fileExtention == ".sam") :
            print ("****** SAM extention not explicit ******")
        if not (check.startswith("@")) :
            print ("*********** Error SAM format ***********")
            stop=1

# Parsing SAM File (print Header and split line in list

QNAME=[]
FLAG=[]
RNAME=[]
POS=[]
MAPQ=[]
CIGAR=[]
RNEXT=[]
PNEXT=[]
TLEN=[]
SEQ=[]
QUAL=[]

with open(samFileName, "r") as sam_file, open(samOutputName, "w"):
    for line in sam_file:
        #Check Sam format
        if stop == 1 :
            break

        if line.startswith("@"): # Header motif
            sam_header = line.strip()
            print(sam_header)
        else:
            samLine = line.split("\t") # Fractionnement de la line suivant \t
            QNAME.append(samLine[0])
            FLAG.append(samLine[1])
            RNAME.append(samLine[2])
            POS.append(samLine[3])
            MAPQ.append(samLine[4])
            CIGAR.append(samLine[5])
            RNEXT.append(samLine[6])
            PNEXT.append(samLine[7])
            TLEN.append(samLine[8])
            SEQ.append(samLine[9])
            QUAL.append(samLine[10])

'''
for i in range (0,3) :
    print("\n")
    print (QNAME[i])
    print (FLAG[i])
    print (CIGAR[i])
    print (SEQ[i])
'''







## samOutputName.write(add_line[0] 
