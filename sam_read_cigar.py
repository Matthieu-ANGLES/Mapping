#!/usr/bin/python3
# -*- coding : utf-8 -*-

import sys,re

# Arguments
pyScriptName = sys.argv[0]
samFileName = sys.argv[1]
samOutputName = sys.argv[2]

# Functions readCigar : Make a distionary (key=mut, value=nb)
def readCigar (cigar): # Make a directory with lists key, value
    ext = re.findall('\w',cigar) # split cigar 
    key=[]      # arg string in key list
    value=[]    # arg int in value list
    val=""

    for i in range (0,len(ext)) :
        if (ext[i]=='M' or ext[i]=='I' or ext[i]=='D' or ext[i]=='S' or ext[i]=='H' 
                        or ext[i]=="N" or ext[i]=='P'or ext[i]=='X'or ext[i]=='=') :
            key.append(ext[i])
            value.append(val)
            val=""
        else :
            val=""+val+ext[i]
    
    dico={}
    n=0
    for k in key:   # Dictionnary contruction in range size lists              
        if k not in dico.keys():    # for each key, insert int value
            dico[k]=int(value[n])   # if key not exist, create and add value
            n+=1
        else:
            dico[k]+=int(value[n])  # inf key exist add value
            n+=1
    return (dico)

# Functions readCigar : Give percentage for each mutation in readCigar dictionnary
def percentMutation (dico):

    totalValue = 0 # Calculated total mutation number
    for v in dico :
        totalValue += dico[v]

    mutList = ['M','I','D','S','H','N','P','X','=']
    res = ""
    for mut in mutList : # Calculated percent of mutation if mut present in the dictionnary, else, percent of mut = 0
        if mut in dico.keys() :
            res += (mut+str(round((dico[mut]*100)/totalValue,2))+";")
        else :
            res += (mut+"0.00"+";")
    return (res)

# Parsing SAM File (print Header and write in output file ReadNames;cigarR1,percent of each R1 mutations;cigarR2,percent of each R2 mutations on output file)
with open(samFileName, "r") as sam_file, open(samOutputName, "w") as output_sam_cigar:
    cpt = 1 # Paired-and read counter
    for line in sam_file:
        if line.startswith("@"):
            sam_header = line.strip()
            print(sam_header)
        else:
            readLine = line.split("\t") # Split args line
            if cpt == 1:
                add_line = [readLine[0], readLine[5]] # Creation line with Read and cigar R1
                cpt +=1
            elif cpt == 2:
                add_line.append(readLine[5]) # Add cigar R2 into line
                c1=readCigar(add_line[1])
                c2=readCigar(add_line[2])
                output_sam_cigar.write(add_line[0] # Write in output file Read;cigarR1;percent all mutations represent;cigarR2;percent all mutations represent;
                                    +";"+add_line[1]+";"+percentMutation(c1)
                                    +add_line[2]+";"+percentMutation(c1)+"\n")
                cpt = 1