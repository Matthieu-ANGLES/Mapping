#!/usr/bin/python3
# -*- coding : utf-8 -*-

import sys,re

# Arguments
pyScriptName = sys.argv[0]
samFileName = sys.argv[1]
samOutputName = sys.argv[2]

# Fonctions readCigar : Make a distionary (key=Op, value=nb)

def readCigar (cigar):
    ext = re.findall('\w',cigar) # re fonction, donne une liste chaque éléments individualisés
    value=[]
    key=[]
    val=""

    for i in range (0,len(ext)) :
        if (ext[i]=='M' or ext[i]=='D' or ext[i]=='I' or ext[i]=='S' or ext[i]=='H' 
                        or ext[i]=="=" or ext[i]=='X'or ext[i]=='N'or ext[i]=='P') :
            key.append(ext[i])
            value.append(val)
            val=""
        else :
            val=""+val+ext[i]
    
    dico={}
    n=0
    for k in key:                   
        if k not in dico.keys():
            dico[k]=value[n]
            n+=1
        else:
            dico[k]+=value[n]
            n+=1
    return (dico)
  
# Fonctions readCigar : Tanslate into string a distionary
def textCigar (dico):
    res=""
    for i in dico.keys():
        #print(i)
        if i == 'M':
            res += " Alignment mapped "+dico[i]
        elif i=='D':
            res += " Deletion(s) "+dico[i]
        elif i=='I':
            res += " Insertion(s) "+dico[i]
        elif i=='S':
            res += " Soft clipping "+dico[i]
        elif i=='H':
            res += " Hard clipping "+dico[i]
        elif i=='=':
            res += " Sequence match "+dico[i]
        elif i=='X':
            res += " Sequence mismatch "+dico[i]
        elif i=='N':
            res += " skipped region "+dico[i]
        elif i=='P':
            res += " Padding (silent del) "+dico[i]

    return(res)


# Parsing SAM File (print Header and write ReadNames, cigar R1, R2 and text on output file)

with open(samFileName, "r") as sam_file, open(samOutputName, "w") as output_sam_cigar:
    cpt = 1 # Initialisation d'un compteur
    for line in sam_file:
        if line.startswith("@"): # A modifier
            sam_header = line.strip()
            print(sam_header)
        else:
            #cpt = 1 
            cigar = line.split("\t") # Fractionnement de la line suivant \t
            if cpt == 1:
                add_line = [cigar[0], cigar[5]]
                #print(add_line)
                cpt +=1
            elif cpt == 2:
                add_line.append(cigar[5]) # ajout de cigar2
                new_line = " \t ".join([add_line[0], add_line[1], add_line[2]])
                #print(new_line)
                c1=readCigar(add_line[1])
                c2=readCigar(add_line[2])
                output_sam_cigar.write(add_line[0] + ";" + add_line[1] + ";" + textCigar(c1) + ";" 
                        + add_line[2] + ";" + textCigar(c2) + ";" + "\n") #str(c1) + str(c2)+
                cpt = 1