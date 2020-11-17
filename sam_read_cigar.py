#!/usr/bin/python3
# -*- coding : utf-8 -*-

import sys,re

# Arguments
pyScriptName = sys.argv[0]
samFileName = sys.argv[1]
samOutputName = sys.argv[2]

def readCigar (cigar):
#cigar = input(str("Entrer une valeur de cigar :"))
# exemple cigar 3M1I3M1D5M
# autre exemple 312M41I38M1D669M

    ext = re.findall('\w',cigar) # re fonction, donne une liste chaque éléments individualisés
    #print (ext)

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
  

# Parsing fichier SAM génère un fichier txt en sortie (nomRead;cigarR1;cigarR2)

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
                print(new_line)
                c1=readCigar(add_line[1])
                c2=readCigar(add_line[2])
                print(c1)
                print(c2)
                output_sam_cigar.write(add_line[0] + ";" + add_line[1] + ";" + add_line[2] 
                        +str(c1) + str(c2)+"\n") #+ c1 + c2+
                cpt = 1

