#!/usr/bin/python3
# -*- coding : utf-8 -*-

import sys,re

# Arguments
pyScriptName = sys.argv[0]
samFileName = sys.argv[1]
samOutputName = sys.argv[2]


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
                print(add_line)
                cpt +=1
            elif cpt == 2:
                add_line.append(cigar[5]) # ajout de cigar2
                new_line = " \t ".join([add_line[0], add_line[1], add_line[2]])
                print(new_line)
                output_sam_cigar.write(add_line[0] + ";" + add_line[1] + ";" + add_line[2] + "\n")
                cpt = 1
