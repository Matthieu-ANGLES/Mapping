#!/usr/bin/python
# -*- coding : utf-8 -*-

import os, sys,re

# Arguments
pyScriptName = sys.argv[0]
samFileName = sys.argv[1]
samOutputName = sys.argv[2]

# Read sam file
with open(samFileName, "r") as sam_file, open(samOutputName, "w") as output_sam_flag:
    cpt = 1 # Initialisation d'un compteur
    for line in sam_file:
        if line.startswith("@"): # A modifier
            sam_header = line.strip()
            print(sam_header)
        else:
            #cpt = 1 
            flag = line.split("\t") # Fractionnement de la line suivant \t
            if cpt == 1:
                add_line = [flag[0], flag[1]]
                #print(add_line)
                cpt +=1
            elif cpt == 2:
                add_line.append(flag[1]) # ajout de flag2
                new_line = " \t ".join([add_line[0], add_line[1], add_line[2]])
                #print(new_line)
                output_sam_flag.write(add_line[0] + ";" + add_line[1] + ";" + add_line[2] + "\n")
                cpt = 1

