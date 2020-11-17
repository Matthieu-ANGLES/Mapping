#!/usr/bin/python
# -*- coding : utf-8 -*-

import sys, re

cigar = input(str("Entrer une valeur de cigar :"))
# exemple cigar 3M1I3M1D5M
# autre exemple 312M41I38M1D669M

ext = re.findall('\w',cigar) # re fonction, donne une liste chaque éléments individualisés
print (ext)

valeur=[]
lettre=[]
val=""

for i in range (0,len(ext)) :
    if (ext[i]=='M' or ext[i]=='D' or ext[i]=='I' or ext[i]=='S' or ext[i]=='H' 
                    or ext[i]=="=" or ext[i]=='X'or ext[i]=='N'or ext[i]=='P') :
        lettre.append(ext[i])
        valeur.append(val)
        val=""
    else :
        val=""+val+ext[i]

#print (lettre)
#print (valeur)

for i in range (0, len(lettre)):
    print (valeur[i],lettre[i])

