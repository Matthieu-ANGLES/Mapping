#!/usr/bin/python
#-*- coding : utf-8 -*-

"""
SamReader

SamReader is a parser program for samtool files (.sam). 

How to launch SamReader?

./SamReader.py argument 1 argument 2

argument 1: samtool file (.sam)
argument 2: output table
 
"""

__authors__ = ("Benoit Aliaga", "Matthieu Angles")
__contact__ = ("aliaga.benoit@gmail.com", "matthieu.angles@hotmail.fr")
__version__ = "0.0.1"
__date__ = "11/14/2020"

import os, sys, re, getopt
from FlagDescription import FLAG

def usage():
    sys.stderr.write('''
    
    NAME:
         SamReader

    VERSION:
         0.0.1

    FUNCTION: 
         SamReader is a parser program for samtools files (.sam). 

    OPTION LIST:
         -h or --help : help information
         -i or --input: input file (.sam)
         -o or --output: output name files

    EXAMPLES:
         python SamReader.py -h # Launch the help
         python SamReader.py -i <file> # Launch SamReader to analyse a samtools file (.sam).
         python SamReader.py -i <file> -o <name>

    AUTHORS:
         Benoit ALIAGA (aliaga.benoit@gmail.com)
         Matthieu ANGLES (matthieu.angles@hotmail.fr)

    MORE INFORMATIONS:
         https://github.com/wanaga3166/Mapping
    
    ''')
    sys.exit(-1)

def openSam(argv):
    """
    Read a sam file. 
    """

    sam_header = []
    sam_line = []
    
    with open(argv, "r") as sam_file:
        for line in sam_file:
            if line.startswith("@"):
                header_line = line.strip()
                sam_header.append(header_line)
            else:
                #print("Line en cours")
                sam_line.append(line)
        return sam_line
            
def outFile(argv):
    """
    Output file name
    """
    os.remove("table_flag.txt")
    os.rename("Final_Flag_table.txt", argv)
    
def parseFlag(sam_line):
    """
    Docstring
    """
    cpt = 1 # Initialisation d'un compteur
    with open("table_flag.txt", "w") as output_sam_flag:
        for line in sam_line:
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

def countFlag():
    """
    Docstring
    """
    with open("table_flag.txt", "r") as table_flag, open("Final_Flag_table.txt", "w") as output_flag:
        dico = {}
        for line in table_flag:
            flag = line.rstrip("\n").rsplit(";")
            inter = flag[1] + "-" + flag[2]
            if inter not in dico.keys():
                dico[inter] = 1
            else:
                dico[inter] += 1
            #print(dico)
        for key in dico.keys():
            output_flag.write(key + ";" + str(dico[key]) + "\n")

def parseCigar():
    """
    Docstring
    """
    pass

def parseMAPQ():
    """
    Docstring
    """
    pass

def main(argv):
    """
    Docstring
    """
    short_options = "hi:o:"
    long_options = ["help", "input=", "output="]
    
    try:
        # Parsing argument
        arguments, values = getopt.getopt(argv, short_options, long_options)

        # Evaluate given options
        for current_argument, current_value in arguments:
            if current_argument in ("-h", "--help"):
                usage()
            if current_argument in ("-i", "--input"):
                print("Read the file")
                resSam = openSam(current_value)
                parseFlag(resSam)
                countFlag()
            if current_argument in ("-o", "--output"):
                print("Ouput the file")
                outFile(current_value)
    except getopt.error as err:
        # Output error, and return with an error code
        print (str(err))
        usage()
        sys.exit(2)
        
if __name__ == "__main__":
    main(sys.argv[1:])
