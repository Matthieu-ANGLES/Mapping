#!/usr/bin/python3
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
__date__ = "11/21/2020"

import os, sys, re, getopt
#from FlagDescription import FLAG

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

def checkSam(sam_header):
    """
    docstring
    """
    pass



def openSam(argv):
    """
    Read a sam file. 
    """
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

    with open(argv, "r") as sam_file :
        for line in sam_file:
            #Check Sam format
            #if stop == 1 :
            #    break

            if line.startswith("@"): # Header motif
                # checkSam()
                sam_header = line.strip()
                print(sam_header)
                # parse header
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


        superList = []
        superList.append(QNAME)
        superList.append(FLAG)
        superList.append(RNAME)
        superList.append(POS)
        superList.append(MAPQ)
        superList.append(CIGAR)
        superList.append(RNEXT)
        superList.append(PNEXT)
        superList.append(TLEN)
        superList.append(SEQ)
        superList.append(QUAL)
        #print(superList)
    #keySam = {}
    
    #keySam["qname"] = QNAME
    #keySam["flag"] = FLAG
    #keySam["rname"] = RNAME
    #keySam["pos"] = POS
    #keySam["mapq"] = MAPQ
    #keySam["cigar"] = CIGAR 
    #keySam["rnext"] = RNEXT
    #keySam["pnext"] = PNEXT
    #keySam["tlen"] = TLEN
    #keySam["seq"] = SEQ
    #keySam["qual"] = QUAL

    return superList
            
def outFile(argv):
    """
    Output file name
    """
    os.remove("parse_flag_table.txt")
    os.remove("count_flag_table.txt")
    os.rename("Final_Flag_table.txt", argv)
    
def parseFlagCigar(superList):
    """
    Function which create a table with flag and cigar statistics.
    """
    # Parsing SAM File (print Header and write in output file ReadNames;cigarR1,percent of each R1 mutations;cigarR2,percent of each R2 mutations on output file)
    with open(sam_file, "r") as sam_file, open(samOutputName, "w") as output_sam_cigar:
        cpt = 1 # Paired-and read counter
        for line in sam_file:
            if line.startswith("@"):
                sam_header = line.strip()
                print(sam_header)
            else:
                readLine = line.split("\t") # Split args line
                if cpt == 1:
                    add_line = [readLine[0],readLine[1],readLine[5]] # Creation line with Read and cigar R1
                    cpt +=1
                elif cpt == 2:
                    add_line.append(readLine[1]) # Add flag R2
                    add_line.append(readLine[5]) # Add cigar R2 into line
                    c1=readCigar(add_line[3])
                    c2=readCigar(add_line[5])
                    print (add_line)
                    '''
                    output_sam_cigar.write(add_line[0] # Write in output file Read;cigarR1;percent all mutations represent;cigarR2;percent all mutations represent;
                                    +";"+add_line[1]+";"+percentMutation(c1)
                                    +add_line[2]+";"+percentMutation(c2)+"\n")
                    '''
                    cpt = 1




def countFlag():
    """
    Docstring
    """
    with open("parse_flag_table.txt", "r") as table_flag, open("count_flag_table.txt", "w") as output_flag:
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

def total():
    with open("count_flag_table.txt", "r") as output_sam_flag:
        total = 0
    
        for line in output_sam_flag:
            #print(line)
            number = line.rstrip("\n").split(";")
            numberi = int(number[1])
            #print(numberi)
            total = numberi + total
        #print(total)
        return total

def percFlag(total):
    percList = []
    with open("count_flag_table.txt", "r") as output_sam_flag, open("Final_Flag_table.txt", "w") as FinalFlag:
        for line in output_sam_flag:
            number2 = line.rstrip("\n").split(";")
            perc = (int(number2[1])*100) / total
            newline = ";".join([number2[0],number2[1],str(round(perc,4))])
            FinalFlag.write(number2[0] + ";" + number2[1] + ";" + str(round(perc,4)) + "\n")
            
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
                print("Read the file.")
                resSam = openSam(current_value)
                parseFlagCigar(resSam)

                print("Parse the flag.")
                #parseFlag(resSam)

                print("Count the number of read for each flag combinations.")
                #countFlag()

                print("Calculate the percentage for each flag combination.")
                #numberReads = total()
                #percFlag(numberReads)
                
            if current_argument in ("-o", "--output"):
                print("Ouput the file.")
                outFile(current_value)

        print("Analyse finished.")
                
    except getopt.error as err:
        # Output error, and return with an error code
        print (str(err))
        usage()
        sys.exit(2)
        
if __name__ == "__main__":
    main(sys.argv[1:])


