#!/usr/bin/python3
#-*- coding : utf-8 -*-

"""
SamReader

 
"""

__authors__ = ("Benoit Aliaga", "Matthieu Angles")
__contact__ = ("aliaga.benoit@gmail.com", "matthieu.angles@hotmail.fr")
__version__ = "0.0.1"
__date__ = "12/21/2020"



############### IMPORT MODULES ###############

import os, sys, re, getopt


############### FUNCTIONS ############### 

#### Read, store and parse a sam file. ####

def usage():
    """
      This function show the help menu and exit the programm.
    """
    
    sys.stderr.write('''
    
    NAME:
        SamReader

    VERSION:
        0.0.1

    DATE: 
        12/21/2020

    FUNCTION:
        SamReader is a parser program for samtools files (.sam).

    OPTION LIST:
        -h or --help : help information
        -i or --input: input file (.sam)
        -o or --output: output name files

    EXAMPLES:
        python SamReader.py -h # Launch the help.
        python SamReader.py -i <file> # Launch SamReader to analyse a samtools file (.sam).
        python SamReader.py -i <file> -o <name>
        python SamReader.py -i <file> -f <option> # explain

    TROUBLESHOOTINGS:
        If you encounter one or several malfunctions when you execute this software, 
        please contact us (see our e-mail addresses below).

    AUTHORS:
        Benoit ALIAGA (aliaga.benoit@gmail.com)
        Matthieu ANGLES (matthieu.angles@hotmail.fr)

    MORE INFORMATIONS:
        https://github.com/wanaga3166/Mapping
    
    LICENCE:
        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.
        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU General Public License for more details.
        You should have received a copy of the GNU General Public License
        along with this program. If not, see <https://www.gnu.org/licenses/>.
    
    ''')
    sys.exit(-1)

def checkSamFormat(argv): 
    """
      Check only two lines for header and the first line of read
    """

    # Check samFilename (= argv)
    fileName, fileExtention = os.path.splitext(argv)
    if not (fileExtention == ".sam") :
        print ("****** SAM extention not explicit ******")

    sam_header = []
    sam_line = []
    check = 0
    
    with open(argv, "r") as sam_file :
        for line in sam_file:                   
            if line.startswith("@"):    # Check header motif  
                header_line = line.strip()
                sam_header.append(header_line)
                check = 1
            else :
                sam_line.append(line)
                break

    # Check sam_header (Two lines only)
    if (check == 0) or (len(sam_header)< 2):
        print ("*********** Error SAM format ***********")
        sys.exit()

    # Check sam_line (First line only)
    if ("\t" not in sam_line[0]) :
        print ("*********** Error SAM format ***********")
        sys.exit()

    nbElem = 0
    line = sam_line[0]
    line = line.split("\t")
    for elem in line :
        nbElem += 1

    if nbElem < 10 :
        print ("*********** Error SAM format ***********")
        sys.exit()

    return fileName

def openSamHeader(argv):
    """
      Extract header from a sam file. 
    """

    sam_header = []
    
    with open(argv, "r") as sam_file, open("summary_header.txt", "w") as summary_header:
        for line in sam_file:
            if line.startswith("@"): # Header motif
                header_line = line.strip()
                sam_header.append(header_line)  
            else:
                break
        
        for line in sam_header:                
            if line.startswith("@SQ"):
                ligne = line.split("\t")
                for elem in ligne :
                    if elem.startswith("SN"):
                        SN = elem[3:]
                    if elem.startswith("LN"):
                        LN = elem[3:]
                    if elem.startswith("AS"):
                        AS = elem[3:]
                    if elem.startswith("M5"):
                        M5 = elem[3:]
                    if elem.startswith("SP"):
                        SP = elem[3:]
                    if elem.startswith("UR"):
                        UR = elem[3:]

            if line.startswith("@PG"):
                ligne = line.split("\t")
                for elem in ligne :
                    if elem.startswith("ID"):
                        ID = elem[3:]
                    if elem.startswith("PN"):
                        PN = elem[3:]
                    if elem.startswith("VN"):
                        VN = elem[3:]
                    if elem.startswith("CL"):
                        CL = elem[3:]
      
        summary_header.write("======================= HEADER INFORMATION ======================="+"\n")    
        summary_header.write("Reference sequence name : "+SN+"\n")
        summary_header.write("Reference sequence lenght : "+LN+"\n")
        summary_header.write("Programme record identifier : "+ID+"\n") 
        summary_header.write("Programme name : "+PN+"\n")
        summary_header.write("Programme version : "+VN+"\n")
        summary_header.write("Command Line : "+CL+"\n")
        #summary_header.write("=================================================================="+"\n")
        
    return sam_header

def openSam(argv):
    """ Open a sam file and get the lines after the header. 

    Parameters
    ----------
    samtools file (.sam) as argv.

    Returns
    -------
    A list of line from a sam file.
    """

    sam_header = []
    sam_line = []
    
    with open(argv, "r") as sam_file :
        for line in sam_file:
            if line.startswith("@"):
               pass
            else:
                sam_line.append(line)
        return sam_line

#### Parse function ####


def parseSamLine(sam_line):
    """ Parse function of a list which containt the sam file lines. 

    Parameters
    ----------
    sam_line (a list of )

    Returns
    -------
    None
    """

    cpt = 1 # Initialisation d'un compteur
    nbReads = 0
    globalGC = 0
    with open("parse_flag_table.txt", "w") as output_sam_flag, open("outpuTable_cigar.txt", "w") as output_sam_cigar, open("outpuTable_GC_percent.txt", "w") as output_GCpercent:
        for line in sam_line:
            parse = line.split("\t") # Fractionnement de la line suivant \t
            if cpt == 1:
                add_line = [parse[0], parse[1], parse[5], parse[9]] # ajout de flag1 et cigar1 et seq1
                cpt +=1
            elif cpt == 2:
                add_line.append(parse[1]) # ajout de flag2 et de cigar 2 et seq2
                add_line.append(parse[5]) 
                add_line.append(parse[9]) 

                output_sam_flag.write(add_line[0] + ";" + add_line[1] + ";" + add_line[4] + "\n")

                c1 = readCigar(add_line[2])
                c2 = readCigar(add_line[5])
                output_sam_cigar.write(add_line[0] + ";" + add_line[2] + ";" + percentMutation(c1) + add_line[5] + ";" + percentMutation(c2) + "\n")

                GC1 = percentGC(add_line[3])
                output_GCpercent.write(add_line[0] + ";" + str(GC1) + "\n")
                GC2 = percentGC(add_line[6])
                output_GCpercent.write(add_line[0] + ";" + str(GC2) + "\n")

                cpt = 1

#### Calculus functions ####

def flagBinary(flag) :
    """ Convert the flag into binary.

    Parameters
    ----------
    flag from a sam file (samtools) (integer)

    Returns
    -------
    a list of integer.
    """

    flagB = bin(int(flag)) 
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 

    if len(flagB) < 12: # Size adjustement to 12 (normalized size)
        add = 12 - len(flagB) 
        for t in range(add):
            flagB.insert(0,'0')

    return flagB

def unmapped(sam_line):
    """
        Analyse the reads which are unmapped (not paired).
    """
    
    unmapped_count = 0
    
    with open ("only_unmapped.fasta", "a+") as unmapped_fasta, open("summary_unmapped.txt", "w") as summary_file:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1])

            if int(flag[-3]) == 1:
                unmapped_count += 1
                unmapped_fasta.write(toStringOutput(line))

        summary_file.write("================== FLAG INFORMATIONS ===================" + "\n")
        summary_file.write("\n")
        summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") 

        return unmapped_count

def partiallyMapped(sam_line):
    """
        Analyse the read which are only partially mapped.
    """

    partially_mapped_count = 0

    with open ("only_partially_mapped.fasta", "a+") as partillay_mapped_fasta, open("summary_partially_mapped.txt", "w") as summary_file:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1])

            if int(flag[-2]) == 1:
                if col_line[5] != "100M":
                    partially_mapped_count += 1
                    partillay_mapped_fasta.write(toStringOutput(line))

        summary_file.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n") 
        return partially_mapped_count

def pairedMapped1_Partial2(sam_line):
    """
        Analyze the paired read which are mapped for the first one and partially mapped for the second one.
        A fasta file will be written at the end of the script.
    """

    pairedMapped1Partially2_count = 0
    cpt = 1 # Initialisation d'un compteur
    nameRead1 = ''    

    with open ("read1_mapped_read2_partially.fasta", "a+") as mapped1_partially2_fasta, open("summary_paired_mapped_partially.txt", "w") as summary_file:

        # Option 1: the first read is correctly mapped and the second is partially mapped.
        for line in sam_line:
            col_line = line.split("\t") # Fractionnement de la line suivant \t
            flag = flagBinary(col_line[1])

            if cpt == 1:
                if int(flag[-2]) == 1:
                    if col_line[5] != '100M':
                        nameRead1 = col_line[0]
                        mapped1_partially2_fasta.write(toStringOutput(line))
                cpt += 1
            elif cpt == 2:
                #print(cpt)
                if col_line[0] == nameRead1:
                    if col_line[5] == '100M':
                        mapped1_partially2_fasta.write(toStringOutput(line))
                        pairedMapped1Partially2_count += 1
                nameRead1 = ''
                cpt = 1

        # Option 2: the first read is partially mapped and the second read is correctly mapped.
        for line in sam_line:
            col_line = line.split("\t") # Fractionnement de la line suivant \t
            flag = flagBinary(col_line[1])

            if cpt == 1:
                if int(flag[-2]) == 1:
                    if col_line[5] == '100M':
                        nameRead1 = col_line[0]
                        Read1 = line
                cpt +=1
            elif cpt == 2:
                #print("2")
                nameRead2 = col_line[0]
                #print(nameRead2)
                if nameRead2 == nameRead1:
                    #print("It's the same !")
                    if col_line[5] != '100M':
                        pairedMapped1Partially2_count += 1
                        Read2 = line
                        mapped1_partially2_fasta.write(toStringOutput(Read1))
                        mapped1_partially2_fasta.write(toStringOutput(Read2))
                nameRead1 = ''
                nameRead2 = ''
                cpt = 1
                Read1 = ''
                Read2 = ''

        summary_file.write("Total mapped reads and partially mapped reads (paired): " + str(pairedMapped1Partially2_count) + "\n")          
    return pairedMapped1Partially2_count

def mapped1_unmapped2(sam_line):
    """
        Analyze the paired read which are mapped for the first one and unmapped for the second one.
        A fasta file will be written at the end of the script.
    """

    pairedmapped1Unmapped2_count = 0
    cpt = 1 
    nameRead1 = ''

    with open ("read1_mapped_read2_unmapped.fasta", "a+") as mapped1_unmapped2_fasta, open("summary_mapped1_unmmaped2.txt", "a+") as summary_file:

        # Option 1: the first read is correctly mapped and the second read is unmapped.
        for line in sam_line:
            col_line = line.split("\t") # Fractionnement de la line suivant \t
            flag = flagBinary(col_line[1])

            if cpt == 1:
                if (int(flag[-3]) == 1) and not (int(flag[-4]) == 1) : # unmapped
                    nameRead1 = col_line[0]
                    Read1 = line
                cpt += 1
            elif cpt == 2:
                nameRead2 = col_line[0]
                if nameRead1 == nameRead2:
                    Read2 = line
                    pairedmapped1Unmapped2_count += 1
                    mapped1_unmapped2_fasta.write(toStringOutput(Read1))
                    mapped1_unmapped2_fasta.write(toStringOutput(Read2))
                cpt = 1
                Read1 = ''
                Read2 = ''

        # Option 2: the first read is unmapped and the second read is correctly mapped.
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1])

            if cpt == 1:
                nameRead1 = col_line[0]
                Read1 = line
                if (int(flag[-3]) == 0) and (int(flag[-4]) == 1): #
                    cpt += 1
            elif cpt == 2:
                nameRead2 = col_line[0]
                Read2 = line
                if nameRead1 == nameRead2:   
                    if (int(flag[-3]) == 1) and (int(flag[-4]) == 0):
                        Read2 = line
                        pairedmapped1Unmapped2_count += 1
                        mapped1_unmapped2_fasta.write(toStringOutput(Read1))
                        mapped1_unmapped2_fasta.write(toStringOutput(Read2))
                cpt = 1

        summary_file.write("One read is mapped, the other one is not: {0} (paired) {1} (not paired) \n".format(pairedmapped1Unmapped2_count, (pairedmapped1Unmapped2_count * 2)))
    return pairedmapped1Unmapped2_count

def pairedunmapped(sam_line):
    """ Compute the all the paired reads which are unmapped.
    
    Parameters
    ----------
    A line (list) from a sam file (samtools) previous parsed with  

    Returns
    -------
    None
        
    """

    pairedUnmapped1Unmapped2_count = 0
    cpt = 1 
    nameRead1 = ''

    with open("read1_unmapped_read2_unmapped.fasta", "a+") as unmapped1_unmapped2_fasta, open("Summary_paired_unmapped.txt", "a+") as summary_file:

        # The first read is unmapped and the second read is unmapped.
        for line in sam_line:
            col_line = line.split("\t") # Fractionnement de la line suivant \t
            flag = flagBinary(col_line[1])

            if cpt == 1:
                if (int(flag[-3]) == 1) and (int(flag[-4]) == 1) : # unmapped
                    nameRead1 = col_line[0]
                    Read1 = line
                    pairedUnmapped1Unmapped2_count += 1
                cpt+= 1
            elif cpt == 2:
                nameRead2 = col_line[0]
                if nameRead1 == nameRead2: # Control if the read name is the same than the previous line
                    Read2 = line
                    pairedUnmapped1Unmapped2_count += 1
                    unmapped1_unmapped2_fasta.write(toStringOutput(Read1))
                    unmapped1_unmapped2_fasta.write(toStringOutput(Read2))
                cpt = 1
                    
        summary_file.write("The both reads are unmapped: {0} (reads) \n".format(pairedUnmapped1Unmapped2_count))
        summary_file.write("\n")
        summary_file.write("For more information about the flag combinations and their percentage present in this sam file, please check: Final_Flag_table.txt. Be carefull if you use the output option (-o), you changed the name of your file which is: newname_Final_Flag_table.txt.")
    
    return pairedUnmapped1Unmapped2_count

def countFlag():
    """
      Count the flag number in the sam file.
    """
    with open("parse_flag_table.txt", "r") as table_flag, open("count_flag_table.txt", "w") as output_flag:
        dico = {}
        for line in table_flag:
            flag = line.rstrip("\n").rsplit(";")
            inter = flag[1] + "-" + flag[2]
            #inter = flag[1] + "-" + flag[3]
            if inter not in dico.keys():
                dico[inter] = 1
            else:
                dico[inter] += 1
            #print(dico)
        for key in dico.keys():
            output_flag.write(key + ";" + str(dico[key]) + "\n")

def perceFlag(): #total
    """
      This function compute the percentage for each combination flag and return a table with the combinations, 
      their reads numbers and their percentages.
    """

    total = 0
    with open("count_flag_table.txt", "r") as output_sam_flag1 :
        for line in output_sam_flag1:
            number = line.rstrip("\n").split(";")
            numberi = int(number[1])
            total = numberi + total
    
    percList = []
    with open("count_flag_table.txt", "r") as output_sam_flag2, open("Final_Flag_table.txt", "w") as FinalFlag:
        for line2 in output_sam_flag2:
            number2 = line2.rstrip("\n").split(";")
            perc = (int(number2[1]) * 100) / total
            FinalFlag.write(number2[0] + ";" + number2[1] + ";" + str(round(perc,4)) + "\n")
 
def readCigar(cigar): # Make a directory with lists key, value
    """
    Translate a cigar into un dictionnary (Keys = each mutation, Values = number of each mutation)
    Use in parseCigar and tranlate in standar output percentagecin percentMutation
    """
    
    ext = re.findall('\w',cigar) # split cigar 
    key=[]      # arg string in key list
    value=[]    # arg int in value list
    val=""

    for i in range(0,len(ext)) :
        if (ext[i] == 'M' or ext[i] == 'I' or ext[i] == 'D' or ext[i] == 'S' or ext[i] == 'H' 
                        or ext[i] == "N" or ext[i] == 'P'or ext[i] == 'X'or ext[i] == '=') :
            key.append(ext[i])
            value.append(val)
            val = ""
        else :
            val = "" + val + ext[i]
    
    dico = {}
    n = 0
    for k in key:   # Dictionnary contruction in range size lists              
        if k not in dico.keys():    # for each key, insert int value
            dico[k] = int(value[n])   # if key not exist, create and add value
            n += 1
        else:
            dico[k] += int(value[n])  # inf key exist add value
            n += 1
    return dico

def percentMutation(dico):
    """
      Give percentage for each mutation in readCigar dictionnary
      Standardize the values under the same format of all possible mutations
      Use in parseCigar with dicionnary gived by readCigar
    """
        
    totalValue = 0 # Calculated total mutation number
    for v in dico :
        totalValue += dico[v]

    mutList = ['M','I','D','S','H','N','P','X','=']
    res = ""
    for mut in mutList : # Calculated percent of mutation if mut present in the dictionnary, else, percent of mut = 0
        if mut in dico.keys() :
            res += (str(round((dico[mut] * 100) / totalValue, 2)) + ";")
        else :
            res += ("0.00" + ";")
    return res

def globalPercentCigar():
    """
      Global representation of cigar distribution.
    """
    
    with open ("outpuTable_cigar.txt","r") as outpuTable, open("Final_Cigar_table.txt", "w") as FinalCigar:
        nbReads, M, I, D, S, H, N, P, X, Egal = [0 for n in range(10)]

        for line in outpuTable :
            mutValues = line.split(";")
            nbReads += 2
            M += float(mutValues[2])+float(mutValues[12])
            I += float(mutValues[3])+float(mutValues[13])
            D += float(mutValues[4])+float(mutValues[14])
            S += float(mutValues[5])+float(mutValues[15])
            H += float(mutValues[6])+float(mutValues[16])
            N += float(mutValues[7])+float(mutValues[17])
            P += float(mutValues[8])+float(mutValues[18])
            X += float(mutValues[9])+float(mutValues[19])
            Egal += float(mutValues[10])+float(mutValues[20])

        FinalCigar.write("Global cigar mutation observed :"+"\n"
                        +"Alignlent Match : "+str(round(M/nbReads,2))+"\n"
                        +"Insertion : "+str(round(I/nbReads,2))+"\n"
                        +"Deletion : "+str(round(D/nbReads,2))+"\n"
                        +"Skipped region : "+str(round(S/nbReads,2))+"\n"
                        +"Soft Clipping : "+str(round(H/nbReads,2))+"\n"
                        +"Hard Clipping : "+str(round(N/nbReads,2))+"\n"
                        +"Padding : "+str(round(P/nbReads,2))+"\n"
                        +"Sequence Match : "+str(round(Egal/nbReads,2))+"\n"
                        +"Sequence Mismatch : "+str(round(X/nbReads,2))+"\n")

def percentGC(seq):
    """ Compute the GC percentage in a sequence (see the formula below). 
    Formula :  ((G+C) / (A+T+G+C) * 100)

    Parameters
    ----------
    A sequence (seq) extracted from a sam file (samtools).

    Returns
    -------
    a float
    
    """
    countGC = 0
    countAT = 0
    # Count nucleotides :
    for n in seq :
        if (n == 'G') or (n == 'C') :
            countGC += 1
        if (n == 'A') or (n == 'T') :
            countAT += 1
    # Calcul percentage :
    return((countGC/(countAT+countGC))*100)

def countGC():
    """ Count the GC

    This function doesn't have a parameter and write a text file which summarize the GC% in the file.

    Parameters
    ----------
    None

    Returns
    -------
    Nonce
    """
    with open("outpuTable_GC_percent.txt", "r") as table_GCpercent, open("Final_GC_table.txt", "w") as FinalGC:
        dico = {"70-100":0,"60-70":0,"50-60":0,"40-50":0,"30-40":0,"0-30":0}
        nbReads = 0
        for line in table_GCpercent:
            nbReads += 1
            parse = line.rstrip("\n").rsplit(";")
            percentGC = float(parse[1])
            if percentGC >= 70 and percentGC <= 100 :
                dico["70-100"] += 1
            if percentGC >= 60 and percentGC < 70 :
                dico["60-70"] += 1
            if percentGC >= 50 and percentGC < 60 :
                dico["50-60"] += 1
            if percentGC >= 40 and percentGC < 50 :
                dico["40-50"] += 1
            if percentGC >= 30 and percentGC < 40 :
                dico["30-40"] += 1
            if percentGC >= 0 and percentGC < 30 :
                dico["0-30"] += 1
        
        dicoPercent = {"70-100":0,"60-70":0,"50-60":0,"40-50":0,"30-40":0,"0-30":0}
        for k in dico :
            dicoPercent[k] = (dico[k]*100)/nbReads

        FinalGC.write("GC%: >70%: "+str(dico["70-100"])+" ("+str(round(dicoPercent["70-100"],2))+"%)\n"
                        +"GC%: >60%: "+str(dico["60-70"])+" ("+str(round(dicoPercent["60-70"],2))+"%)\n"
                        +"GC%: >50%: "+str(dico["50-60"])+" ("+str(round(dicoPercent["50-60"],2))+"%)\n"
                        +"GC%: >40%: "+str(dico["40-50"])+" ("+str(round(dicoPercent["40-50"],2))+"%)\n"
                        +"GC%: >30%: "+str(dico["30-40"])+" ("+str(round(dicoPercent["30-40"],2))+"%)\n"
                        +"GC%: <30%: "+str(dico["0-30"])+" ("+str(round(dicoPercent["0-30"],2))+"%)\n")

#### Output functions ####

def toStringOutput(read_line):
    """
    Write function for output
    """
    line = read_line.split("\t")
    qname = str(line[0])
    flag = str(line[1])
    mapQ = str(line[4])
    cigar = str(line[5])
    seq = str(line[9])
    return ("> " + qname + " | flag:"+flag+" | cigar:"+ cigar +" | mapQ:"+mapQ+" | GC:"+str(round(percentGC(seq),2))+"%"+"\n"+seq+"\n"+"\n")

def computeSummary(fileName):
    """ computeSummary will write a summary file.

    This function need the name of the flag

    Parameters
    ----------
    file name 

    Returns
    -------
    None

    """

    with open("summary_header.txt", "r") as F_Head, open("Final_Cigar_table.txt", "r") as F_Cigar, open("Final_GC_table.txt", "r") as F_GC, open("summary_unmapped.txt", "r") as F_UN, open("summary_partially_mapped.txt", "r") as F_PA, open("summary_paired_mapped_partially.txt", "r") as F_M1P2, open("summary_mapped1_unmmaped2.txt", "r") as F_M1U2, open("Summary_paired_unmapped.txt", "r") as F_PU, open("summary.txt", "a+") as F_Sum:
        F_Sum.write("==================================================================\n"+
                    "  >  Summary " + fileName + "\n")
        for line in F_Head :
            l = line.rstrip("\n")
            F_Sum.write(line)
        F_Sum.write("\n")
        F_Sum.write("\n")
        for line in F_UN:
            l =line.rstrip("\n")
            F_Sum.write(line)
        for line in F_PA :
            l = line.rstrip("\n")
            F_Sum.write(line)
        for line in F_M1P2:
            l = line.rstrip("\n")
            F_Sum.write(line)
        for line in F_M1U2:
            l = line.rstrip("\n")
            F_Sum.write(line)
        for line in F_PU:
            l = line.rstrip("\n")
            F_Sum.write(line)
        F_Sum.write("\n")
        F_Sum.write("======================= PARSE INFORMATIONS =======================\n")
        #for line in F_flag :
        #    l = line.rstrip("\n")
        #    F_Sum.write(line)
        F_Sum.write("\n")    
        for line in F_Cigar :
            l = line.rstrip("\n")
            F_Sum.write(line)
        F_Sum.write("\n")
        for line in F_GC :
            l = line.rstrip("\n")
            F_Sum.write(line)

def outFile(newname): 
    """ Output file name

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    os.remove("summary_header.txt")
    os.remove("parse_flag_table.txt")
    os.remove("count_flag_table.txt")
    os.remove("outpuTable_cigar.txt")
    os.remove("outpuTable_GC_percent.txt")
    os.remove("Final_Cigar_table.txt")
    os.remove("Final_GC_table.txt")
    os.remove("summary_unmapped.txt")
    os.remove("summary_partially_mapped.txt")
    os.remove("summary_paired_mapped_partially.txt")
    os.rename("Final_Flag_table.txt", newname + "_final_Flag_table.txt")
    os.rename("summary.txt", newname + "_summary.txt")

def summaryPrint():
    """ Print the summary.txt in the terminal.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    with open("summary.txt", "r") as summary_read:
        for line in summary_read:
            l = line.rstrip("\n")
            print(l)

#### Main function ####

def main(argv):
    """ Main function.

    This function use the getopt module for the argument options. 
    For option list and their command line, please watch the usage function for more 
    information or README file (README.txt or README.md in github page.

    Parameters
    ----------
    None

    Returns
    -------
    None


    """
    short_options = "hi:o:f:"
    long_options = ["help", "input=", "output="]
    
    try:
        # Parsing argument
        arguments, values = getopt.getopt(argv, short_options, long_options)

        # Evaluate given options
        for current_argument, current_value in arguments:
            if current_argument in ("-h", "--help"):
                usage()
                
            if current_argument in ("-i", "--input"):
                print("Check Sam format.")
                fileName = checkSamFormat(current_value)
                
                print("Read the file.")
                resSamHeader = openSamHeader(current_value)
                resSam = openSam(current_value)

                print("Parsing.")
                parseSamLine(resSam)
                #parseCigar(resSam)
                #flagBin2(resSam)

                print("Analyze the only the unmapped reads.")
                unmapped(resSam)

                print("Analyze the only the partially mapped reads.")
                partiallyMapped(resSam)

                print("Analyze the reads which are mapped in one read and partially mapped in the other read.")
                pairedMapped1_Partial2(resSam)

                print("Analyze the reads which are mapped in one read and unmapped in the other read.")
                mapped1_unmapped2(resSam)

                print("Analyze the reads which are paired and unmapped.")
                pairedunmapped(resSam)

                print("Count the number of read for each flag combinations.")
                countFlag()

                print("Calculate the percentage for each flag combination.")
                #numberReads = total()
                #percFlag(numberReads)
                perceFlag()

                print("Calculate the global percentage mutation of cigars")
                globalPercentCigar()

                print("Calculate the percentage of GC contain")
                countGC()

                print("Compute Summary")
                computeSummary(fileName)
                summaryPrint()

            if current_argument in ("-o", "--output"):
                print("Ouput the file.")
                # voir pour options ?!!
                outFile(current_value) #  mis en dernier car affichage et renommage summary plus bas

        print("Analyze finished.")
                
    except getopt.error as err:
        # Output error, and return with an error code
        print (str(err))
        usage()
        sys.exit(2)

############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main(sys.argv[1:])
