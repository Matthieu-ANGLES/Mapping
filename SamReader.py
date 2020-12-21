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
      Check only two ligne for header and the first line of read
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

    return (fileName)

def openSamHeader(argv):
    """
      Extract header from a sam file. 
    """

    sam_header = []
    
    with open(argv, "r") as sam_file, open("Summary.txt", "w") as summary_header:
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

        summary_header.write("********************************************************" + "\n")
        summary_header.write("*" + "                                                      " +"*" +  "\n")
        summary_header.write("*" + "                    SUMMARY                           " +"*" +  "\n")
        summary_header.write("*" + "                                                      " +"*" +  "\n")
        summary_header.write("********************************************************" + "\n")
        summary_header.write("\n")
        summary_header.write("================== HEADER INFORMATION ==================" + "\n")
        summary_header.write("\n")
        summary_header.write("Reference sequence name : " + SN + "\n")
        summary_header.write("Reference sequence lenght : "+ LN + "\n")
        summary_header.write("Programme record identifier : " + ID + "\n") 
        summary_header.write("Programme name : " + PN + "\n")
        summary_header.write("Programme version : " + VN + "\n")
        summary_header.write("Command Line : " + CL + "\n")
        summary_header.write("\n")
        
    return sam_header

def openSam(samFile):
    """
      Open a sam file and get the lines after the header.   
    """

    sam_header = []
    sam_line = []
    
    with open(samFile, "r") as sam_file :
        for line in sam_file:
            if line.startswith("@"):
               pass
            else:
                sam_line.append(line)
        return sam_line

    
#### Parse function ####

def parseSam(sam_line):
    """
        Docstring
    """

    cpt = 1 # Start the counter at 1

    with open("parse_flag_table.txt", "w") as output_sam_flag:
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

def flagBinary(flag):
    """
        Docstring
    """
    flagB = bin(int(flag)) 
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 

    if len(flagB) < 12: # Size ajustement to 12 (normalized size)
        add = 12 - len(flagB) 
        for t in range(add):
            flagB.insert(0,'0')

    return flagB

def unmapped(sam_line):
    """
        Analyse the reads which are unmapped (not paired).
    """
    
    unmapped_cout = 0
    
    with open ("only_unmapped.fasta", "a+") as unmapped_fasta:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1])

            if int(flag[-3]) == 1:
                unmapped_cout += 1
                unmapped_fasta.write(writeFasta(line))

        print(unmapped_cout)
        return unmapped_cout

def partiallyMapped(sam_line):
    """
        Analyse the read which are only partially mapped.
    """

    partially_mapped_count = 0

    with open ("only_partially_mapped.fasta", "a+") as partillay_mapped_fasta:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1])

            if int(flag[-2]) == 1:
                if col_line[5] != "100M":
                    partially_mapped_count += 1
                    partillay_mapped_fasta.write(writeFasta(line))

        print(partially_mapped_count)
        return partially_mapped_count

def pairedMapped1_Partial2(sam_line):
    """
        Analyze the paired read which are mapped for the first one and partially mapped for the second one.
        A fasta file will be written at the end of the script.
    """

    pairedMapped1Partially2_count = 0
    cpt = 1 # Initialisation d'un compteur
    nameRead1 = ''    

    with open ("read1_mapped_read2_partially.fasta", "a+") as mapped1_partially2_fasta:

        # Option 1: the first read is correctly mapped and the second is partially mapped.
        for line in sam_line:
            col_line = line.split("\t") # Fractionnement de la line suivant \t
            flag = flagBinary(col_line[1])

            if cpt == 1:
                if int(flag[-2]) == 1:
                    if col_line[5] != '100M':
                        nameRead1 = col_line[0]
                        mapped1_partially2_fasta.write(writeFasta(line))
                cpt += 1
            elif cpt == 2:
                #print(cpt)
                if col_line[0] == nameRead1:
                    if col_line[5] == '100M':
                        mapped1_partially2_fasta.write(writeFasta(line))
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
                        mapped1_partially2_fasta.write(writeFasta(Read1))
                        mapped1_partially2_fasta.write(writeFasta(Read2))
                nameRead1 = ''
                nameRead2 = ''
                cpt = 1
                Read1 = ''
                Read2 = ''
                    
    print("Paired reads", pairedMapped1Partially2_count)
    #return mapped1_partially2

def mapped1_unmapped2(sam_line):
    """
        Analyze the paired read which are mapped for the first one and unmapped for the second one.
        A fasta file will be written at the end of the script.
    """

    pairedMapped1Unmapped2_count = 0
    cpt = 1 # Initialisation d'un compteur
    nameRead1 = ''

    with open ("read1_mapped_read2_unmapped.fasta", "a+") as mapped1_unmapped2_fasta:

        # Option 1: the first read is correctly mapped and the second read is unmapped.
        for line in sam_line:
            col_line = line.split("\t") # Fractionnement de la line suivant \t
            flag = flagBinary(col_line[1])

            if cpt == 1:
                if int(flag[-2]) == 1: # Correctly mapped
                    nameRead1 = col_line[0]
                    Read1 = line
                    print(nameRead1)
                    #print(Read1)
                cpt += 1
            elif cpt == 2:
                nameRead2 = col_line[0]
                if nameRead2 == nameRead1:
                    print(nameRead2)
                    print("It's the same !")
                    if int(flag[-3]) == 1:
                        Read2 = line
                        print(Read2)
                        #mapped1_unmapped2_fasta.write(writeFasta(Read1))
                        #mapped1_unmapped2_fasta.write(writeFasta(Read2))               
                nameRead1 = ''
                nameRead2 = ''
                cpt = 1
                Read1 = ''
                Read2 = ''

        # Option 2: the first read is unmapped and the second read is correctly mapped.

def flagBin(sam_line):
    """
        Provides binary traduction of the Flag.
        Fasta outputs options for Read Unmapped or others informations.
    """

    with open("Summary.txt", "a+") as summary_file:
        flag_output = {}
        read_count = 0
        for l in sam_line: 
            read_count += 1       
            line = l.split("\t")
            flag = str(line[1]) # Collect the flag for each line
        
            binflag  = bin(int(flag))
            binflag = binflag[2:] # Remove '0b' Example: '0b1001101' > '1001101'
            binflag = list(binflag) 

            if len(binflag) < 12: # Size ajustement to 12 (normalized size)
                add = 12 - len(binflag) 
                for t in range(add):
                    binflag.insert(0,'0')

            # Flag traduction
            if 1 == int(binflag[-1]): # "Read paired"

                # Count the read paired flag
                if "Read paired" not in flag_output.keys():
                    flag_output["Read paired"] = 1
                else:
                    flag_output["Read paired"] += 1

            if 1 == int(binflag[-2]): # "Read mapped in proper pair"
                if "Read mapped in proper pair" not in flag_output.keys():
                    flag_output["Read mapped in proper pair"] = 1
                else:
                    flag_output["Read mapped in proper pair"] += 1
                  
            if (1 == int(binflag[-3])): # "Read unmapped"
                #with open ("Reads_unmapped_only.txt","a+") as outputUMO:
                #outputUMO.write(toStringOutput(l))
                if "Read unmapped" not in flag_output.keys():
                    flag_output["Read unmapped"] = 1
                else:
                    flag_output["Read unmapped"] += 1
                
            if 1 == int(binflag[-4]): # "Mate unmapped"
                if "Mate unmapped" not in flag_output.keys():
                    flag_output["Mate unmapped"] = 1
                else:
                    flag_output["Mate unmapped"] += 1
                
            if 1 == int(binflag[-5]): # "Read reverse strand"
                if "Read reverse strand" not in flag_output.keys():
                    flag_output["Read reverse strand"] = 1
                else:
                    flag_output["Read reverse strand"] += 1
                
            if 1 == int(binflag[-6]): # "Mate reverse strand"
                if "Mate reverse strand" not in flag_output.keys():
                    flag_output["Mate reverse strand"] = 1
                else:
                    flag_output["Mate reverse strand"] += 1
            
            if 1 == int(binflag[-7]): # "First in pair"
                if "First in pair" not in flag_output.keys():
                    flag_output["First in pair"] = 1
                else:
                    flag_output["First in pair"] += 1
                
            if 1 == int(binflag[-8]): # "Second in pair"
                if "Second in pair" not in flag_output.keys():
                    flag_output["Second in pair"] = 1
                else:
                    flag_output["Second in pair"] += 1

            if 1 == int(binflag[-9]): # "Not primary alignment"
                if "Not primary alignment" not in flag_output.keys():
                    flag_output["Not primary alignment"] = 1
                else:
                    flag_output["Not primaxry alignment"] += 1
                
            if 1 == int(binflag[-10]): # "Read fails platform/vendor quality checks"
                if "Read fails platform/vendor quality checks" not in flag_output.keys():
                    flag_output["Read fails platform/vendor quality checks"] = 1
                else:
                    flag_output["Read fails platform/vendor quality checks"] += 1
                
            if 1 == int(binflag[-11]): # "Read is PCR or optical duplicate"
                #with open ("Reads_is_PCR_or_optical_duplicate.txt","a+") as outputOD:
                #   outputOD.write(toStringOutput(l))
                if "Read is PCR or optical duplicate" not in flag_output.keys():
                    flag_output["Read is PCR or optical duplicate"] = 1
                else:
                    flag_output["Read is PCR or optical duplicate"] += 1
                
            if 1 == int(binflag[-12]): # "Supplementary alignment"
                if "Supplementary alignment" not in flag_output.keys():
                    flag_output["Supplementary alignment"] = 1
                else:
                    flag_output["Supplementary alignment"] += 1
                    # with open ("Supplementary_alignment.txt","a+") as outputSA:
                    #    outputOD.write(toStringOutput(l))
        
        summary_file.write("================== FLAG INFORMATIONS ===================" + "\n")
        summary_file.write("\n")
        summary_file.write("Total read counted: " + str(read_count) + "\n") 
        for k in flag_output:
            summary_file.write(k + ": " + str(flag_output[k]) + "\n")
        summary_file.write("\n")

        return flag_output

                
#### Calculus functions ####

def countFlag():
     """
     Count the flag number in the sam file.
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

def percFlag(total):
     """
     Docstring
     """
     percList = []
     with open("count_flag_table.txt", "r") as output_sam_flag, open("Final_Flag_table.txt", "w") as FinalFlag:
         for line in output_sam_flag:
             number2 = line.rstrip("\n").split(";")
             perc = (int(number2[1]) * 100) / total
             newline = ";".join([number2[0],number2[1],str(round(perc,4))])
             FinalFlag.write(number2[0] + ";" + number2[1] + ";" + str(round(perc,4)) + "\n")

def total():
     """
     Docstring
     """
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

             
#### Output functions ####

def writeFasta(read_line):
    """
    Write function for output
    """

    line = read_line.split("\t")
    qname = str(line[0])
    flag = str(line[1])
    mapq = str(line[4])
    cigar = str(line[5])
    seq = str(line[9])
    return ("> " + qname + " | flag:" + flag + " | cigar: " + cigar + " | mapQ: " + mapq + "\n" + seq + "\n" + "\n")

def readSummary():
    with open ("Summary.txt", "r") as final_summary:
        for line in final_summary:
            l = line.rstrip("\n")
            print (l)

def outFile():
    """
    Output file name
    """
    os.remove("parse_flag_table.txt")
    os.remove("count_flag_table.txt")
    #os.rename("Final_Flag_table.txt", argv)
            
#### Main function ####

def main(argv):
    """
      Main function.

      This function use the getopt module for the argument options. 

      For option list and their command line, please watch the usage function for more 
      information or README file (README.txt or README.md in github page.
    """

    short_options = "hi:o:f:"
    long_options = ["help", "input=", "output=", "fasta="]
    
    try:
        # Parsing argument
        arguments, values = getopt.getopt(argv, short_options, long_options)

        # Evaluate given options
        for current_argument, current_value in arguments:
            if current_argument in ("-h", "--help"):
                usage()
                
            if current_argument in ("-i", "--input"):
                print("Check Sam format.")
                checkSamFormat(current_value)
                
                print("Read the file.")
                openSamHeader(current_value)
                resSam = openSam(current_value)

                print("Parse the file.")
                test1 = flagBin(resSam)
                parseSam(resSam)
                
                print("Analyze the only the unmapped reads.")
                unmapped(resSam)

                print("Analyze the only the partially mapped reads.")
                partiallyMapped(resSam)

                print("Analyze the reads which are mapped in one read and partially mapped in the other read.")
                pairedMapped1_Partial2(resSam)

                #print("Analyze the reads which are mapped in one read and unmapped in the other read.")
                #mapped1_unmapped2(resSam)

                print("Construct the final summary table.")
                countFlag()
                numberReads = total()
                percFlag(numberReads)

                # Remove the intermediate tables.
                outFile()

        print("Analyse finished. \n")

        readSummary()
                
    except getopt.error as err:
        # Output error, and return with an error code
        print (str(err))
        usage()
        sys.exit(2)

############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main(sys.argv[1:])
