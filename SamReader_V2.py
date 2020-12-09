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
         python SamReader.py -i <file> -f <option> # explain
    AUTHORS:
         Benoit ALIAGA (aliaga.benoit@gmail.com)
         Matthieu ANGLES (matthieu.angles@hotmail.fr)
    MORE INFORMATIONS:
         https://github.com/wanaga3166/Mapping
    
    ''')
    sys.exit(-1)

def checkSamFormat (argv): 
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
            if line.startswith("@"):             # Check header motif  
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
   

def openSamHeader(argv):
    """
    Read a sam file to give header
    """
    sam_header = []
    with open(argv, "r") as sam_file, open("summary_header.txt", "w") as summary_header :
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
            '''
            if line.startswith("@RG"):   on fait aussi ?
            if line.startswith("@CO"):
            if line.startswith("@HD"):
            '''

        print ("======================= HEADER INFORMATION =======================")
        print ("Reference sequence name :",SN)
        print ("Reference sequence lenght :",LN)
        print ("Programme record identifier :",ID)
        print ("Programme name :",PN)
        print ("Programme version :",VN)
        print ("Command Line :",CL)
        print ("==================================================================")
        
        summary_header.write("======================= HEADER INFORMATION ======================="+"\n")    
        summary_header.write("Reference sequence name : "+SN+"\n")
        summary_header.write("Reference sequence lenght : "+LN+"\n")
        summary_header.write("Programme record identifier : "+ID+"\n") 
        summary_header.write("Programme name : "+PN+"\n")
        summary_header.write("Programme version : "+VN+"\n")
        summary_header.write("Command Line : "+CL+"\n")
        summary_header.write("=================================================================="+"\n")
        
    return sam_header

def openSam(argv):
    """
    Read a sam file to give reads
    """

    sam_header = []
    sam_line = []
    
    with open(argv, "r") as sam_file, open("summary.txt", "w") as summary :
        for line in sam_file:
            if line.startswith("@"):
               pass
            else:
                sam_line.append(line)
        return sam_line
        
def outFile(argv):
    """
    Output file name
    """
    os.remove("parse_flag_table.txt")
    os.remove("count_flag_table.txt")
    os.rename("Final_Flag_table.txt", argv)

    
def parseFlagCigar(sam_line):
    """
    Docstring
    """
    cpt = 1 # Initialisation d'un compteur
    with open("parse_flag_table.txt", "w") as output_sam_flag:
        for line in sam_line:
            flag = line.split("\t") # Fractionnement de la line suivant \t
            if cpt == 1:
                add_line = [flag[0], flag[1], flag[5]] # ajout de flag1 et cigar1
                #print(add_line)
                cpt +=1
            elif cpt == 2:
                add_line.append(flag[1])
                add_line.append(flag[5]) # ajout de flag2 et de cigar 2
                new_line = "\t".join([add_line[0], add_line[1], add_line[2], add_line[3], add_line[4]])
                #print(new_line)
                
                output_sam_flag.write(add_line[0] + ";" + add_line[1] + ";" + add_line[2] + ";" + add_line[3] + ";" + add_line[4] + "\n")
                #output_sam_flag.write(
                
                cpt = 1

def countFlag():
    """
    Docstring
    """
    with open("parse_flag_table.txt", "r") as table_flag, open("count_flag_table.txt", "w") as output_flag:
        dico = {}
        for line in table_flag:
            flag = line.rstrip("\n").rsplit(";")
            inter = flag[1] + "-" + flag[3]
            if inter not in dico.keys():
                dico[inter] = 1
            else:
                dico[inter] += 1
            #print(dico)
        for key in dico.keys():
            output_flag.write(key + ";" + str(dico[key]) + "\n")

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

def parseCigar(sam_line):
    cpt = 1
    with open("outpuTable_cigar.txt", "w") as output_sam_cigar:
        for line in sam_line:
            #print(line)
            cigar = line.split("\t") # Split args line
            if cpt == 1:
                add_line = [cigar[0], cigar[5]] # Creation line with Read and cigar R1
                cpt +=1
            elif cpt == 2:
                add_line.append(cigar[5]) # Add cigar R2 into line
                c1 = readCigar(add_line[1])
                #print(c1)
                c2 = readCigar(add_line[2])
                #print(c2)
                output_sam_cigar.write(add_line[0] + ";" + add_line[1] + ";" + percentMutation(c1) + add_line[2] + ";" + percentMutation(c2) + "\n")
                cpt = 1
                

def readCigar(cigar): # Make a directory with lists key, value
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

def percentMutation (dico):
    """
    Functions readCigar : Give percentage for each mutation in readCigar dictionnary
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

def percentGC (seq) :
    countGC = 0
    total = 0
    # Count nucleotides :
    for n in seq :
        total += 1
        #print (n)
        if (n == 'G') or (n == 'C') :
            countGC += 1
    # Calcul percentage :
    percentGC = (countGC * 100) / total
    return(percentGC)


def parseMAPQ():
    """
    Docstring
    """
    pass

def toStringOutput (sam_line):
    """
    Write function for output
    """
    line = sam_line.split("\t")
    qname = str(line[0])
    flag = str(line[1])
    seq = str(line[9])
    return ("> " + qname + " / flag:"+flag+"\n"+seq+"\n"+"\n")

def flagBin (sam_line, toto): # VOIR POUR ECRITURE DANS SUMMARY OU AUTRE UTILISATION
    """
    Provides binary translation of the Flag
    Fasta outputs options for Read Unmapped or others informations
    """
    for l in sam_line:          
        line = l.split("\t") 
        flag = str(line[1])  # Line Flag retrieval

        te = bin(int(flag)) #### entre 0 et 4095 => mettre une securite !!!!!!
        te = te[2:] # Cleaning regex '0b'  >'0b1001101'
        te = list(te) 

        if len(te) < 12: # size ajustement to 12 (normalyzed size)
            add = 12 - len(te) 
            for t in range(add):
                te.insert(0,'0')

        # Flag traduction
        if 1 == int(te[-1]): # "Read paired"
            pass

        if 1 == int(te[-2]): # "Read mapped in proper pair"
            pass
            # faire un comptage ici ??

        if (1 == int(te[-3]) and toto == "umo"): # "Read unmapped"
            with open ("Reads_unmapped_only.txt","a+") as outputUMO:
                outputUMO.write(toStringOutput(l))
                # faire un comptage ici ??
                
        if 1 == int(te[-4]): # "Mate unmapped"
            pass
                
        if 1 == int(te[-5]): # "Read reverse strand"
            pass
                
        if 1 == int(te[-6]): # "Mate reverse strand"
            pass
                
        if 1 == int(te[-7]): # "First in pair"
            pass
                
        if 1 == int(te[-8]): # "Second in pair"
            pass

        if 1 == int(te[-9]): # "Not primary alignment"
            pass
                
        if 1 == int(te[-10]): # "Read fails platform/vendor quality checks"
            pass
                
        if 1 == int(te[-11]): # "Read is PCR or optical duplicate"
            with open ("Reads_is_PCR_or_optical_duplicate.txt","a+") as outputOD:
                outputOD.write(toStringOutput(l))

        if 1 == int(te[-12]): # "Supplementary alignment"
            with open ("Supplementary_alignment.txt","a+") as outputSA:
                outputOD.write(toStringOutput(l))
    
    #return(?)


def fastaOutput(sam_line, toto):
    """
    Fasta outputs options for Commons Flags
    """

    One_of_the_reads_is_unmapped = ('73','133','89','121','165','181','101','117','153','185','69','137')
    Both_reads_are_unmapped = ('77','141')
    Mapped_within_the_insert_size_and_in_correct_orientation = ('99','147','83','163')
    Mapped_within_the_insert_size_but_in_wrong_orientation = ('67','131','115','179')
    Mapped_uniquely_but_with_wrong_insert_size = ('81','161','97','145','65','129','113','177')
    test1 = ('83','163','77','141')

#### voir pour expatrier les listes ci-dessus en lecture de flag binaire ####

    for l in sam_line:
        line = l.split("\t")
        flag = str(line[1])
        
        if (flag in One_of_the_reads_is_unmapped and toto == "oum"):                            
            with open("OneUnMapped.fasta", "a+") as Output1:    
                Output1.write(toStringOutput(l))

        if (flag in Both_reads_are_unmapped and toto == "bum"):
            with open("BothUnMapped.fasta", "a+") as Output2:
                Output2.write(toStringOutput(l))

        if (flag in Mapped_within_the_insert_size_and_in_correct_orientation and toto == "co"):
            with open("CorrectOrientation.fasta", "a+") as Output3:
                Output3.write(toStringOutput(l))

        if (flag in Mapped_within_the_insert_size_but_in_wrong_orientation and toto == "wo"):
            with open("WrongOrientation.fasta", "a+") as Output4:
               Output4.write(toStringOutput(l))

        if (flag in Mapped_uniquely_but_with_wrong_insert_size and toto == "wi"):
            with open("WrongInsertSize.fasta", "a+") as Output5:
                Output5.write(toStringOutput(l))

        if (flag in test1 and toto == "test"):
            with open("test.fasta", "a") as Output6:
                Output6.write(toStringOutput(l))

def main(argv):
    """
    Docstring
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
                resSamHeader = openSamHeader(current_value)
                resSam = openSam(current_value)

                print("Parse the flag.")
                parseFlagCigar(resSam)
                parseCigar(resSam)

                print("Count the number of read for each flag combinations.")
                countFlag()

                print("Calculate the percentage for each flag combination.")
                numberReads = total()
                percFlag(numberReads)
                
            if current_argument in ("-o", "--output"):
                print("Ouput the file.")
                outFile(current_value)

            if current_argument in ("-f", "--fasta"):
                print("Output the fasta file.")
                fastaOutput(resSam, current_value) # for commons Flags
                flagBin(resSam, current_value) # for Read Unmapped Only

        print("Analyse finished.")
                
    except getopt.error as err:
        # Output error, and return with an error code
        print (str(err))
        usage()
        sys.exit(2)
        
if __name__ == "__main__":
    main(sys.argv[1:])