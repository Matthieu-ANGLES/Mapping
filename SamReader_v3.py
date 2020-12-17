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
    """
      Open a sam file and get the lines after the header.   
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
    """
    Docstring
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

def flagBinary (flag) :
    flagB = bin(int(flag)) 
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 

    if len(flagB) < 12: # Size ajustement to 12 (normalized size)
        add = 12 - len(flagB) 
        for t in range(add):
            flagB.insert(0,'0')

    return (flagB)

def flagBin2(sam_line): # VOIR POUR ECRITURE DANS SUMMARY OU AUTRE UTILISATION
    """
      Provides binary traduction of the Flag.
      Fasta outputs options for Read Unmapped or others informations.
    """
    with open("bin_Flag_table.txt", "w") as output_bin_flag:
        
        flag_counts = {"TotalReads":0}
        ReadName1 = "Je suis le read1 et j'en suis fier"
        
        for l in sam_line:
            flag_counts["TotalReads"] += 1
            line = l.split("\t")
            ReadName = str(line[0])

            if ReadName != ReadName1 :
                ReadName1 = ReadName    # il est le read1
                lineRead1 = l
                #print (lineRead1)
            else :                      # il est le read2
                lineRead2 = l  
                #print (lineRead2)

                line1 = lineRead1.split("\t")
                line2 = lineRead2.split("\t")
                flagB1 = flagBinary(line1[1])
                flagB2 = flagBinary(line2[1])
                mapq1 = int(line1[4])
                mapq2 = int(line2[4])
                
                # Voir si utile :
                #cigar1 = readCigar(line1[5])
                #cigar2 = readCigar(line2[5])
                #GC1 = percentGC(line1[9])
                #GC2 = percentGC(line2[9])


                # Cas où les deux reads sont mappés :
                if int(flagB1[-2]) == 1 and int(flagB2[-2]) == 1 :
                    if mapq1 >= 60 and mapq2 >= 60 : # Qualité >+ 60 (ajouter condition cigar aussi ?)
                        # Comptage reads mappés bonnes qualités +2
                        if "Reads_correctement_mappés_et_bonne_qualite" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] = 2
                        else:
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] += 2
                        # Comptage combanaisons mappés mauvaises qualités +1 
                        if "Paires_mappés" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] = 1
                        else:
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] += 1
                    else:
                        # Comptage reads mappés mauvaises qualités +2
                        if "Reads_correctement_mappés_mais_mauvaise_qualite" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] = 2
                        else:
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] += 2
                        # Comptage combinaisons mauvaises qualités +1
                        if "Paires_mappés" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] = 1
                        else:
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] += 1

                # Cas où les deux reads ne sont pas mappés :
                if int(flagB1[-3]) == 1 and int(flagB2[-3]) == 1 :
                        # Comptage reads non-mappés +2
                        if "Reads_non_mappés" not in flag_counts.keys():
                            flag_counts["Reads_non_mappés"] = 2
                        else:
                            flag_counts["Reads_non_mappés"] += 2
                        # Comptage combainaison non-mappés +1
                        if "Paires_non_mappés" not in flag_counts.keys():
                            flag_counts["Paires_non_mappés"] = 1
                        else:
                            flag_counts["Paires_non_mappés"] += 1

                # Cas où le read1 est mappé et le read2 est non-mappé :
                if int(flagB1[-2]) == 1 and int(flagB2[-3]) == 1 :
                    # Comptage reads non-mappés +1
                    if "Reads_non_mappés" not in flag_counts.keys():
                        flag_counts["Reads_non_mappés"] = 1
                    else:
                        flag_counts["Reads_non_mappés"] += 1

                    if mapq1 >= 60 : # Qualité >+ 60 (ajouter condition cigar aussi ?)
                        # Comptage reads mappés bonnes qualités +1                      
                        if "Reads_correctement_mappés_et_bonne_qualite" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] = 1
                        else:
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] += 1     
                        # Comptage combinaison Un Read mappés bonne qualité et Autre Read non mappés +1
                        if "UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_non_mappés" not in flag_counts.keys():
                            flag_counts["UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_non_mappés"] = 1
                        else:
                            flag_counts["UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_non_mappés"] += 1
                    else:
                        # Comptage reads mappés mauvaise qualités +1 
                        if "Reads_correctement_mappés_mais_mauvaise_qualite" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] = 1
                        else:
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] += 1
                        # Comptage combinaison UnRead mappés mauvaise qualité et Autre Read non mappés +1
                        if "UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_non_mappés" not in flag_counts.keys():
                            flag_counts["UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_non_mappés"] = 1
                        else:
                            flag_counts["UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_non_mappés"] += 1

                # Cas où le read1 est non-mappé et le read2 est mappé :
                if int(flagB1[-2]) == 1 and int(flagB2[-3]) == 1 :
                    # Comptage reads non-mappés +1
                    if "Reads_non_mappés" not in flag_counts.keys():
                        flag_counts["Reads_non_mappés"] = 1
                    else:
                        flag_counts["Reads_non_mappés"] += 1

                    if mapq2 >= 60 : # Qualité >+ 60 (ajouter condition cigar aussi ?)
                        # Comptage reads mappés bonnes qualités +1                      
                        if "Reads_correctement_mappés_et_bonne_qualite" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] = 1
                        else:
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] += 1     
                        # Comptage combinaison Un Read mappés bonne qualité et Autre Read non mappés +1
                        if "UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_non_mappés" not in flag_counts.keys():
                            flag_counts["UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_non_mappés"] = 1
                        else:
                            flag_counts["UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_non_mappés"] += 1
                    else:
                        # Comptage reads mappés mauvaise qualités +1 
                        if "Reads_correctement_mappés_mais_mauvaise_qualite" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] = 1
                        else:
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] += 1
                        # Comptage combinaison UnRead mappés mauvaise qualité et Autre Read non mappés +1
                        if "UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_non_mappés" not in flag_counts.keys():
                            flag_counts["UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_non_mappés"] = 1
                        else:
                            flag_counts["UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_non_mappés"] += 1

                # Cas où le read1 est mappé et le read2 est dans le mauvais sens :
                if int(flagB1[-2]) == 1 and int(flagB2[-5]) == 1 :
                    # Comptage reads mauvais sens +1
                    if "Reads_mauvais_sens" not in flag_counts.keys():
                        flag_counts["Reads_mauvais_sens"] = 1
                    else:
                        flag_counts["Reads_mauvais_sens"] += 1

                    if mapq1 >= 60 : # Qualité >+ 60 (ajouter condition cigar aussi ?)
                        # Comptage reads mappés bonnes qualités +1                      
                        if "Reads_correctement_mappés_et_bonne_qualite" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] = 1
                        else:
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] += 1     
                        # Comptage combinaison Un Read mappés bonne qualité et Autre Read dans le mauvais sens +1
                        if "UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_mauvais_sens" not in flag_counts.keys():
                            flag_counts["UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_mauvais_sens"] = 1
                        else:
                            flag_counts["UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_mauvais_sens"] += 1
                    else:
                        # Comptage reads mappés mauvaise qualités +1 
                        if "Reads_correctement_mappés_mais_mauvaise_qualite" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] = 1
                        else:
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] += 1
                        # Comptage combinaison UnRead mappés mauvaise qualité et Autre Read mauvais sens +1
                        if "UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_mauvais_sens" not in flag_counts.keys():
                            flag_counts["UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_mauvais_sens"] = 1
                        else:
                            flag_counts["UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_mauvais_sens"] += 1

                # Cas où le read1 est dans le mauvais sens et le read2 est mappé :
                if int(flagB1[-5]) == 1 and int(flagB2[-2]) == 1 :
                    # Comptage reads mauvais sens +1
                    if "Reads_mauvais_sens" not in flag_counts.keys():
                        flag_counts["Reads_mauvais_sens"] = 1
                    else:
                        flag_counts["Reads_mauvais_sens"] += 1

                    if mapq2 >= 60 : # Qualité >+ 60 (ajouter condition cigar aussi ?)
                        # Comptage reads mappés bonnes qualités +1                      
                        if "Reads_correctement_mappés_et_bonne_qualite" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] = 1
                        else:
                            flag_counts["Reads_correctement_mappés_et_bonne_qualite"] += 1     
                        # Comptage combinaison Un Read mappés bonne qualité et Autre Read dans le mauvais sens +1
                        if "UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_mauvais_sens" not in flag_counts.keys():
                            flag_counts["UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_mauvais_sens"] = 1
                        else:
                            flag_counts["UnRead_correctement_mappés_et_bonne_qualite_et_AutreRead_mauvais_sens"] += 1
                    else:
                        # Comptage reads mappés mauvaise qualités +1 
                        if "Reads_correctement_mappés_mais_mauvaise_qualite" not in flag_counts.keys():
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] = 1
                        else:
                            flag_counts["Reads_correctement_mappés_mais_mauvaise_qualite"] += 1
                        # Comptage combinaison UnRead mappés mauvaise qualité et Autre Read mauvais sens +1
                        if "UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_mauvais_sens" not in flag_counts.keys():
                            flag_counts["UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_mauvais_sens"] = 1
                        else:
                            flag_counts["UnRead_correctement_mappés_mais_mauvaise_qualite_et_AutreRead_mauvais_sens"] += 1

        #return flag_counts
        for k in flag_counts :
            output_bin_flag.write(k+" : "+str(flag_counts[k])+"\n")

# Version Benoit -----------------------------------------------------------------
def flagBin(sam_line, out): # VOIR POUR ECRITURE DANS SUMMARY OU AUTRE UTILISATION
    """
      Provides binary traduction of the Flag.
      Fasta outputs options for Read Unmapped or others informations.
    """
    #with open("bin_Flag_table.txt", "w") as output_bin_flag:
    flag_description = {}
    
    for l in sam_line:          
        line = l.split("\t") 
        flag = str(line[1])  # Line Flag retrieval
        
        te = bin(int(flag)) #### entre 0 et 4095 => mettre une securite !!!!!!
        te = te[2:] # Remove '0b' Example: '0b1001101' > '1001101'
        te = list(te) 

        if len(te) < 12: # Size ajustement to 12 (normalized size)
            add = 12 - len(te) 
            for t in range(add):
                te.insert(0,'0')

        # Flag traduction
        if 1 == int(te[-1]): # "Read paired"
            if "Read paired" not in flag_description.keys():
                flag_description["Read paired"] = 1
            else:
                flag_description["Read paired"] += 1

        if 1 == int(te[-2]): # "Read mapped in proper pair"
            if "Read mapped in proper pair" not in flag_description.keys():
                flag_description["Read mapped in proper pair"] = 1
            else:
                flag_description["Read mapped in proper pair"] += 1
                  
        if (1 == int(te[-3]) and out == "umo"): # "Read unmapped"
            #with open ("Reads_unmapped_only.txt","a+") as outputUMO:
            #    outputUMO.write(toStringOutput(l))
                
                if "Read unmapped" not in flag_description.keys():
                    flag_description["Read unmapped"] = 1
                else:
                    flag_description["Read unmapped"] += 1
                
        if 1 == int(te[-4]): # "Mate unmapped"
            if "Mate unmapped" not in flag_description.keys():
                flag_description["Mate unmapped"] = 1
            else:
                flag_description["Mate unmapped"] += 1
                
        if 1 == int(te[-5]): # "Read reverse strand"
            if "Read reverse strand" not in flag_description.keys():
                flag_description["Read reverse strand"] = 1
            else:
                flag_description["Read reverse strand"] += 1
                
        if 1 == int(te[-6]): # "Mate reverse strand"
            if "Mate reverse strand" not in flag_description.keys():
                flag_description["Mate reverse strand"] = 1
            else:
                flag_description["Mate reverse strand"] += 1
            
        if 1 == int(te[-7]): # "First in pair"
            if "First in pair" not in flag_description.keys():
                flag_description["First in pair"] = 1
            else:
                flag_description["First in pair"] += 1
                
        if 1 == int(te[-8]): # "Second in pair"
            if "Second in pair" not in flag_description.keys():
                flag_description["Second in pair"] = 1
            else:
                flag_description["Second in pair"] += 1

        if 1 == int(te[-9]): # "Not primary alignment"
            if "Not primary alignment" not in flag_description.keys():
                flag_description["Not primary alignment"] = 1
            else:
                flag_description["Not primary alignment"] += 1
                
        if 1 == int(te[-10]): # "Read fails platform/vendor quality checks"
            if "Read fails platform/vendor quality checks" not in flag_description.keys():
                flag_description["Read fails platform/vendor quality checks"] = 1
            else:
                flag_description["Read fails platform/vendor quality checks"] += 1
                
        if 1 == int(te[-11]): # "Read is PCR or optical duplicate"
            #with open ("Reads_is_PCR_or_optical_duplicate.txt","a+") as outputOD:
         #   outputOD.write(toStringOutput(l))
            if "Read is PCR or optical duplicate" not in flag_description.keys():
                flag_description["Read is PCR or optical duplicate"] = 1
            else:
                flag_description["Read is PCR or optical duplicate"] += 1
                
        if 1 == int(te[-12]): # "Supplementary alignment"
            if "Supplementary alignment" not in flag_output.keys():
                flag_description["Supplementary alignment"] = 1
            else:
                flag_description["Supplementary alignment"] += 1
                # with open ("Supplementary_alignment.txt","a+") as outputSA:
                #    outputOD.write(toStringOutput(l))

    for k in flag_description.keys():
        output_bin_flag.write(k + ":" + str(flag_description[k]) + "\n")
            
    return flag_description

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
            newline = ";".join([number2[0],number2[1],str(round(perc,4))])
            #print(newline)
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
    """
      Formula :  ( (G+C) / (A+T+G+C) * 100 )
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
    """
    Docstring
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

def toStringOutput (read_line):
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

def computeSummary (fileName):
    with open("summary_header.txt", "r") as F_Head, open("bin_Flag_table.txt", "r") as F_flag, open("Final_Cigar_table.txt", "r") as F_Cigar, open("Final_GC_table.txt", "r") as F_GC, open("summary.txt", "a+") as F_Sum :
        F_Sum.write("==================================================================\n"+
                    "  >  Summary "+fileName+"\n")
        for line in F_Head :
            l = line.rstrip("\n")
            F_Sum.write(line)
        F_Sum.write("\n")
        F_Sum.write("======================= PARSE INFORMATIONS =======================\n")
        for line in F_flag :
            l = line.rstrip("\n")
            F_Sum.write(line)
        F_Sum.write("\n")    
        for line in F_Cigar :
            l = line.rstrip("\n")
            F_Sum.write(line)
        F_Sum.write("\n")
        for line in F_GC :
            l = line.rstrip("\n")
            F_Sum.write(line)

def outFile(fileName): 
    """
    Output file name
    """

    os.remove("summary_header.txt")
    os.remove("parse_flag_table.txt")
    os.remove("count_flag_table.txt")
    os.remove("outpuTable_cigar.txt")
    os.remove("outpuTable_GC_percent.txt")
    #os.remove("Final_Flag_table.txt")
    os.remove("Final_Cigar_table.txt")
    os.remove("Final_GC_table.txt")
    os.rename("summary.txt", fileName + "_summary.txt")

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
                fileName = checkSamFormat(current_value)
                
                print("Read the file.")
                resSamHeader = openSamHeader(current_value)
                resSam = openSam(current_value)

                print("Parsing.")
                parseSamLine(resSam)
                #parseCigar(resSam)
                flagBin2(resSam)

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

            if current_argument in ("-o", "--output"):
                print("Ouput the file.")
                # voir pour options ?!!
                #outFile(fileName)  mis en dernier car affichage et renommage summary plus bas

            if current_argument in ("-f", "--fasta"):
                print("Compute the fasta file.")
                fastaOutput(resSam, current_value) # for commons Flags
                flagBin(resSam, current_value) # for Read Unmapped Only


        print("Analyse finished.")

        with open ("summary.txt") as summary :
            for line in summary :
                l = line.rstrip("\n")
                print (l)
                
    except getopt.error as err:
        # Output error, and return with an error code
        print (str(err))
        usage()
        sys.exit(2)

############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main(sys.argv[1:])
