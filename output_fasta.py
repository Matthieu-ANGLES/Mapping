import sys, re, os


# Arguments
pyScriptName = sys.argv[0]
samFileName = sys.argv[1]
fileName, fileExtention = os.path.splitext(samFileName)
samOutputName = fileName.split('/')

#samOutputName = samOutputName[-1]+".fasta"
samOutputName1 = samOutputName[-1]+"_One_of_the_reads_is_unmapped.fasta"
samOutputName2 = samOutputName[-1]+"_Both_reads_are_unmapped.fasta"
samOutputName3 = samOutputName[-1]+"_Mapped_within_the_insert_size_and_in_correct_orientation.fasta"
samOutputName4 = samOutputName[-1]+"_Mapped_within_the_insert_size_but_in_wrong_orientation.fasta"
samOutputName5 = samOutputName[-1]+"_Mapped_uniquely_but_with_wrong_insert_size.fasta"


test1 = ('83','163')
test2 = ('77','141')
One_of_the_reads_is_unmapped = ('73','133','89','121','165','181','101','117','153','185','69','137')
Both_reads_are_unmapped = ('77','141')
Mapped_within_the_insert_size_and_in_correct_orientation = ('99','147','83','163')
Mapped_within_the_insert_size_but_in_wrong_orientation = ('67','131','115','179')
Mapped_uniquely_but_with_wrong_insert_size = ('81','161','97','145','65','129','113','177')



with open(samFileName, "r") as sam_file, open(samOutputName1, "w") as Output1, open(samOutputName2, "w") as Output2, open(samOutputName3, "w") as Output3, open(samOutputName4, "w") as Output4, open(samOutputName5, "w") as Output5:                   
    for line in sam_file:
        if line.startswith("@"): # Header motif
            sam_header = line.strip()
            print(sam_header)
        else:
            

            # PARTIE UTILSE POUR FONCTION


            ligne = line.split("\t")
            qname = str(ligne[0])
            flag = str(ligne[1])
            seq = str(ligne[9])
            #print (qname, flag, seq)
            if (flag in One_of_the_reads_is_unmapped):
                Output1.write("> "+qname+" / flag:"+flag+"\n"+seq+"\n"+"\n")
            if (flag in Both_reads_are_unmapped):
                Output2.write("> "+qname+" / flag:"+flag+"\n"+seq+"\n"+"\n")
            if (flag in Mapped_within_the_insert_size_and_in_correct_orientation):
                Output3.write("> "+qname+" / flag:"+flag+"\n"+seq+"\n"+"\n")
            if (flag in Mapped_within_the_insert_size_but_in_wrong_orientation):
                Output4.write("> "+qname+" / flag:"+flag+"\n"+seq+"\n"+"\n")
            if (flag in Mapped_uniquely_but_with_wrong_insert_size):
                Output5.write("> "+qname+" / flag:"+flag+"\n"+seq+"\n"+"\n")
            
            if (flag in test1):
                Output1.write("> "+qname+" / flag:"+flag+"\n"+seq+"\n"+"\n")
            if (flag in test2):
                Output1.write("> "+qname+" / flag:"+flag+"\n"+seq+"\n"+"\n")





## samOutputName.write(add_line[0] 