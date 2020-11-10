#!/usr/bin/python3
# -*- coding : utf-8 -*-

import sys,re

# Arguments
pyScriptName = sys.argv[0]
samFileName = sys.argv[1]
samOutputName = sys.argv[2]

# Empty list
res = []
res2 = []

# Read sam file
with open(samFileName, "r") as sam_file, open(samOutputName, "w") as output_sam_flag:
    for line in sam_file:
        if line.startswith("@"):
            sam_header = line.strip()
            print(sam_header)
        else:
            flag = line.split("\t")
            res = [str(flag[0]), int(flag[1])]
            print(res)
            print(type(res))
            if res[0] == str(flag[0]):
                print(res[0] , " = " , str(flag[0]))
                res2.append(res[0])
                res2.append(res[1])
                res2.append(int(flag[1]))
                print(res2)
                output_sam_flag.write(res2[0] + ";" + str(res2[1]) + ";" + str(res2[2]) + "\n")
                res = []
                

#print(res2)
