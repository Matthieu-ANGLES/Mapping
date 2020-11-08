#!/usr/bin/python3
# -*- coding : utf-8 -*-

import sys,re

# Arguments
pyScriptName = sys.argv[0]
samFileName = sys.argv[1]

# Empty list
sam_header = []
flag = []
res = []

# Read sam file
with open(samFileName, "r") as sam_file:
    for line in sam_file:
        if line.startswith("@"):
            sam_header = line.strip()
        else:
            flag = line.split("\t")
            print(flag)
            res = str(flag[0]), int(flag[1])
            if res[0] == str(flag[0]):
                print(res[0] + " = " + str(flag[0]))
                print(res[0] + res[1] + int(flag[1])


