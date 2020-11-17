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

import sys, re, getopt

def usage():
    sys.stderr.write('''
    
    NAME:
         SamReader

    VERSION:
         0.0.1

    FUNCTION: 
         SamReader is a parser program for samtool files (.sam). 

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
    with open(argv, "r") as sam_file:
        for line in sam_file:
            if line.startswith("@"):
                sam_header = line.strip()
                print(sam_header)
            else:
                print(line)

def parseFlag():
    """
    Docstring
    """
    pass

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
            elif current_argument in ("-i", "--input"):
                print("Read the file")
                openSam(current_value)
            elif current_argument in ("-o", "--output"):
                print("Ouput the file")
    except getopt.error as err:
        # Output error, and return with an error code
        print (str(err))
        usage()
        sys.exit(2)
        
if __name__ == "__main__":
    main(sys.argv[1:])
