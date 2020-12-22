# Mapping Project
## Projet mapping HMIN113 BCD

NAME:

    SamReader

VERSION:

    0.1.1

DESCRIPTION:

    SamReader is a parser program for samtools files (extension: .sam).

    This software will return several files:

     - Final_Flag_table.txt
        
        This text file contains a table which summarize the whole flag combination (read 1 and read2). Their numbers and percentages are computed. 

     - 5 fasta files are returned by SamReader:

        * only_partially_mapped.fasta This file contains the reads (not paired) which are partially mapped.

        * only_unmapped.fasta This file contains the reads (not paired) which are unmapped.

        * read1_unmapped_read2_unmapped.fasta This file contains the both reads (Read 1 and Read 2) which are unmapped. 

        * read1_mapped_read2_partially.fasta
        This file contains one read (Read 1 or Read 2) which are mapped and the other one (Read 1 or Read 2) which are partially mapped.

        * read1_mapped_read2_unmapped.fasta This file contains one read (Read1 or Read 2) which are mapped and the other one (Read 1 or Read 2) which are unmapped.

    - A summary file is returned too. In this document, more information is available as mapper programm, version ...

    - This software analyze only Paired-End (PE) dataset. Contact us if you want a new version with the Single-End (Authors).

    More information about SAM file, please visit this website: https://samtools.github.io/hts-specs/SAMv1.pdf

REQUIREMENT:

    Python3 is necessary to launch the programm. If you don't have Python 3, you can download here:  https://www.python.org/downloads/

    This software work with Linux. We are not sure this programm works with windows.

INSTALLATION AND LAUNCH THE SOFTWARE:

    Download SamReader.py from our Github page (https://github.com/wanaga3166/Mapping/tree/main). 

    Move the SamReader.py in your tool directory.

    Open your terminal, go to your tool directory. 

    Use this command:

    python3 SamReader.py - h 

    Please see the help menu as showed above to read the full option list or check the next paragraph below.

OPTION LIST:

    -h or --help : help information
    -i or --input: input file (.sam)
    -o or --output: output name files

EXAMPLES:

    python3 SamReader.py -h # Launch the help
    python3 SamReader.py -i <file> # Launch SamReader to analyse a samtools file (.sam).
    python3 SamReader.py -i <file> -o <name>
    python3 SamReader.py --input <file>
    python3 SamReader.py --help
    python3 SamReader.py --input <file> --output <name>  

TROUBLESHOOTINGS:

    If you encounter one or several malfunctions when you execute this software, please contact us (see our e-mail addresses below).

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
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.

