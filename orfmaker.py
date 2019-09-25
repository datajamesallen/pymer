""" simple program to parse out the open reading frame from a fasta file """

import sys
import os.path

def makeorf(fastafile, start, end, outfile = "orf.fasta"):
    """ returns the sequence of only the Coding Domain Sequence / Open Reading Frame
        you can find this information on the NCBI NM_ ref seq under CDS"""
    with open(fastafile, 'r') as fasta:
        head = fasta.readline()
        seq = fasta.read()

        # check to make sure that the first line is the fasta comment,
        # if not then add the first line into the sequence, as it is likely DNA bases
        if not head[0] == ">":
            seq = head + seq

        # get rid of whitespace
        seq = seq.replace('\n','')
        orf = seq[start-1:end]
    fasta.close()
    with open(outfile, 'w') as orffile:
        orffile.write(orf)
    orffile.close()
    return(None)

if __name__ == "__main__":
    args = sys.argv
    #print(args)
    print(args[2])
    print(args[3])
    try:
        START = int(args[2])
        END = int(args[3])
    except:
        sys.exit("Requires argument 2 and 3 to be integers start and end of CDS/ORF\n example: python3 orfmaker.py GRIN2A.txt 311 4705.")
    if not os.path.isfile(args[1]):
        sys.exit("Argument 1 must be a filename\nexample: python3 orfmaker.py GRIN2A.txt 311 4705")
    # create default file name if argument 4 is not set (orf.fasta)
    try:
        OUTFILE = args[4]
        makeorf(args[1], start = START, end = END, outfile = OUTFILE)
    except:
        OUTFILE = "orf.fasta"
        makeorf(args[1], start = START, end = END)
    print("orf file \"" + OUTFILE + "\" written")


