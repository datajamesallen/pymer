"""
QCL Primer Designing tool 

James Allen 2019

Creates mutagenic primers for the Quik Change Lightning mutagenic PCR protocol.

Primer Specifications - from QuikChange Lightning Protocol
---------------------

Summary:

1.  Both mutagenic primers must contain the desired mutation and anneal to the same sequence on opposite strands of the plasmid
2.  Pimers should be between 25 and 45 bases in length. 
    Longer than 45 bases may be used but using longer pimers increases the likelihood of secondary structure formation, 
    which may affect the efficiency of the mutagenesis reaction.
3.  The following formula is commonly used for estimating the Tm of primers:
        Tm = 81.5 + 0.41(%GC) - (675/N) - %mismatch
            % mismatch -> 1 for SNP
            N          -> length of the primer
            % GC       -> percentage of G or C in sequence vs A or T
    For calculating Tm for primers intended to introduce insertions or deletions, use this modified version
        Tm = 81.5 + 0.41(%GC) - (675/N)
            where N does not include the bases which are being inserted or deleted.
4. The desired mutation (deletion or insertion) should be in the middle of the primer
     with ~10-15 bases of correct sequence on both sides
5. The primers optimally should have a minimum GC content of 40% and should terminate in one or more C or G bases
6. Primers need not be 5' phosphorylated

"""

def reverse_complement(seq):
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

class primer:
    """
    primer class
    """
    def __init__(self, sequence, transcript, position, bp_original,
                 bp_change, len_side_1, len_side_2):
        self.sequence = sequence
        self.sequence_r = reverse_complement(sequence)
        self.transcript = transcript
        self.position = position
        self.bp_original = bp_original
        self.bp_change = bp_change
        self.len_side_1 = len_side_1
        self.len_side_2 = len_side_2
        self.GC = False
        self.Tm = False
        self.GCend = False

    def __str__(self):
        printvar = "Sequence (f)  : " + self.sequence + "\n" \
                 + "Sequence (r)  : " + self.sequence_r + "\n" \
                 + "Base change   : c." + str(self.position) + self.bp_original + ">" + self.bp_change + "\n" \
                 + "GC content    : " + str(self.GC) + " %" + "\n" \
                 + "Melting temp  : " + str(self.Tm) + " degrees C" + "\n" \
                 + "Length base   : " + str(len(self.sequence)) + "\n"
        return printvar

    def update_rev(self):
        self.sequence_r = reverse_complement(self.sequence)
        return(self)

    def check_GCend(self):
        start = self.sequence[0]
        end = self.sequence[-1]
        if start not in ["G","C"] or end not in ["G","C"]:
            self.GCend = False
        else:
            self.GCend = True
        return self
    def inc_side_1(self):
        self.len_side_1 += 1
        additionbase_side_1 = self.transcript[self.position-1-self.len_side_1]
        self.sequence = additionbase_side_1 + self.sequence
        return self

    def inc_side_2(self):
        self.len_side_2 += 1
        additionbase_side_2 = self.transcript[self.position-1+self.len_side_2-1]
        self.sequence = self.sequence + additionbase_side_2
        return self

    def dec_side_1(self):
        self.len_side_1 -= 1
        self.sequence = self.sequence[1:]
        return self

    def inc(self):
        self.len_side_1 += 1
        self.len_side_2 += 1
        additionbase_side_1 = self.transcript[self.position-1-self.len_side_1]
        additionbase_side_2 = self.transcript[self.position-1+self.len_side_2-1]
        self.sequence = additionbase_side_1 + self.sequence
        self.sequence = self.sequence + additionbase_side_2
        return self

    def dec(self):
        self.len_side_1 -= 1
        self.len_side_2 -= 1
        self.sequence = self.sequence[1:-1]
        return self

    def check_GC(self):
        gc_count = 0
        at_count = 0
        for item in self.sequence:
            if item in ["A","T"]:
                at_count += 1
            if item in ["G","C"]:
                gc_count += 1
        self.GC = (gc_count/(at_count + gc_count)) * 100
        return self

    def check_Tm(self):
        """
        #    Calculate the Melting Temperature as:
        #    Tm = 81.5 + 0.41(%GC) - (675/N) - %mismatch
        #    % mismatch -> 1 for SNP
        #    N          -> length of the primer
        #    % GC       -> percentage of G or C in sequence vs A or T
        """
        self.check_GC()
        N = len(self.sequence)
        # percent missmatched (should be 1 most cases, test others later)
        per_mm = len(self.bp_change) 
        self.Tm = 81.5 + (0.41 * round(self.GC)) - (675/N) - per_mm
        return self

def makeprimer(transcript, result):
    """
    Creates and returns a mutagenic primer object given supplied args
    where transcript is the base template sequence string, and
    result = [{position (int)}, {original base pair (1 char string)}, {mutated base pair (1 char string)}]
    """
    # 1. add bases on each side until the primer is 23 bases in length
    # result from getinput() / validate_cDNA
    position = result[0]
    bp_original = result[1]
    bp_change = result[2]
    # create the initial primer 23 bases long
    len_side_1 = 11 # how much transcript is on each side
    len_side_2 = 11
    pm_start = position-1-11
    if pm_start < 0:
        sys.exit("transcript sequence too short to make a primer there!")
    pm_end = position+11
    if pm_end > len(transcript):
        sys.exit("transcript sequence too short to make a primer there!")
    # intial primer
    primer_initial_seq = transcript[pm_start:position-1] + bp_change + transcript[position:pm_end]
    print(primer_initial_seq)

    # this is not the best primer, but it is for now...
    primer_best = primer(sequence = primer_initial_seq, transcript = transcript,
                         position = position, bp_original = bp_original,
                         bp_change = bp_change, len_side_1 = len_side_1, len_side_2 = len_side_2)
    #while not primer_best.GCend:
        #primer_best.inc()
    #    primer_best.check_GCend()
    #    primer_best.inc()
    #primer_best.check_GC
    while True:
        # update data on the current primer:
        primer_best.check_GCend()
        primer_best.check_GC()
        primer_best.check_Tm()
        primer_best.update_rev()
        print(primer_best)
        # check all required conditions, if they are met, then break
        # the GC content must be over 40%, melting temperature above 78 degrees C
        # and it must terminate in G or C bases
        if all([primer_best.GC >= 40, 
                   primer_best.GCend,
                   primer_best.Tm >= 78]):
            break
        # primers also should not be more than 45 base pairs
        # so here is another way to exit the loop
        if len(primer_best.sequence) > 45:
            break
        print(primer_best.len_side_1)
        print(primer_best.len_side_2)
        if primer_best.len_side_1 > primer_best.len_side_2:
            print("extend right")
            primer_best.dec_side_1()
            primer_best.inc_side_2()
        elif primer_best.len_side_1 < primer_best.len_side_2:
            print("extend left")
            primer_best.inc_side_1()
        elif primer_best.len_side_1 == primer_best.len_side_2:
            print("extend left")
            primer_best.inc_side_1()
    print("Final primer: ", primer_best)
    return(primer_best)

# 3. Evalutate GC content. Each time save the information for the primer that had the maximum GC content
# 4. If GC content is less than 40%, add bases until G or C is reached, repeat until GC >= 40%, or the primer is 45 bases long.
# 5. If 45 bases is reached, then stop. Take the maxiumum GC content reached.

###--------------------------###
### Interface and User Input ###
###--------------------------###

import sys

from tkinter import Tk
from tkinter.filedialog import askopenfilename
from re import search

def gettranscript():
    """ opens a file dialoge and returns a fasta format transcript """
    print("Please specify a fasta file in the file dialogue")
    Tk().withdraw() # dont want the full GUI
    fastafile = askopenfilename(initialdir = ".",
                                    title = "Select a valid fasta file")
    Tk().update()
    with open(fastafile, 'r') as fasta:
        head = fasta.readline()
        seq = fasta.read()
    fasta.close()

    # check to make sure that the first line is the fasta comment,
    # if not then add the first line into the sequence, as it is likely DNA bases
    if not head[0] == ">":
        seq = head + seq

    # get rid of whitespace
    seq = seq.replace('\n','')
    seqres = search("[ATGC]+", seq)
    if not seqres:
        print("Invlaid fasta file, please fix the file and try again")
        sys.exit()
    else:
        return(seqres.group())

def getinput(transcript):
    """
    get the user input for the cDNA change specification 

    returns the validated result from validate_cDNA
    """
    while True:
        print("Please enter the cDNA change you would like to make")
        print("Example: c.1345A>T")
        cDNA_change = input()
        result = validate_cDNA(cDNA_change, transcript)
        if not result:
            print("Invalid cDNA change entered")
            continue
        else:
            print(result)
            return(result)


def validate_cDNA(cDNA_change, transcript):
    #clean_change = re.sub(r'[c.]', '', cDNA_change)
    position = search(r'[0-9]+', cDNA_change)
    change = search(r'[ATGC]>[ATGC]', cDNA_change)
    if not change or not position:
        return(False) # quit with error, see getinput()
    position = int(position.group())
    if position > len(transcript):
        return(False)
    change = change.group()
    bp_original = change[0]
    bp_change = change[-1]
    bp_original_on_transcript = transcript[position-1]
    if bp_original != bp_original_on_transcript:
        print("The supplied cDNA change does not match the supplied transcript")
        return(False)
    else:
        result = [position, bp_original, bp_change]
        return(result)

def checkprimers():
    # open the file
    # read and parse out those that have Gene = GRIN2A
    # for those entries, feed the cDNA change to validate_cDNA(), then makeprimer()
    transcript = open("ORF/orf-GRIN2A.fasta",'r').readlines()
    with open("QC_primers_all.csv") as qcfile:
        result = qcfile.readlines()
        num = 0
        counter = 0
        for row in result:
            row = row.rstrip()
            row = row.split(",")
            print(row)
            if row[0] == "GRIN2A" and row[6] == "forward":
                valid = validate_cDNA(row[1], transcript)
                primer = makeprimer(transcript, result)
                if primer.sequence == row[3]:
                    num += 1
                counter += 1
        print(num/counter)
    return None

def main():
    transcript = gettranscript()
    #print(transcript)
    print("transcript length:", len(transcript))
    result = getinput(transcript)
    makeprimer(transcript, result)

#def checker():


if __name__ == "__main__":
    main()
