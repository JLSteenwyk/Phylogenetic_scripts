#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo

def report_result(
    RCV, treeness
    ):
    """
    calculates and prints out treeness/RCV

    Parameters
    ----------
    argv: RCV
        RCV value
    argv: treeness
        treeness value
    """

    # calculate treeness over RCV and save result as ToverR
    ToverR = treeness/RCV
    print(ToverR)

def RCV(
    alignment, treeness
    ):
    """
    calculate RCV

    Parameters
    ----------
    argv: alignment
        alignment fasta file
    argv: treeness
        treeness value
    """

    # initialize holder for RCV 
    RCV = ''

    # string to hold all sequences
    concatSeq = ''
    # initialize a counter for the number of sequences in the input fasta file
    numRecords = 0

    # read in records of a fasta file
    records = list(SeqIO.parse(alignment, "fasta"))
    
    # for each record join concatSeq string and sequence as well as keeping track 
    # of the number of records
    for record in records:
        concatSeq  += record.seq
        numRecords += 1

    # dictionary to hold the average occurence of each sequence letter
    averageD = {}
    # loop through the different sequences that appear in the fasta file
    # population dictionary with how many times that sequence appears
    for seq in set(concatSeq):
    	averageD[seq] = (concatSeq.count(seq)/numRecords)

    # determine the length of the alignment
    alignmentLen = 0
    alignmentLen = len(records[0].seq)

    # intiailize list to hold the RCV values per ith taxa 
    # that will later be summed
    indivRCVvalues = []

    # loop through records again and calculate RCV for 
    # each taxa and append to indivRCVvalues
    for record in records:
        # temp holds a temporary value of the numerator before appending
        # to numeratorRCVvalues and then is reassigned to 0 when it goes
        # through the loop again
        temp = 0
        # calculates the absolute value of the ith sequence letter minus the average
        for seqLetter in set(concatSeq):
            temp += abs(record.seq.count(seqLetter)-averageD[seqLetter])
        indivRCVvalues.append(temp/(numRecords*alignmentLen))

    # sum of all RCV values
    RCV = (sum(indivRCVvalues))

    report_result(
        RCV, treeness
        )

def calculate_treeness(
    alignment, tree
    ):
    """
    calculates treeness of a newick tree file
    
    Parameters
    ----------
    argv: alignment
        alignment fasta file
    argv: tree
        newick tree file
    """

    # initialize variable to calculate treeness
    treeness = 0

    # read in tree
    tree = Phylo.read(tree, 'newick')

    # initialize variables for terminal branch length
    interLen = float(0.0)
    # determine internal branch lengths
    for interal in tree.get_nonterminals():
        # only include if a branch length value is present
        if interal.branch_length != None:
            interLen += interal.branch_length

    # initialize variable for total branch length
    totalLen = float(0.0)
    # determine total branch length
    totalLen = tree.total_branch_length()

    treeness = (float(interLen/totalLen))

    RCV(alignment, treeness)

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    alignment = ''
    tree      = ''

    try:
        opts, args = getopt.getopt(argv, "hi:t:")
    except getopt.GetoptError:
        # error message
        print("Error\nFor help use -h argument\n")
        sys.exit(2)
    # if no arguments are used print a help message
    if len(opts) == 0:
        # error message
        print("\nNo arguments provided...")
        print("For help use -h argument\n")
        sys.exit(2)
    # test for arguments
    for opt, arg in opts:
        if opt == '-h':
            ## explanation
            print("\nCalculate treeness and relative nucleotide composition variability (RCV) as defined by")
            print("Phillips and Penny 2003, in the manuscript titled 'The root of the mammalian")
            print("tree inferred from whole mitochondrial genomes. Mol Phylogenet Evol. 28:171–185.'")
            print("https://www.ncbi.nlm.nih.gov/pubmed/12878457\n")
            print("The formula for RCV is as follows:")
            print("the sum of i=1^n of (|Ai - A*| + |Ti - T*| + |Ci - C*| + |Gi - G*|)/n•t")
            print("where, Ai, Ti, Ci, and Gi are the number of nucleotides for the ith taxon,")
            print("A*, T*, C*, and G* are the averages of each nucleotide across the n taxa")
            print("t is the number of sites and n is the number of taxa.\n")
            print("The formula to calculate treeness sum of internal branch lengths divided by total tree length.")
            print("\nThe output will a number that is the treeness/RCV value for a given alignment and tree.\n")
            ## options
            # alignment files list
            print("-i <alignment fasta file>")
            print("\tA multi-fasta sequence alignment file")
            # tree file
            print("-t <newick tree file>")
            print("\tNewick tree file associated with the alignment")
            print("")
            sys.exit()

        elif opt == '-i':
            if os.path.isfile(arg):
                alignment = arg
            else:
                # error message
                print("\n\nThe specified alignment list (-i) file does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to read_config parameter
    calculate_treeness(
        alignment, tree
        )

if __name__ == '__main__':
    main(sys.argv[1:])