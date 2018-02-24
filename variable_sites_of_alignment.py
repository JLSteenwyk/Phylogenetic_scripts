#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pprint import pprint

def variable_sites(
    fasta, sequence_type
    ):
    """
    Reads the alignment list and taxa list
    
    Parameters
    ----------
    argv: fasta
        input fasta file
    argv: sequence_type
        specifies sequence type
    """

    # read in fasta file to alignment variable
    format    = "fasta"
    handle    = open(fasta)
    alignment = list(SeqIO.parse(handle, format))

    # create a variable with all taxa names
    taxa_list = []
    for record in alignment:
        taxa_list.append(record.id)

    # initialize dictionary to hold sequences for each taxa
    seqDict = {}
    # populate seqDict with individual id (key) and sequence (value)
    for indiv in alignment:
        if indiv.id in taxa_list:
            seqDict[indiv.id] = indiv.seq
    
    # determine length of sequence
    length   = len(alignment[0].seq)
    # set step and window
    step     = 1
    window   = 1
    # intialize variable to hold number of variable sites
    varSites = 0

    # loop through fasta alignment file one position at a time
    for i in range(0, (int(length)+1) - int(step), int(step)):
        # list to hold sequence at that position
        positionSeq = ''
        # loop through individuals in sequence and append the sequence
        # at the position in the loop
        for k, v in seqDict.items():
            positionSeq += (v[i:i+int(window)])
        # extract only the sequence from the sequence object
        positionSeq = positionSeq._data
        # remove gaps (-?) and ambiguous sites (NX)
        positionSeq = re.sub('[-?NX]', '', positionSeq)
        # if there are multiple sequences in the position, add 1 to varSites
        if len(set(positionSeq)) > 1:
        	varSites += 1
    print("{}\t{}".format(varSites, varSites/length))

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    alignment_list  = ''
    taxa_list       = ''
    prefix          = ''

    try:
        opts, args = getopt.getopt(argv, "hi:t:s:")
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
            print("\nDetermine the number of variable sites in a alignment fasta file.")
            print("This script is designed to help facilitate creating phylogenomic matrix")
            print("subsets. The output has two columns where the first column is the number")
            print("of variable sites and the second column is the percentage.")
            ## options
            # fasta file
            print("\n-i <input fasta file>:")
            print("\tSingle column file of the alignment files that will be concatenated.")
            # specify if nucleotide or protein files list
            print("\n-t <fasta files are either nucleotide or protein fastas>:")
            print("\tArgument can be either 'prot' or 'nucl' to specify if the")
            print("\talignments contain protein or nucleotide sequences.\n")
            sys.exit()

        elif opt == '-i':
            if os.path.isfile(arg):
                fasta = arg
            else:
                # error message
                print("\n\nThe specified fasta file (-i) does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if arg in ('prot', 'nucl'):
                sequence_type = arg 
            else:
                # error message
                print("\n\nThe specified input for (-t) should be 'prot' or 'nucl'.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()

    # pass to variable_sites function
    variable_sites(
        fasta, sequence_type
        )

if __name__ == '__main__':
    main(sys.argv[1:])
