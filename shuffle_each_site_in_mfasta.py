#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random

def shuffle(
    fasta, sequence_type, seed
    ):
    """
    Reads the alignment list and taxa list
    
    Parameters
    ----------
    argv: fasta
        input fasta file
    argv: sequence_type
        specifies sequence type
    argv: seed
    	pseudorandomizes sequences
    """

    # read in fasta file to alignment variable
    format    = "fasta"
    handle    = open(fasta)
    alignment = list(SeqIO.parse(handle, format))

    # initialize step and window size
    step   = 1
    window = 1

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
    length  = len(alignment[0].seq)

    # initialize output dictionary and populate with taxa names
    outDict = {}
    for indiv in alignment:
        if indiv.id in taxa_list:
            outDict[indiv.id] = ''

    # set seed
    random_numbers = []
    random.seed(seed)
    for x in range(length):
        random_numbers.append(random.randint(1,9999999999))

    # loop through aligned sequence using step size
    for i, randoNum in zip(range(0, (int(length)+1) - int(step), int(step)), random_numbers):
        # initialize variables to hold sequence at position
        positionSeq = []
        # loop through all individuals and save 
        # window being analyzed into a list
        for k, v in seqDict.items():
            positionSeq.append(v[i:i+int(window)])
        
        # shuffle the position of the nucleotides
        random.seed(randoNum)
        random.shuffle(positionSeq)
        
        # append the shuffled sequence to the value string
        for indiv, sequence in zip(taxa_list, positionSeq):
            outDict[indiv]+=sequence

    # print results to a fasta file format
    for indiv in taxa_list:
        print(">"+indiv+"\n"+outDict[indiv])

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
            print("\nShuffle the nucleotides or amino acids at each site in a fasta file.")
            print("Because this script shuffles, the relative number of each nucleotide or")
            print("amino acid is maintained after shuffling. This script is specifically")
            print("designed to facilitate randomization tests. For example, the tree length")
            print("test for clonality -- see figure 4 in http://jcm.asm.org/content/41/2/703.full")
            ## options
            # fasta file
            print("\n-i <input fasta file>:")
            print("\tSingle column file of the alignment files that will be concatenated.")
            # specify if nucleotide or protein files list
            print("\n-t <fasta files are either nucleotide or protein fastas>:")
            print("\tArgument can be either 'prot' or 'nucl' to specify if the")
            print("\talignments contain protein or nucleotide sequences.")
            # seed integer
            print("\n-s <seed integer>:")
            print("\tset seed for psuedorandomization of sequences\n")
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
        elif opt == '-s':
            if int(arg):
                seed = arg
            else:
                # error message
                print("\n\nThe specified input for (-s) should be an integer.\n")

    # pass to shuffle function
    shuffle(
        fasta, sequence_type, seed
        )

if __name__ == '__main__':
    main(sys.argv[1:])