#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import SeqIO


def compare_lengths(
    alignmentM, alignmentT
    ):
    """
    compares the lengths of two aligned fasta files

    Parameters
    ----------
    argv: alignmentM
        alignmentM fasta file
    argv: alignmentT
        alignmentT fasta file
    """

    # initialize variables to hold sequence lengths
    Mlen     = ''
    Tlen     = ''

    # read in records of the fasta files
    recordsM = list(SeqIO.parse(alignmentM, "fasta"))
    recordsT = list(SeqIO.parse(alignmentT, "fasta"))

    # determine the sequence lengths
    Mlen     = len(recordsM[0].seq)
    Tlen     = len(recordsT[0].seq)

    # print the comparison
    print("{}\t{}\t{}".format(Tlen, Mlen, (Tlen/Mlen)*100))


def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    alignmentT = ''
    alignmentM = ''

    try:
        opts, args = getopt.getopt(argv, "hm:t:")
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
            print("\nCompare the alignment length between two alignments. This script was designed")
            print("for phylogenomics projects. More specifically, it is good to filter out genes who's")
            print("trimmed and aligned length is less than 50 percent the aligned gene's length.")
            print("The output is the length of the -t alignment, the length of the -m alignment,")
            print("followed by the percentage of the -t alignment divided by the -m alignment.")
            print("The logic of the -t and -m denomination is the -m represents a mafft alignment")
            print("and -t represents a trimAl alignment.\n")
            ## options
            # alignment files list
            print("\n-t <sequence alignment>")
            print("\tA multi-fasta sequence alignment file")
            print("\n-m <sequence alignment>")
            print("\tA multi-fasta sequence alignment file")
            sys.exit()

        elif opt == '-m':
            if os.path.isfile(arg):
                alignmentM = arg
            else:
                # error message
                print("\n\nThe specified alignment list (-m) file does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                alignmentT = arg
            else:
                # error message
                print("\n\nThe specified alignment list (-t) file does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()

    # pass to read_config parameter
    compare_lengths(
        alignmentM, alignmentT
        )

if __name__ == '__main__':
    main(sys.argv[1:])