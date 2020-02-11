#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio import Phylo
import statistics as stat
import numpy as np

def print_BS(
    tree
    ):
    """
    Prints out bootstrap values from a phylogenetic tree in newick format

    Parameters
    ----------
    argv: tree
        newick tree file
    """

    # read in tree
    tree = Phylo.read(tree, 'newick')

    # initialize list to hold bootstrap values
    bs_vals = []

    # populate bs_vals with bootstrap values
    for terminal in tree.get_nonterminals():
        # only include if a bootstrap value is present
        if terminal.confidence != None:
            bs_vals.append(terminal.confidence)

    # calculate various statistics
    mean          = stat.mean(bs_vals)
    median        = stat.median(bs_vals)
    mini         = np.min(bs_vals)
    maxi           = np.max(bs_vals)
    twenty_fifth  = np.percentile(bs_vals, 25)
    seventy_fifth = np.percentile(bs_vals, 75)

    print("{}\t{}\t{}\t{}\t{}\t{}".format(mean, median, mini, maxi, twenty_fifth, seventy_fifth))

    
def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree  = ''

    try:
        opts, args = getopt.getopt(argv, "ht:")
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
            # general script explanation
            print("\nPrints out various statistics regarding the bootstrap support of a phylogenetic tree.")
            print("Statistics printed include the mean, median, standard deviation, variance, and the 25th")
            print("and 75th percentile. The output is printed in a tab delimited format following the same")
            print("order such as that")
            print("col 1: mean")
            print("col 2: median")
            print("col 3: min")
            print("col 4: max")
            print("col 5: 25th percentile")
            print("col 6: 75th percentile")
            sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    print_BS(
        tree
        )

if __name__ == '__main__':
    main(sys.argv[1:])