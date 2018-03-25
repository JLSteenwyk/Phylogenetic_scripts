#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio import Phylo
import re
import numpy

def print_IC(
    tree
    ):
    """
    Prints out internode certainty (IC) values from a phylogenetic tree in newick format

    Parameters
    ----------
    argv: tree
        newick tree file
    """

    # read in tree
    tree = Phylo.read(tree, 'newick')

    # populate loop through internal nodes
    for terminal in tree.get_nonterminals():
        # only include if a IC value is present
        if terminal.comment != None:
            # print IC values
            print(float(re.sub(r',[^,]*$', '', terminal.comment)))
    
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
            print("\nPulls out the internode certainty (IC) values of a phylogenetic tree outputed from RAxML.")
            print("Specify the tree with the -n parameter")
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
    print_IC(
        tree
        )

if __name__ == '__main__':
    main(sys.argv[1:])