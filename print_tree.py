#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import Phylo

def print_tree(
    tree
    ):
    """
    prints newick tree in ASCII format
    
    Parameters
    ----------
    argv: tree
        newick tree file
    """

    # read in tree
    tree = Phylo.read(tree, 'newick')

    Phylo.draw_ascii(tree)

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree  = ''

    try:
        opts, args = getopt.getopt(argv, "hi:")
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
            print("\nPrint an ASCI version of a phylogenetic tree in newick format.")
            print("Only argument is -i <newick.tree>.\n")
            sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    print_tree(
        tree
        )

if __name__ == '__main__':
    main(sys.argv[1:])