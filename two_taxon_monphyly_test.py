#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

def determine_monophyly(
    tree, first, second
    ):
    """
    Determines if the two taxa are monophyletic
    
    Parameters
    ----------
    argv: tree
        newick tree file
    argv: first
        first taxa in tree
    argv: second
        second taxa in tree
    """

    # read in tree
    tree = Phylo.read(tree, 'newick')

    # read in taxa names
    clade      = [first, second]
    # initialize list to hold taxa names
    cladeCheck = []

    # loop through terminal branch
    for term in tree.get_terminals():
        # if outgroup taxa is present, append it to outPres
        if term.name in clade:
            cladeCheck.append(term.name)

    # check that length of cladeCheck length is 2 to ensure all
    # taxa are present
    if 2 == len(cladeCheck):
        1
    else:
        print("Both taxa are not in the tree\nExiting now...")
        sys.exit()

    for term in tree.common_ancestor(first, second):
        print(len(term.clades))
        print(vars(term))


    
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
        opts, args = getopt.getopt(argv, "ht:f:s:")
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
            print("\nDetermines if two species are monophyletic.")
            sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-f':
            if arg:
                first = arg
            else:
                # error message
                print("\n\nThe first taxa (-f) is incorrect.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-s':
            if arg:
                second = arg
            else:
                # error message
                print("\n\nThe second taxa (-s) is incorrect.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    determine_monophyly(
        tree, first, second
        )

if __name__ == '__main__':
    main(sys.argv[1:])
