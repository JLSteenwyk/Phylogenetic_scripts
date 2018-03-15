#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio import Phylo
from itertools import islice

def compare_trees(
    t1, t2
    ):
    """
    compares and determine which internodes conflict between
    a reference tree (t1) and a second tree (t2)
    
    Parameters
    ----------
    argv: t1
        reference tree
    argv: t2
        second tree
    """

    # read in tree
    t1 = Phylo.read(t1, 'newick')
    t2 = Phylo.read(t2, 'newick')

    for clade in t1.find_clades():
        print(clade)


def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree  = ''

    try:
        opts, args = getopt.getopt(argv, "hr:t:")
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
            print("\nCompares two trees in newick format and reports differences")
            sys.exit()
        elif opt == '-r':
            if os.path.isfile(arg):
                t1 = arg
            else:
                # error message
                print("\n\nThe specified reference tree file (-r) does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                t2 = arg
            else:
                # error message
                print("\n\nThe specified tree file (-t) does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    compare_trees(
        t1, t2
        )

if __name__ == '__main__':
    main(sys.argv[1:])