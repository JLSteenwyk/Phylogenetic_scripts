#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo.BaseTree import Clade
from Bio import Phylo
import re
from pprint import pprint

def modTree(
    tree
    ):
    """
    removes comments 
    so that the tree can be visualized
    in strap
    
    Parameters
    ----------
    argv: tree
        newick tree file
    """

    filename=tree+".mcmc.tre"

    # read in tree
    tree = Phylo.read(tree, 'nexus')
    for i in tree.get_nonterminals():
        if i.comment is None:
            1
        elif i.comment is not None:
            i.comment = None
        #pprint(i)
        for ii in i.clades:
            if ii.comment is None:
                1
            elif ii.comment is not None:
                ii.comment = None
            #print(vars(ii))

    Phylo.write(tree, filename, 'newick')

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
            print("\nDivergence time estimates as outputted from MCMCtree are modified so that 95% CI intervals are removed.")
            print("The input tree is specified with the -t parameter. The tree can then be easily viewed using the strap R")
            print("package or other similar software.")
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
    modTree(
        tree
        )

if __name__ == '__main__':
    main(sys.argv[1:])