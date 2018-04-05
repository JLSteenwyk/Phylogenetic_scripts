#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo.BaseTree import Clade
from Bio import Phylo

def modTree(
    tree
    ):
    """
    creates numerical identifiers for internodes
    
    Parameters
    ----------
    argv: tree
        newick tree file
    """

    # create file name 
    filename=tree+".topology.tree"

    # read in tree
    tree = Phylo.read(tree, 'newick')
    # set branch lengths to None
    for i in tree.get_nonterminals():
        i.branch_length=None
    for i in tree.get_terminals():
        i.branch_length=None

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
            print("\nRemoves distances from a tree in newick format. This can be used if you are only interested")
            print("in the topology of the tree or want to reestimate branch lengths after making a consensus tree.")
            print("The input newick tree file is specified with the -t parameter.")
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