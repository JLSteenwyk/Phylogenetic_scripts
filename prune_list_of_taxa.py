#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

def prune(
    tree, taxa
    ):
    """
    prune taxa from a newick tree
    
    Parameters
    ----------
    argv: tree
        newick tree file
    argv: taxa
        taxa to prune from tree
    """

    # create output file name
    filename = tree+".pruned"

    # read in tree
    tree = Phylo.read(tree, 'newick')

    # save taxa into a list
    taxaL  = []
    taxaL = [line.rstrip('\n') for line in open(taxa)]

    # prune taxa
    for taxon in taxaL:
        tree.prune(taxon)

    # write tree
    Phylo.write(tree, filename, 'newick')

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree = ''
    taxa = ''

    try:
        opts, args = getopt.getopt(argv, "hi:l:")
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
            print("\nPrunes list of a taxa from a newick tree file.")
            print("Arguments include -i <newick.tree> and -l <taxa.list>\n")
            print("Where the -i argument specifies a tree in newick format and")
            print("the -l option specifies a single column file with taxa to prune.")
            print("Output is a tree with the same name as the input tree with")
            print(".pruned appended to the end of the name.")
            sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-l':
            if os.path.isfile(arg):
                taxa = arg
            else:
                # error message
                print("\n\nThe specified list file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    prune(
        tree, taxa
        )

if __name__ == '__main__':
    main(sys.argv[1:])