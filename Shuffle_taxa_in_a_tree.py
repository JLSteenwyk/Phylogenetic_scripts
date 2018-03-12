#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import Phylo
import random 

def shuffle_leaves(
    tree, seed 
    ):
    """
    shuffles taxa in a tree
    
    Parameters
    ----------
    argv: tree
        newick tree file
    argv: seed
        seed for pseudorandom shuffle
    """

    # create name for output tree file
    outTree = tree+".shuffled.tre"

    # read in tree
    tree = Phylo.read(tree, 'newick')

    # initialize list for taxa names
    taxa = []

    # loop through tree and save taxa names to taxa
    for term in tree.get_terminals():
        taxa.append(term.name)

    # shuffle taxa randomly using a seed
    random.Random(seed).shuffle(taxa)

    # loop through tree and replace taxa names
    CNT = 0
    for term in tree.get_terminals():
        term.name=term.name.replace(term.name, taxa[CNT])
        CNT+=1

    for node in tree.get_nonterminals():
        node.confidence = None

    Phylo.write(tree, outTree, 'newick')


def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree  = ''
    seed  = ''

    try:
        opts, args = getopt.getopt(argv, "ht:s:")
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
            print("\nShuffle taxa in a tree")
            sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-s':
            if int(arg):
                seed = arg
            else:
                # error message
                print("\n\nThe seed must be an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    shuffle_leaves(
        tree, seed
        )

if __name__ == '__main__':
    main(sys.argv[1:])