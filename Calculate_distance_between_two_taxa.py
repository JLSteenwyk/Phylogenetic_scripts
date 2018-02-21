#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio.Phylo.BaseTree import TreeMixin
from Bio import Phylo

def readTree(
    inputTree, target1, target2
    ):
    """
    Reads input tree

    Parameters
    ----------
    argv: inputTree
        reads the input tree file
    argv: target1
        target1 taxa 
    argv: target2
        target2 taxa
    """
    # initialize variables to store tree, leaf1, leaf2
    tree  = ''
    leaf1 = ''
    leaf2 = ''
    # read tree
    tree  = Phylo.read(inputTree, 'newick')
    
    # determine is leaves exists
    leaf1 = TreeMixin.find_any(tree, target1)
    if (bool(leaf1)) == False:
        #print(target1, "not on tree\nExiting...")
        sys.exit()
    leaf2 = TreeMixin.find_any(tree, target2)
    if (bool(leaf2)) == False:
        #print(target2, "not on tree\nExiting...")
        sys.exit()

    # determine leaf1 and leaf2 distance
    distance = ''
    distance = TreeMixin.distance(tree, target1, target2)
    print(distance)

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    inputTree = ''
    target1   = ''
    target2   = ''

    try:
        opts, args = getopt.getopt(argv, "hi:o:t:")
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
            # explanation
            print("Explanation:")
            print("Determines the distance between two taxa in a newick tree file")
            # input fasta file
            print("\n-i\ttree file:")
            print("\ttree file must be in newick format")
            # first taxa explanation
            print("\n-o\tfirst target taxa:")
            print("\tfirst taxa of interest")
            # second taxa explanation
            print("\n-t\tsecond target taxa:")
            print("\tsecond taxa of interest")
            print("\ndistance will be calculated between taxa specificed using -o and -t arguments")
            sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                inputTree = arg
            else:
                # error message
                print("\n\nInput tree file does not exist!\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-o':
            if arg:
                target1 = arg
                #print("\nTarget1: {}".format(target1))
            else:
                # error message
                print("\n\nInvalid target1\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if arg:
                target2 = arg
                #print("\nTarget2: {}".format(target2))
            else:
                # error message
                print("\n\nInvalid target2\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
    # pass to readTree function
    readTree(inputTree, target1, target2)

if __name__ == '__main__':
    main(sys.argv[1:])