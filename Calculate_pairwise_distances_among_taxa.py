#!/usr/bin/env python

import sys
import getopt
import os.path
import re
from Bio.Phylo.BaseTree import TreeMixin
from Bio import Phylo
import itertools

def readTree(
    inputTree, perms
    ):
    """
    Reads input tree

    Parameters
    ----------
    argv: inputTree
        reads the input tree file
    argv: perms
        all permutations of taxa list with length 2
    """
    # initialize variables to store tree, leaf1, leaf2
    tree    = ''
    target1 = ''
    target2 = ''
    leaf1   = ''
    leaf2   = ''
    # read tree
    tree    = Phylo.read(inputTree, 'newick')
    
    # loop through the pairs
    for pair in perms:
        # initialize loop specific variable of distance between 
        # two taxa
        distance = ''

        # set pair[0] to target1 and pair[1] to target2
        target1=pair[0]
        target2=pair[1]

        # determine if leaves exists
        leaf1 = TreeMixin.find_any(tree, target1)
        if (bool(leaf1)) == False:
            print(target1, "not on tree\nExiting...")
            sys.exit()
        leaf2 = TreeMixin.find_any(tree, target2)
        if (bool(leaf2)) == False:
            print(target2, "not on tree\nExiting...")
            sys.exit()

        # determine leaf1 and leaf2 distance
        distance = TreeMixin.distance(tree, target1, target2)
        #
        print('{}\t{}\t{}'.format(target1, target2, distance))

def create_perms(
    inputTree, taxaList
    ):
    """
    Reads input file into a list and creates 
    every unique permutation of length 2

    Parameters
    ----------
    argv: inputTree
        reads the input tree file
    argv: taxaList
        single col file of taxa
    """

    # initialize variables 
    taxaL = []
    perms = []
    # read file into list
    taxaL = [line.rstrip('\n') for line in open(taxaList)]
    # create permutations
    perms = list(itertools.combinations(taxaL, 2))

    readTree(
        inputTree, perms
        )

def main(
    argv
    ):
    """
    Reads arguments and passes to create_perms variable
    """

    # initialize argument variables
    inputTree = ''
    taxaList  = ''

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
            # explanation
            print("Explanation:")
            print("Determines the distance between two taxa in a newick tree file")
            # input fasta file
            print("\n-i\ttree file:")
            print("\ttree file must be in newick format")
            # first taxa explanation
            print("\n-l\tlist of target taxa:")
            print("\tsingle column file of taxa of interest to")
            print("\tcalculate pairwisedistances for")
            sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                inputTree = arg
            else:
                # error message
                print("\n\nInput tree file does not exist!\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-l':
            if arg:
                taxaList = arg
                #print("\nTarget1: {}".format(target1))
            else:
                # error message
                print("\n\nInvalid target1\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to create_perms function
    create_perms(
        inputTree, taxaList
        )

if __name__ == '__main__':
    main(sys.argv[1:])