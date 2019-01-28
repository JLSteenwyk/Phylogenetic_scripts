#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo.BaseTree import Clade
from Bio import Phylo
import re


def modTree(
    tree
    ):
    """
    reports IC as internode support 
    so that the tree can be visualized
    in figturee
    
    Parameters
    ----------
    argv: tree
        newick tree file
    """

    filename=tree+".ic.tree"

    # read in tree
    tree = Phylo.read(tree, 'newick')
    for i in tree.get_nonterminals():
        # if there is no comment, continue
        if i.comment is None:
            1
        # if there is a comment switch confidence for IC values
        elif i.comment is not None:
            i.confidence=float(re.sub(r',[^,]*$', '', i.comment))
            i.comment = None
        # Repeat the previous if statement logic but loop through
        # the remainder of the tree
        for ii in i.clades:
            if ii.comment is None:
                1
            elif ii.comment is not None:
                ii.confidence = float(re.sub(r',[^,]*$', '', ii.comment))
                ii.comment = None

    # write to a new tree
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
            print("\nInternode certainty scores as outputted from RAxML are modified to be reported branch support labels.")
            print("The input tree is specified with the -t parameter. Internode certainty values can then be easily viewed")
            print("in FigTree of other similar software.\n\n")
            print("The output file will have the same name as the input file with the suffix 'ic.tree' added to it")
            sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to modTree function
    modTree(
        tree
        )

if __name__ == '__main__':
    main(sys.argv[1:])