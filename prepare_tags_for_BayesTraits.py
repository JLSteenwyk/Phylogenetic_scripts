#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio.Phylo.BaseTree import Clade
from Bio import Phylo
from pprint import pprint

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
    filename=tree+".internodeLabels.tree"

    # internode label starter
    label=1

    # read in tree and create tags
    tree = Phylo.read(tree, 'nexus')
    # loop through internal nodes
    for i in tree.get_nonterminals():
        # create temp array to hold 'addtag, nodexyz, and children of node names'
        temp=[]
        temp.append("addtag")
        nodeID="node"+str(label)
        temp.append(nodeID)
        # for each internal node, get the children tips in the tree and append them to the temp list
        for ii in i.get_terminals():
            temp.append(ii.name)
        # prints lines for bayesTraits
        print(*temp, sep=' ')
        print("addMRCA", nodeID, nodeID, sep=' ')
        # replace the confidence value with nodeID in the phylogeny.
        # This is for the additional newick tree that gets written out
        i.confidence=nodeID
        # add one for the next internode label
        label+=1
    #Phylo.write(tree, filename, 'nexus')

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
            print("\nCreates tags for BayesTraits run/parameters file by looping through all internal")
            print("nodes and identifying children tips from there. Nodes are labeled according to")
            print("how they are encoutnered in the tree file. The tree file should be a nexus file")
            print("and be the same one as is used as input for BayesTraits. To specify a tree file, use")
            print("option -t.")
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