#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo.BaseTree import Clade
from Bio import Phylo

def root(
    tree, clade, filename
    ):
    """
    roots tree in newick format
    on a single column list of outgroup
    clade names
    
    Parameters
    ----------
    argv: tree
        newick tree file
    argv: clade
        single column file of outgroup taxa
    argv: filename
        output file name
    """

    # read in tree
    tree = Phylo.read(tree, 'newick')

    # initialize variables for terminal branch length
    clade = [line.rstrip('\n') for line in open(clade)]

    outgroup = [{'name': taxon_name} for taxon_name in clade]

    tree.root(outgroup)

    Phylo.draw_ascii(tree)

    Phylo.write(tree, filename, 'newick')

    

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree  = ''
    clade = ''

    try:
        opts, args = getopt.getopt(argv, "ht:o:")
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
            print("\nRoots tree on specified clade of individuals. The individuals to root the tree on are")
            print("specified in a separate file that contains a single column file list of taxa. This file")
            print("is specified using the -o parameter. The input tree is specified with the -t parameter.")
            print("The output is a newick tree file with the same name as the original tree file with")
            print(".rooted.tre appended to the end of the file name.")
            sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                tree = arg
                filename=tree+".rooted.tre"
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-o':
            if os.path.isfile(arg):
                clade = arg
            else:
                # error message
                print("\n\nThe specified outgroup clade file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    root(
        tree, clade, filename
        )

if __name__ == '__main__':
    main(sys.argv[1:])