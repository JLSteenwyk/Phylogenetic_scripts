#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio import Phylo
import statistics

def determine_long_sequences(
    tree, times
    ):
    """
    determines if terminal tree branches
    are 'times' longer than the median 
    branch length
    
    Parameters
    ----------
    argv: tree
        newick tree file
    argv: times
        the number of times sequences have
        to be longer the median branch to be
        reported
    """

    # read in tree
    tree = Phylo.read(tree, 'newick')

    # initialize list of branch lengths, dictionary to store 
    # terminal names and the branch length, a variable to 
    # hold the threshold value of terminal branch length,
    # and the median value of terminal branch lengths
    allLens    = []
    LenDict    = {}
    threshold  = ''
    median     = ''

    # append terminal branch lengths to allLens
    for terminal in tree.get_terminals():
        allLens.append((terminal.branch_length))
        LenDict[terminal.name] = terminal.branch_length

    # determine median and threshold terminal branch length
    median    = statistics.median(allLens)
    threshold = median*float(times)

    # if terminal branch is 'times' longer than the 
    # median terminal branch, print the sequence name
    for name, length in LenDict.items():
        if length >= threshold:
            print("{}\t{}\t{}\t{}".format(name, length, threshold, median))



def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree  = ''
    times = ''

    try:
        opts, args = getopt.getopt(argv, "ht:n:")
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
            print("\nDetermines the branch length of each branch (internal and terminal) in")
            print("a newick tree file. Then it determines the median branch length and reports")
            print("any terminal branches that are equal to or greater than 'xx' the median terminal")
            print("branch length. 'xx' is defined by the -n parameter. This script was designed")
            print("specifically to identify potentially spurious sequences and subsequently remove")
            print("them from a single gene tree.\n")
            print("Output is formatted in the following way:")
            print("Col 1: Terminal branch name that violates branch length test")
            print("Col 2: Branch length of the taxa mentioned in Col 1")
            print("Col 3: The threshold value of a terminal branch to length")
            print("Col 4: The median branch length across all terminal branches\n")
            print("This filtering method was applied by Xing-xing Shen did in his 332 Saccharomycotina")
            print("project. The project should be cited if this general method is used.\n")
            sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified tree file (-t) does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-n':
            if float(arg):
                times = arg
            else:
                # error message
                print("\n\nThe specified tree file (-t) does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    determine_long_sequences(
        tree, times
        )

if __name__ == '__main__':
    main(sys.argv[1:])