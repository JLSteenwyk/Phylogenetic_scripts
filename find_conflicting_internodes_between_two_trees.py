#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio import Phylo
from pprint import pprint
import numpy


def tabulate_names(tree):
    
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            clade.name = '%d_%s' % (idx, clade.name)
        else:
            clade.name = str(idx)
        names[clade.name] = clade
    return(names)

def compare_trees(
    t1, t2
    ):
    """
    compares and determine which internodes conflict between
    a reference tree (t1) and a second tree (t2)
    
    Parameters
    ----------
    argv: t1
        reference tree
    argv: t2
        second tree
    """

    # read in tree
    t1 = Phylo.read(t1, 'newick')
    t2 = Phylo.read(t2, 'newick')
    t1.rooted=True

    tabulate_names(t1)
    tabulate_names(t2)
    
    for I1, I2 in zip(t1.get_nonterminals(), t2.get_nonterminals()):
        if I1.clades[0] != I2.clades[0]:
            print("conflict")
            print(vars(I1.clades)) #, I1.clades, I2.clades)
        break
    

    # print("t1")
    # Phylo.draw_ascii(t1)
    # print("t2")
    # Phylo.draw_ascii(t2)

    # for I1, I2 in zip(t1.find_clades(terminal=False, order='level'), t2.find_clades(terminal=False, order='level')):
    #     for I1child, I2child in zip(I1.clades, I2.clades):
    #         #print(vars(I1child))
    #         if I1child.name != I2child.name:
    #             print("{}\t{}\t{}".format("C", I1child.name, I2child.name))
    #         else:
    #             print("{}\t{}\t{}".format("P", I1child.name, I2child.name))

    #         #adjmat[lookup[parent], lookup[child]] = 1
    #     #for child in parent.clades:
    #         #adjmat[lookup[parent], lookup[child]] = 1(t1)



def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree  = ''

    try:
        opts, args = getopt.getopt(argv, "hr:t:")
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
            print("\nCompares two trees in newick format and reports differences")
            sys.exit()
        elif opt == '-r':
            if os.path.isfile(arg):
                t1 = arg
            else:
                # error message
                print("\n\nThe specified reference tree file (-r) does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                t2 = arg
            else:
                # error message
                print("\n\nThe specified tree file (-t) does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    compare_trees(
        t1, t2
        )

if __name__ == '__main__':
    main(sys.argv[1:])