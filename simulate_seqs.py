#!/usr/bin/env python

import sys
import getopt
import os.path
import pyvolve

####################################################################
### Master execute Function                                      ###
### This function executes the main functions and calls other    ###
### subfunctions to calculate pairwise correlations of relative  ###
### evolutionary rates (represented by treeXBL, where X is       ###
### tree1 or tree2)                                              ###
####################################################################
def execute(
    tree, 
    model,
    length,
    out,
    numSim
    ):

    # read in model, tree, and define partition
    pyvolveModel     = pyvolve.Model(model)
    pyvolveTree      = pyvolve.read_tree(file = tree)
    pyvolvePartition = pyvolve.Partition(models = pyvolveModel, size = int(length))

    # create evolver
    my_evolver = pyvolve.Evolver(tree = pyvolveTree, partitions = pyvolvePartition)
    my_evolver()

    print("Simulating sequences...")
    # create simluated sequences
    for i in range(int(numSim)):
        print(str(out) + "." + str(i) + ".fa")
        my_evolver(seqfile = str(out) + "." + 
        	str(model) + "-" + str(i) + ".fa")


####################################################################
### END Functions that read the input files and create output    ###
####################################################################

## Function to print the help message
def help_message(
    ):
    """
    Prints help message
    -------------------
    argv: NA
    """ 

    print("\n\nThis script simulates sequence evolution along a phylogeny.\n")
    print("Usage:")
    print("python script.py -t newick.tre -m model")
    print("\t-l sequence length -n number of simulations")
    print("\t-o output name\n\n")
    print("Nucleotide models include:")
    print("\tNucleotide Models")
    print("\t-----------------")
    print("\tnucleotide = equal equilibrium frequencies and equal mutation rates")
    print("\n\tAmino acid models")
    print("\t-----------------")
    print("\tJTT, WAG, LG, AB, DAYHOFF, MTMAM, or MTREV24")
    print("\n\tMechanistic (dN/dS) codon models")
    print("\t--------------------------------")
    print("\tGY, MG, or codon\n\n")

## Main function that reads in user arguments
def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree    = ''
    model   = ''
    length  = ''
    output  = ''
    numSim  = ''

    try:
        opts, args = getopt.getopt(argv, "ht:m:l:o:n:")
    except getopt.GetoptError:
        # error message
        print("Error\nFor help use -h argument\n")
        sys.exit(2)
    # if no arguments are used print a help message
    if len(opts) == 0:
        # error message
        help_message()
        sys.exit(2)
    # test for arguments
    for opt, arg in opts:
        if opt == '-h':
            ## explanation
            help_message()
            sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified file of arg -t does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-m':
            if arg:
                model = arg
            else:
                # error message
                print("\n\nThe specified input of arg -m does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-l':
            if int(arg):
                length = arg
            else:
                # error message
                print("\n\nThe -l argument must be an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-o':
            if arg:
                output = arg
            else:
                # error message
                print("\n\nPlease provide a prefix for the output files.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-n':
            if int(arg):
                numSim = arg
            else:
                # error message
                print("\n\nThe -n argument must be an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()


    # pass to read inputs function
    execute(
        tree, 
        model,
        length,
        output,
        numSim
        )

####################################################################
### END Functions that print messages and read arguments         ###
####################################################################

if __name__ == '__main__':
    main(sys.argv[1:])
