#!/usr/bin/env python

import sys
from sys import stdout
import getopt
import os.path
import itertools
from Bio import SeqIO
import random


def create_combos(
    fasta, length, prefix, replicates, seed
    ):
    """
    creates n replicates of taxa combinations based on the seed and length

    Parameters
    ----------
    argv: taxa
        single column file with list of taxa
    argv: fasta
        fasta file to subset
    """
    
    # create list of taxa sequences
    fastaL = list(SeqIO.parse(fasta, "fasta"))
    taxaL = []
    for i in range(0,len(fastaL)):
        taxaL.append(fastaL[i].name)

    # print log message
    print("\n"+"-"*(len("- General features -")))
    print("| General features |")
    print("-"*(len("- General features -")))
    print("Total number of taxa:", len(taxaL))
    print("Number of taxa per subset file:", length)
    print("Number of replicates:", replicates)
    print("Seed set:", seed)

    # set seed
    random.seed(int(seed))
    # initialize taxa array
    taxa_arr = []
    # create an array of replicates
    for replicate in range(int(replicates)):
        # initialize list to hold the random seed
        tempSeed = []
        # generate a random seed
        tempSeed = random.randint(1, 10000000001) #10,000,000,000
        # set random seed
        random.seed(int(tempSeed))
        # intialize list to hold sample
        temp_arr = []
        # generate sample
        temp_arr = random.sample(taxaL, int(length))
        # append sample to taxa array
        taxa_arr.append(temp_arr)

    # save fasta file to dictionary
    fastaD = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    ## create replicate fasta files
    # initialize counter
    cnt = 1
    # loop through replicates
    print("\n")
    for replicate in taxa_arr:
        # create list of records
        records = []
        # create output file name
        out = prefix + "." + str(cnt) + ".fa"
        # for individual fasta entry header in replicate
        for indiv in replicate:
            records.append(fastaD[indiv])
        # print to stdout
        percent = (cnt/replicates)*100
        sys.stdout.write("\rWriting {} of {} ({:.2f}%) fasta files".format(cnt, replicates, percent))
        # write output fasta file
        SeqIO.write(records, out, "fasta")
        sys.stdout.flush()
        cnt+=1
    
    print("\n")


def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    fasta  = ''
    length = ''
    replicates = ''
    prefix = ''
    seed = ''

    try:
        opts, args = getopt.getopt(argv, "hi:l:p:r:s:")
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
            print("\nThis script is designed to make subsets of a given fasta file. Subsets are created based")
            print("on taxa and not the alignment such that a fraction of the total taxa represented in the input")
            print("fasta file are used to create the resulting fasta files. This script samples without replacement.\n")
            print("Required arguments include:")
            print("-i: input fasta file")
            print("-l: subset length (i.e., how many taxa will be in each subset")
            print("-r: number of replicates/subset files to create")
            print("-s: seed")
            print("-p: output prefix\n")
            print("For example, the command")
            print("python create_taxon_bootstrap_alignments.py -i input.fa -l 50 -r 100 -s 1 -p subset")
            print("will subset 50 taxa from input.fa a total of 100 times. Output files will be subset.{1..100}.fa")
            print("\n")
            sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                fasta = arg
            else:
                # error message
                print("\n\nThe specified fasta file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()   
        elif opt == '-l':
            try:
                if arg:
                    length = int(arg)
                else:
                    # error message
                    print("\n\nThe specified length argument must be an integer.\n")
                    print("For detailed explanation use -h argument\n")
                    sys.exit()
            except ValueError:
                print("\n\nThe specified length argument must be an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-p':
            if arg:
                prefix = arg
            else:
                # error message
                print("\n\nThe required prefix argument was not provided.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-r':
            try:
                if arg:
                    replicates = int(arg)
                else:
                    # error message
                    print("\n\nThe specified replicates argument must be an integer.\n")
                    print("For detailed explanation use -h argument\n")
                    sys.exit()
            except ValueError:
                print("\n\nThe specified replicates argument must be an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-s':
            try:
                if arg:
                    seed = int(arg)
                else:
                    # error message
                    print("\n\nThe specified seed argument must be an integer.\n")
                    print("For detailed explanation use -h argument\n")
                    sys.exit()
            except ValueError:
                print("\n\nThe specified seed argument must be an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()      

    # pass to modTree function
    create_combos(
        fasta, length, prefix, replicates, seed 
        )

if __name__ == '__main__':
    main(sys.argv[1:])