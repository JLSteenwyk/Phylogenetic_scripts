#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
import re

def calcDesirability(
        alignments, characteristics_arr
        ):
    """
    determines if input fasta is nucleotide or protein

    Parameters
    ----------
    argv: alignments
        list of alignment fasta files
    argv: characteristics_arr
        characteristics of the alignment in an array
    """

    # intiailize list to hold RCVs, alignment lengths, and percentParsimonyInformativeSites
    RCVs=[]
    alignmentLens=[]
    percentParsimonyInformativeSites=[]
    # loop through characteristics_arr and populate empty lists with corresponding values
    for row in characteristics_arr:
        RCVs.append(row[0])
        alignmentLens.append(row[1])
        percentParsimonyInformativeSites.append(row[2])

    ## calculate desirability among RCVs (low values most desirable)
    # initialize list to hold dRCV values
    #print("RCVs")
    dRCVs=[]
    minRCV=min(RCVs)
    maxRCV=max(RCVs)
    scale=1
    for RCV in RCVs:
        dRCV=((RCV - maxRCV)/(minRCV - maxRCV))**scale
        # append dRCV to dRCVs
        dRCVs.append(dRCV)
        #print(dRCV)

    ## calculate desirability for alignment lengths
    # initialize list to hold dRCV values
    #print("Alignment Lengths")
    dAlignmentLens=[]
    minAlignmentLens=min(alignmentLens)
    maxAlignmentLens=max(alignmentLens)
    scale=1
    for alignmentLen in alignmentLens:
        dAlignmentLen=((alignmentLen - minAlignmentLens)/(maxAlignmentLens - minAlignmentLens))**scale
        # append dAlignmentLen to dAlignmentLens
        dAlignmentLens.append(dAlignmentLen)
        #print(dAlignmentLen)

    ## calculate desirability for alignment lengths
    # NOTE percentParsimonyInformativeSites IS PERCENTvariableSITES
    # initialize list to hold dRCV values
    #print("Percent parsimony informative sites")
    dPercentParsimonyInformativeSites=[]
    minPercentParsimonyInformativeSites=min(percentParsimonyInformativeSites)
    maxPercentParsimonyInformativeSites=max(percentParsimonyInformativeSites)
    scale=1
    for percentParsimonyInformativeSite in percentParsimonyInformativeSites:
        dPercentParsimonyInformativeSite=((percentParsimonyInformativeSite - minPercentParsimonyInformativeSites)/(maxPercentParsimonyInformativeSites - minPercentParsimonyInformativeSites))**scale
        # append dAlignmentLen to dAlignmentLens
        dPercentParsimonyInformativeSites.append(dPercentParsimonyInformativeSite)
        #print(dPercentParsimonyInformativeSite)

    # loop through dRCVs and create a new char_arr that has each row as one
    # alignment overall desirability score and then each score as a column
    # in the array
    char_arr=[]
    for i in range(0,len(dRCVs)):
        temp_arr=[]
        val=0
        #dRCVs
        val=dRCVs[i]
        temp_arr.append(val)
        #dAlignmentLengths
        val=dAlignmentLens[i]
        temp_arr.append(val)
        #dPercentParsimonyInformativeSites
        val=dPercentParsimonyInformativeSites[i]
        temp_arr.append(val)
        char_arr.append(temp_arr)

    # calculate arthmetic mean
    overallD=[]
    for row in char_arr:
        val=float(sum(row))/len(row)
        overallD.append(val)

    for alignment, desirabilityScore in zip(alignments, overallD):
        print('{}\t{}'.format(alignment, desirabilityScore))


            



def RCV(
    alignments, state
    ):
    """
    determines if input fasta is nucleotide or protein

    Parameters
    ----------
    argv: alignments
        list of alignment fasta files
    """

    # read alignments into list
    alignments = [line.rstrip('\n') for line in open(alignments)]
    
    # intiailize array to hold characteristics of the alignments
    characteristics_arr=[]

    # loop through alignments
    for alignment in alignments:
        # string to hold all sequences
        concatSeq = ''
        # initialize a counter for the number of sequences in the input fasta file
        numRecords = 0

        # read in records of a fasta file
        records = list(SeqIO.parse(alignment, "fasta"))
        
        # for each record join concatSeq string and sequence as well as keeping track 
        # of the number of records
        for record in records:
            concatSeq  += record.seq
            numRecords += 1

        # dictionary to hold the average occurence of each sequence letter
        averageD = {}
        # loop through the different sequences that appear in the fasta file
        # population dictionary with how many times that sequence appears
        for seq in set(concatSeq):
        	averageD[seq] = (concatSeq.count(seq)/numRecords)

        # determine the length of the alignment
        alignmentLen = 0
        alignmentLen = len(records[0].seq)

        # intiailize list to hold the RCV values per ith taxa 
        # that will later be summed
        indivRCVvalues = []

        # loop through records again and calculate RCV for 
        # each taxa and append to indivRCVvalues
        for record in records:
            # temp holds a temporary value of the numerator before appending
            # to numeratorRCVvalues and then is reassigned to 0 when it goes
            # through the loop again
            temp = 0
            # calculates the absolute value of the ith sequence letter minus the average
            for seqLetter in set(concatSeq):
                temp += abs(record.seq.count(seqLetter)-averageD[seqLetter])
            indivRCVvalues.append(temp/(numRecords*alignmentLen))
        # calculate RCV
        RCV=sum(indivRCVvalues)


        # create a list with all taxa names
        taxa_list = []
        for record in records:
            taxa_list.append(record.id)
        # initialize dictionary to hold sequences for each taxa
        seqDict = {}
        # populate seqDict with individual id (key) and sequence (value)
        for indiv in records:
            if indiv.id in taxa_list:
                seqDict[indiv.id] = indiv.seq

        # set step and window
        step     = 1
        window   = 1
        # intialize variable to hold number of variable sites
        varSites = 0
        # if protein alignment
        if state == 'prot':
            # loop through fasta alignment file one position at a time
            for i in range(0, (int(alignmentLen)+1) - int(step), int(step)):
                # list to hold sequence at that position
                positionSeq = ''
                # loop through individuals in sequence and append the sequence
                # at the position in the loop
                for k, v in seqDict.items():
                    positionSeq += (v[i:i+int(window)])
                # extract only the sequence from the sequence object
                positionSeq = positionSeq._data.upper()
                # remove gaps (-?) and ambiguous sites (NX)
                positionSeq = re.sub('[-?BZJX.]', '', positionSeq)
                # if there are multiple sequences in the position, add 1 to varSites
                if len(set(positionSeq)) > 1:
                    varSites += 1
            # if protein alignment
        elif state == 'nucl':
            # loop through fasta alignment file one position at a time
            for i in range(0, (int(alignmentLen)+1) - int(step), int(step)):
                # list to hold sequence at that position
                positionSeq = ''
                # loop through individuals in sequence and append the sequence
                # at the position in the loop
                for k, v in seqDict.items():
                    positionSeq += (v[i:i+int(window)])
                # extract only the sequence from the sequence object
                positionSeq = positionSeq._data.upper()
                # remove gaps (-?) and ambiguous sites (NX)
                positionSeq = re.sub('[RYWSKMDVHB]', '', positionSeq)
                # if there are multiple sequences in the position, add 1 to varSites
                if len(set(positionSeq)) > 1:
                    varSites += 1

        percentParsimonyInformativeSites=varSites/alignmentLen
        
        # populate characteristics_arr w/ RCV, alignmentLen, percentParsimonyInformativeSites
        temp_arr=[]
        temp_arr.append(RCV)
        temp_arr.append(alignmentLen)
        temp_arr.append(percentParsimonyInformativeSites)
        characteristics_arr.append(temp_arr)

    # pass to calcDesirability
    calcDesirability(alignments, characteristics_arr)


def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    alignment = ''

    try:
        opts, args = getopt.getopt(argv, "hi:c:")
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
            ## explanation
            print("\nCalculate relative nucleotide composition variability (RCV) as defined by Phillips and Penny")
            print("2003, specifically in the manuscript titled 'The root of the mammalian tree inferred from")
            print("whole mitochondrial genomes. Mol Phylogenet Evol. 28:171–185.")
            print("https://www.ncbi.nlm.nih.gov/pubmed/12878457\n")
            print("The formula for RCV is as follows:")
            print("the sum of i=1^n of (|Ai - A*| + |Ti - T*| + |Ci - C*| + |Gi - G*|)/n•t")
            print("where, Ai, Ti, Ci, and Gi are the number of nucleotides for the ith taxon,")
            print("A*, T*, C*, and G* are the averages of each nucleotide across the n taxa")
            print("t is the number of sites and n is the number of taxa.")
            print("\nThe output will a number that is the RCV value for a given alignment.")
            ## options
            # alignment files list
            print("\n-i <alignment fasta file>")
            print("\tA multi-fasta sequence alignment file")
            sys.exit()

        elif opt == '-i':
            if os.path.isfile(arg):
                alignments = arg
            else:
                # error message
                print("\n\nThe specified alignment list (-i) file does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-c':
            if arg:
                state = str(arg)
            else:
                # error message
                print("\n\nMust be nucl or prot for nucleotides or proteins.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()

    # pass to read_config parameter
    RCV(
        alignments, state
        )

if __name__ == '__main__':
    main(sys.argv[1:])