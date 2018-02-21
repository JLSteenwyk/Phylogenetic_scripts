#!/usr/bin/env python

import sys
import getopt
import os.path
import os
import re
from Bio import SeqIO
from datetime import datetime

def create_fas_per_USCO(
    buscoDirs, taxon_occupancy,
    taxaNames, sequence_type, 
    fullTable, final_USCO_list,
    startTime
    ):
    """
    creates fasta files per USCO that pass the 
    taxon occupancy criterion

    Parameters
    ----------
    argv: buscoDirs
        list of busco output dirs
    argv: taxon_occupancy
        value between 0 and 1 for taxon occupancy per USCO
    argv: taxaNames
        list of taxa names
    argv: sequence_type
        prot or nucl
    argv: fullTable
        list of full table busco output file names
    argv: final_USCO_list
        list of USCOs that pass occupancy criteria
    """

    # initialize variable for fastaFiles len
    fastaFilesLength = ''

    # loop through list of USCOs that pass taxon occupancy
    fastaFilesLength = len(taxaNames)
    total      = len(final_USCO_list)

    # print log message
    print("creating fasta files...")
    print("...")

    for USCOid in final_USCO_list:
        # make fasta file to populate
        USCOfasta = USCOid+".fa"
        with open(USCOfasta, "a+") as output_handle:

            # loop through busco output directories and taxa names
            for directory, taxa in zip(buscoDirs, taxaNames):
                
                # create pathway to busco file dependent on whether protein or nucleotide sequences are desired
                indivBUSCO         = ''
                if sequence_type   == 'prot':
                    indivBUSCO     = directory+"/single_copy_busco_sequences/"+USCOid+".faa"
                elif sequence_type == 'nucl':
                    indivBUSCO     = directory+"/single_copy_busco_sequences/"+USCOid+".fna"
                
                # check if indivBUSCO exists and append it to USCOfasta
                if os.path.isfile(indivBUSCO):
                    # open indivBUSCO and append to USCOfasta
                    format     = "fasta"
                    handle     = open(indivBUSCO)
                    fasta      = SeqIO.parse(handle, format)
                    for record in fasta:
                        entry  = ''
                        Seq    = str(record.seq)
                        entry += ">"+taxa+"|"+USCOid+"\n"+Seq+"\n"
                        output_handle.write(entry)
                        break
                # if file does not exist, skip it because it is not a single copy BUSCO
                else:
                    1

    # print log message
    print("Complete!\n")
    print("Total time:", datetime.now() - startTime)
                
def determine_USCOs(
    buscoDirs, taxon_occupancy,
    taxaNames, sequence_type, 
    fullTable, startTime
    ):
    """
    Uses buscoDirs and fullTable to determine which USCOs
    pass the taxon occupancy criterion and creates a list
    of which USCOs pass the occupnacy criterion

    Parameters
    ----------
    argv: buscoDirs
        list of busco output dirs
    argv: taxon_occupancy
        value between 0 and 1 for taxon occupancy per USCO
    argv: taxaNames
        list of taxa names
    argv: sequence_type
        prot or nucl
    argv: fullTable
        list of full table busco output file names
    argv: startTime
        time script was started
    """

    # initialize dictionaries to store USCOs per busco run
    USCOD           = {}
    length_dict     = {}
    final_USCO_list = []
    
    # loop through each full table file and save taxa name to USCOD if complete (i.e., single copy)
    for directory, table, taxa in zip(buscoDirs, fullTable, taxaNames):
        tableFileName = directory+"/"+table
        with open(tableFileName, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    line = re.split(r'\t+', line.rstrip('\n'))
                    if line[1] == "Complete":
                        if line[0] in USCOD:
                            USCOD[line[0]].append([taxa])
                        else:
                            USCOD[line[0]] = [taxa]
    
    # determine how many taxa have USCO in single copy                        
    length_dict = {key: len(value) for key, value in USCOD.items()}

    # create a list BUSCO ids that meet taxon occupancy cut off
    for k, v in length_dict.items():
        if v >= (float(taxon_occupancy)*len(taxaNames)):
            final_USCO_list.append(k)

    print("\n"+"-"*len("- Log updates -"))
    print("| Log updates |")
    print("-"*len("- Log updates -"))
    print(len(taxaNames), "taxa")
    print(len(final_USCO_list), "USCOs pass occupancy")

    create_fas_per_USCO(
        buscoDirs, taxon_occupancy,
        taxaNames, sequence_type, 
        fullTable, final_USCO_list,
        startTime
        )

def read_list(
    busco_list, taxon_occupancy, 
    taxa_list, sequence_type,
    startTime
    ):
    """
    Reads busco directories and taxa names into python lists.
    Also creates a list of full table files
    
    Parameters
    ----------
    argv: busco_list
        single column file of busco output dirs
    argv: taxon_occupancy
        taxon occupancy value between 0 and 1
    argv: taxa_list
        single column file of taxa names
    argv: sequence_type
        prot or nucl
    argv: startTime
        time script was started
    """

    # initialize lists
    buscoDirs  = []
    taxaNames  = []

    # read lists into python lists
    buscoDirs  = [line.rstrip('\n') for line in open(busco_list)]
    taxaNames  = [line.rstrip('\n') for line in open(taxa_list)]

    # initialize lists for full table busco ouput file
    fullTable  = []

    # get names of full table files
    for ii in buscoDirs:
        fullTableName = ''
        ii           += str('.tsv')
        fullTableName = ii.replace("run_","full_table_")
        fullTable.append(fullTableName)

    determine_USCOs(
        buscoDirs, taxon_occupancy,
        taxaNames, sequence_type, 
        fullTable, startTime
        )

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # start time of script
    startTime = datetime.now()

    # initialize argument variables
    busco_list      = ''
    taxon_occupancy = ''
    taxa_list       = ''
    sequence_type   = ''

    try:
        opts, args = getopt.getopt(argv, "hb:c:o:t:")
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
            print("\nCreate fasta files for every busco gene that have sufficient taxon occupancy.")
            print("This script is designed to facilitate phylogenomic analysis using the results")
            print("obtained from BUSCO. More specifically, BUSCOs have been called from genome sequences.")
            print("That is BUSCO run was used with the '-m genome' parameter. BUSCO fastas can be made with")
            print("either the protein or nucleotide files.")
            ## options
            # busco_out_list
            print("\n-b <list of busco dirs>")
            print("\tSingle column file of busco output dirs. Order of busco dirs should match")
            print("\tthe order of the taxa names file specified with the -t parameter.")
            # busco_out_list
            print("\n-c <use nucleotide or protein usco fasta>")
            print("\tArgument can be either 'prot' or 'nucl' to specify if the resulting BUSCO fasta files should")
            print("\tbe made using protein sequences or nucleotide sequences.")
            # taxon occupancy
            print("\n-o <taxon occupancy>")
            print("\tTaxon occupancy for USCOs value should be between 0 and 1. BUSCOs with taxon occupancy")
            print("\tthat fits this parameter will be created. Specifically, if the user specifies 0.75, only")
            print("\tBUSCOs fasta files with at least 75 percent of the taxa will be made.")
            # taxa name list
            print("\n-t <list of taxa names>")
            print("\tSingle column file of taxa names for final usco fasta file. That is, if the user")
            print("\tspecify 'S_cerevisiae' the usco fa will have ''>S_cerevisiae' followed by")
            print("\tthe sequence of S. cerevisiae for that particular usco. Order of taxa should")
            print("\tbe the same as the list of busco dirs. No spaces should be used for this parameter.\n")
            sys.exit()

        elif opt == '-b':
            if os.path.isfile(arg):
                busco_list = arg
            else:
                # error message
                print("\n\nThe specified busco_list (-b) file does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-c':
            if arg in ('prot', 'nucl'):
                sequence_type = arg 
            else:
                # error message
                print("\n\nThe specified input for (-c) should be 'prot' or 'nucl'.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-o':
            if 0 <= float(arg) <= 1:
                taxon_occupancy = arg
            else:
                # error message
                print("\n\nThe specified taxon occupancy is not between 0 and 1.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                taxa_list = arg
            else:
                # error message
                print("\n\nThe specified configuration file does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()

    # print log output
    print("\n"+"-"*(len("- Argument parameters -")))
    print("| Argument parameters |")
    print("-"*(len("- Argument parameters -")))
    print("-b", busco_list)
    print("-c", sequence_type)
    print("-o", taxon_occupancy)
    print("-t", taxa_list)

    # pass to read_config parameter
    read_list(
        busco_list, taxon_occupancy, 
        taxa_list, sequence_type,
        startTime
        )

if __name__ == '__main__':
    main(sys.argv[1:])