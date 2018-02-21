#!/usr/bin/env python

import sys
import getopt
import os.path
import os
import re

def create_output_file(
    taxa_list, dict_of_dicts, output_file, all_buscos
    ):
    """
    goes through dict_of_dicts keys and retrevies value per key in taxa list
    
    Parameters
    ----------
    argv: taxa_list
        list of taxa
    argv: dict_of_dicts
        dictionary of USCO keys with values that are dictionaries of taxa id
    argv: output_file
        output file name
    argv: all_buscos
        list of all busco ids
    """
 
    fin_dict = {}
    buscoCN  = ''
    open(output_file, "w")

    # populate empty dictionary with no values and keys
    # of the various buscos
    fin_dict = {key: [] for key in all_buscos}

    with open(output_file, "a") as out:
        header = "buscoIDs\t"+"\t".join(taxa_list)+"\n"
        out.write(header)
        
        for busco in all_buscos:
            for taxa in taxa_list:
                #print(taxa, busco)
                fin_dict[busco].append(dict_of_dicts[busco][taxa])
        
        for key, value in fin_dict.items():
            buscoCN = '\t'.join(map(str, value))
            entry = key + "\t" + buscoCN +"\n"
            out.write(entry)


def create_USCOs_dod(
    buscoDirs, fullTable, output_file
    ):
    """
    Reads buscos and creates a taxon occupancy matrix
    
    Parameters
    ----------
    argv: busco_list
        single column file of busco output dirs
    """

    occupancyD    = {}
    taxa_list     = []
    dict_of_dicts = {}
    all_buscos    = []

    # loop through each full table file
    for directory, table in zip(buscoDirs, fullTable):
        tableFileName = directory+"/"+table
        with open(tableFileName, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    line = re.split(r'\t+', line.rstrip('\n'))
                    if (len(line) >= 3) and (line[1] != "Fragmented"):
                        line[2] = line[2].split('|',1)[0]
                        if line[0] not in occupancyD:
                            occupancyD[line[0]] = [line[2]]
                        else:
                            occupancyD[line[0]].append(line[2])
                        if line[2] not in taxa_list:
                            taxa_list.append(line[2])

    # save all BUSCO ids to list based on the last tableFileName
    for k,v in occupancyD.items():
        all_buscos.append(k)
    # remove duplicates
    all_buscos=list(set(all_buscos))

    # count occurence of each string in list and create a dictionary of dictionaries
    dict_of_dicts = {key: dict((x,value.count(x)) for x in set(value)) for key, value in occupancyD.items()}
    
    # if taxa isnt present, put a 0 for buscoID
    for key, valueD in dict_of_dicts.items():
        for taxa in taxa_list:
            if taxa not in valueD:
                valueD[taxa] = 0

    create_output_file(
    taxa_list, dict_of_dicts, output_file, all_buscos
    )
    

def read_busco_lists(
    busco_list, output_file
    ):
    """
    Reads busco and fasta list into a python list 
    
    Parameters
    ----------
    argv: busco_list
        single column file of busco output dirs
    """

    # initialize lists
    buscoDirs  = []

    # read lists into python lists
    buscoDirs  = [line.rstrip('\n') for line in open(busco_list)]

    # initialize lists for full table busco ouput file
    fullTable  = []

    # get names of full table files
    for ii in buscoDirs:
        fullTableName = ''
        ii           += str('.tsv')
        fullTableName = ii.replace("run_","full_table_")
        fullTable.append(fullTableName)

    create_USCOs_dod(
        buscoDirs, fullTable, output_file
        )

def main(
    argv
    ):
    """
    Reads arguments
    """
    busco_list  = ''
    output_file = ''
    
    try:
        opts, args = getopt.getopt(argv, "hb:o:")
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
            print("help message")
            sys.exit()
        elif opt == '-b':
            if os.path.isfile(arg):
                busco_list = arg
            else:
                # error message
                print("\n\nProvided file of busco tables does not exist! Exiting... \n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-o':
            if not os.path.isfile(arg):
                output_file = arg
            else:
                # error message
                print("\n\nFile with name of output file already exists! Exiting... \n")
                print("For detailed explanation use -h argument\n")
                sys.exit()


    # pass to read_lists function
    read_busco_lists(
        busco_list, output_file
        )

if __name__ == '__main__':
    main(sys.argv[1:])