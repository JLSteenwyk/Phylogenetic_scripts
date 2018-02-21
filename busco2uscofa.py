#!/usr/bin/env python

import sys
import getopt
import os.path
import os
import subprocess
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import configparser

# def add_missing_taxa_and_write_final_files(
#     buscoDirs, fastaFiles,
#     fullTable, taxon_occupancy,
#     final_USCO_list, mafft_path,
#     trimAl_path, concat_fa, 
#     USCO_partition
#     ):
#     """
#     aligns and trims USCOs

#     Parameters
#     ----------
#     argv: buscoDirs
#         list of busco output dirs
#     argv: fastaFiles
#         list of of fasta files
#     argv: fullTable
#         list of full table busco output file names
#     argv: taxon_occupancy
#         value between 0 and 1 for taxon occupancy per USCO
#     argv: final_USCO_list
#         list of USCO ids that pass taxon occupancy criterion
#     argv: mafft_path
#         pathway to mafft program
#     argv: trimAl_path
#         pathway to trimAl program
#     argv: concat_fa
#         output file name for concat fa
#     argv: USCO_partition
#         output file name for the partition file
#     """

#     # initialize headerNames to store names of all indivs per fa
#     # , concat.fa file, dictionary to populate concat.fa with
#     headerNames = []
#     concatD     = {}
#     firstLen    = 1
#     secondLen   = 0

#     # loop through fasta files, read first line, save id
#     for fastas in fastaFiles:
#         with open(fastas, 'r') as f:
#             first_line = f.readline().strip()
#             first_line = re.sub(r'>', '', first_line)
#             first_line = re.sub(r'\|.*', '', first_line)
#             headerNames.append(first_line)
    
#     # populate empty dictionary with no values and keys
#     # of the various indivs
#     concatD = {key: [] for key in headerNames}

#     open(USCO_partition, "w")
#     # loop through USCO trimal files
#     for USCOid in final_USCO_list:
#         # list for keep track of taxa per USCO
#         USCOtaxa    = []
#         # used to store length of USCO alignment
#         USCOlen     = ''
#         USCOtrimal  = USCOid+".trimAl"
#         # used to store missing sequence to append
#         missing_seq = ''
#         # open trimal file
#         records = list(SeqIO.parse(USCOtrimal, "fasta"))
#         USCOlen = len(records[0].seq)

#         # create lists of USCOtaxa and missing_taxa
#         for record in records:
#             # remove '|' till end of line and '>'
#             record.id = re.sub(r'\|.*', '', record.id)
#             record.description = re.sub(r'\|.*', '', record.id)
#             # append taxa per USCO to USCOtaxa
#             USCOtaxa.append(record.id)
#             # determine missing taxa
#             missing_taxa = list(set(USCOtaxa)^set(headerNames)) 

#         # create record for missing taxa
#         for taxa in missing_taxa:
#             missingSeq = USCOlen*'?'
#             taxaRecord = SeqRecord(Seq(missingSeq,IUPAC.protein), 
#                 id = taxa, name = taxa, description = taxa)
#             concatD[taxa].append(taxaRecord.seq)
#         # create record for present data
#         for record in records:
#             if record.id in USCOtaxa:
#                 concatD[record.id].append(record.seq)
#         # create partition file
#         with open(USCO_partition, "a") as f:
#             # second value in partition file
#             secondLen += USCOlen
#             entry = "AUTO, "+str(USCOid)+"="+str(firstLen)+"-"+str(secondLen)+"\n"
#             f.write(str(entry))
#             # add to first value for partition file
#             firstLen += USCOlen

#     # join seqs of genes in value (list) and write to concat_fa   
#     with open(concat_fa, "w") as final_fasta_file: 
#         for x in concatD:
#             concatenated = Seq("", IUPAC.protein)
#             for s in concatD[x]:
#                 concatenated += s
#             concatD[x] = concatenated
#             entry = '>'+x+"\n"+concatD[x]+"\n"
#             final_fasta_file.write(str(entry))

# def align_and_trim(
#     buscoDirs, fastaFiles,
#     fullTable, taxon_occupancy,
#     final_USCO_list, mafft_path,
#     trimAl_path, concat_fa, 
#     USCO_partition
#     ):
#     """
#     aligns and trims USCOs

#     Parameters
#     ----------
#     argv: buscoDirs
#         list of busco output dirs
#     argv: fastaFiles
#         list of of fasta files
#     argv: fullTable
#         list of full table busco output file names
#     argv: taxon_occupancy
#         value between 0 and 1 for taxon occupancy per USCO
#     argv: final_USCO_list
#         list of USCO ids that pass taxon occupancy criterion
#     argv: mafft_path
#         pathway to mafft program
#     argv: trimAl_path
#         pathway to trimAl program
#     argv: concat_fa
#         output file name for concat fa
#     argv: USCO_partition
#         output file name for the partition file
#     """

#     # initialize USCO fasta, mafft, and trimal variables
#     USCOfasta  = ''
#     USCOmafft  = ''
#     USCOlength = ''

#     # align each USCOid fasta file
#     USCOlength = len(final_USCO_list)
#     for USCOid in final_USCO_list:
#         USCOfasta  = USCOid+".fa"
#         USCOmafft  = USCOid+".mafft"
#         USCOtrimal = USCOid+".trimAl"
#         # align
#         with open(USCOmafft, 'w') as f:
#             subprocess.call([mafft_path, '--reorder', '--bl', '62', 
#                 '--op', '1.0', '--maxiterate', '1000', '--retree', '1', 
#                 '--genafpair', '--quiet', USCOfasta], stdout = f)

#         # check 25 times if USCOs got aligned
#         for i in range(1,26):
#             # if the mafft file has size 0, realign the fa file
#             if os.stat(USCOmafft).st_size == 0:
#                 # align
#                 with open(USCOmafft, 'w') as f:
#                     subprocess.call([mafft_path, '--reorder', '--bl', '62', 
#                         '--op', '1.0', '--maxiterate', '1000', '--retree', '1', 
#                         '--genafpair', '--quiet', USCOfasta], stdout = f)
        
#         # trim
#         with open(USCOtrimal, 'w') as f:
#             subprocess.call([trimAl_path, '-in', USCOmafft, '-out', 
#                 USCOtrimal, '-automated1'], stdout = f)

#         # check 25 times if USCOs got trimmed
#         for i in range(1,26):
#             # if the mafft file has size 0, realign the fa file
#             if os.stat(USCOtrimal).st_size == 0:
#                 # trim
#                 with open(USCOtrimal, 'w') as f:
#                     subprocess.call([trimAl_path, '-in', USCOmafft, '-out', 
#                         USCOtrimal, '-automated1'], stdout = f)
        
#     add_missing_taxa_and_write_final_files(
#         buscoDirs, fastaFiles,
#         fullTable, taxon_occupancy,
#         final_USCO_list, mafft_path,
#         trimAl_path, concat_fa, 
#         USCO_partition
#         )




def create_fas_per_USCO(
    buscoDirs, fastaFiles,
    fullTable, taxon_occupancy,
    final_USCO_list, mafft_path,
    trimAl_path, concat_fa, 
    USCO_partition
    ):
    """
    creates fasta files per USCO that passes the 
    taxon occupancy criterion

    Parameters
    ----------
    argv: buscoDirs
        list of busco output dirs
    argv: fastaFiles
        list of of fasta files
    argv: fullTable
        list of full table busco output file names
    argv: taxon_occupancy
        value between 0 and 1 for taxon occupancy per USCO
    argv: final_USCO_list
        list of USCO ids that pass taxon occupancy criterion
    argv: mafft_path
        pathway to mafft program
    argv: trimAl_path
        pathway to trimAl program
    argv: concat_fa
        output file name for concat fa
    argv: USCO_partition
        output file name for the partition file
    """

    # initialize variable for fastaFiles len
    fastaFilesLength = ''
    # create empty files to populate fasta seqs to
    for USCOid in final_USCO_list:
        USCOfasta = USCOid+".fa"
        open(USCOfasta, "w")

    # loop through list of USCOs that pass taxon occupancy
    fastaFilesLength = len(fastaFiles)
    for directory, fastas, table in zip(buscoDirs, fastaFiles, fullTable):
        # open busco output table
        tableFileName = directory+"/"+table
        # check if tableFileName exists and exit if not
        if os.path.isfile(tableFileName):
            1
        else:
            print(tableFileName)
            print("File does not exist. Check order of busco_out and fasta_files")
            sys.exit()
        # open busco full table and append to USCOfasta if complete
        with open(tableFileName, "r") as f:
            # open fasta file
            format     = "fasta"
            handle     = open(fastas)
            fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, format))
            # loop through busco output table and get gene ID associated with USCO ID
            for line in f:
                if not line.startswith("#"):
                    line = re.split(r'\t+', line.rstrip('\n'))
                    # if USCO id is present in single copy, add it to fasta file of USCOs
                    if (line[0] in final_USCO_list) and (line[1] == "Complete"):
                        USCOfasta = line[0]+".fa"
                        with open(USCOfasta, "a") as output_handle:
                            SeqIO.write(fasta_dict[line[2]], output_handle, format)

    # align_and_trim(
    #     buscoDirs, fastaFiles,
    #     fullTable, taxon_occupancy,
    #     final_USCO_list, mafft_path,
    #     trimAl_path, concat_fa, 
    #     USCO_partition
    #     )                            

def determine_USCOs(
    buscoDirs, fastaFiles,
    fullTable, taxon_occupancy,
    mafft_path, trimAl_path, 
    concat_fa, USCO_partition
    ):
    """
    Uses buscoDirs and fullTable to determine which USCOs
    pass the taxon occupancy criterion and creates a list
    of which USCOs pass the occupnacy criterion

    Parameters
    ----------
    argv: buscoDirs
        list of busco output dirs
    argv: fastaFiles
        list of of fasta files
    argv: fullTable
        list of full table busco output file names
    argv: taxon_occupancy
        value between 0 and 1 for taxon occupancy per USCO
    argv: mafft_path
        pathway to mafft program
    argv: trimAl_path
        pathway to trimAl program
    argv: concat_fa
        output file name for concat fa
    argv: USCO_partition
        output file name for the partition file
    """

    # initialize dictionaries to store USCOs per busco run
    USCOD           = {}
    length_dict     = {}
    final_USCO_list = []
    
    # loop through each full table file and save complete USCOs to completeD
    for directory, fastas, table in zip(buscoDirs, fastaFiles, fullTable):
        tableFileName = directory+"/"+table
        with open(tableFileName, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    line = re.split(r'\t+', line.rstrip('\n'))
                    if line[1] == "Complete":
                        if line[0] in USCOD:
                            USCOD[line[0]].append(line[2])
                        else:
                            USCOD[line[0]] = [line[2]]
    # determine how many taxa have USCO in single copy                        
    length_dict = {key: len(value) for key, value in USCOD.items()}
    
    # create a list BUSCO ids that meet taxon occupancy cut off
    for k, v in length_dict.items():
        if v >= (float(taxon_occupancy)*len(fastaFiles)):
            final_USCO_list.append(k)

    create_fas_per_USCO(
        buscoDirs, fastaFiles,
        fullTable, taxon_occupancy,
        final_USCO_list, mafft_path,
        trimAl_path, concat_fa, 
        USCO_partition
        )

def read_lists(
    busco_list, fasta_list,
    taxon_occupancy, mafft_path,
    trimAl_path, concat_fa, 
    USCO_partition
    ):
    """
    Reads busco and fasta list into a python list 
    
    Parameters
    ----------
    argv: busco_list
        single column file of busco output dirs
    argv: fasta_list
        single column file of fasta files
    argv: taxon_occupancy
        value between 0 and 1 for taxon occupancy per USCO
    argv: mafft_path
        pathway to mafft program
    argv: trimAl_path
        pathway to trimAl program
    argv: concat_fa
        output file name for concat fa
    argv: USCO_partition
        output file name for the partition file
    """

    # initialize lists
    buscoDirs  = []
    fastaFiles = []

    # read lists into python lists
    buscoDirs  = [line.rstrip('\n') for line in open(busco_list)]
    fastaFiles = [line.rstrip('\n') for line in open(fasta_list)]

    # initialize lists for full table busco ouput file
    fullTable  = []

    # get names of full table files
    for ii in buscoDirs:
        fullTableName = ''
        ii           += str('.tsv')
        fullTableName = ii.replace("run_","full_table_")
        fullTable.append(fullTableName)

    determine_USCOs(
        buscoDirs, fastaFiles,
        fullTable, taxon_occupancy,
        mafft_path, trimAl_path,
        concat_fa, USCO_partition
        )

def read_config(
    config_file
    ):
    """
    Reads config file

    Parameters
    ----------
    argv: configuration file
        provides all necessary parameters to run
    """
    busco_list      = ''
    fasta_list      = ''
    taxon_occupancy = ''
    mafft_path      = ''
    trimAl_path     = ''
    concat_fa       = ''
    USCO_partition  = ''
    

    settings = configparser.ConfigParser()
    settings.read(config_file)

    # read and check busco_out path
    busco_list = settings.get('input_files', 'busco_out_list')
    if os.path.isfile(busco_list):
        1
    else:
        print("No or incorrect pathway to busco_out provided in configuration file!\nExiting now...")    
        sys.exit()

    # read and check fasta_files path
    fasta_list = settings.get('input_files', 'fasta_files_list')
    if os.path.isfile(fasta_list):
        1
    else:
        print("No or incorrect pathway to fasta_files provided in configuration file!\nExiting now...")    
        sys.exit()

    # read and check taxon_occupancy is between 0 and 1
    taxon_occupancy = settings.get('parameters', 'occupancy')
    if 0 <= float(taxon_occupancy) <= 1:
        1
    else:
        print("No or incorrect value for occupancy parameter!\nExiting now...")    
        sys.exit()

    # read and check mafft path
    mafft_path = settings.get('programs', 'mafft')
    if os.path.isfile(mafft_path):
        1
    else:
        print("No or incorrect pathway to mafft provided in configuration file!\nExiting now...")    
        sys.exit()

    # read and check trimal path
    trimAl_path = settings.get('programs', 'trimAl')
    if os.path.isfile(trimAl_path):
        1
    else:
        print("No or incorrect pathway to trimal provided in configuration file!\nExiting now...")    
        sys.exit()

    # read output concatenation file path
    concat_fa = settings.get('output_files', 'concat_fasta')
    if concat_fa:
        1
    else:
        print("No concat_fasta specified in configuration file!\nExiting now...")    
        sys.exit()

    # read and check partition_file path
    USCO_partition = settings.get('output_files', 'partition_file')
    if USCO_partition:
        1
    else:
        print("No partition_file specified in configuration file!\nExiting now...")    
        sys.exit()


    # pass to read_lists function
    read_lists(
        busco_list, fasta_list,
        taxon_occupancy, mafft_path,
        trimAl_path, concat_fa,
        USCO_partition
        )



def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    config_file = ''

    try:
        opts, args = getopt.getopt(argv, "hc:t")
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
            # options
            print("\n\n-c configuration file")
            print("\tfile that contains all necessary parameters and file names")
            print("\n-t print configuration template")
            print("\tprints a template configuration file for the user to fill out")
            print("\n")
            print("Contents of configuration file")
            print("------------------------------")
            ## programs
            # mafft
            print("mafft: <pathway to mafft program>")
            print("\tpathway to mafft program. Version(s) tested: 7.294b")
            # trimal
            print("\ntrimAl: <pathway to trimAl program>")
            print("\tpathway to trimal program. Version(s) tested: 1.4.rev11")
            ## input_files explanation
            # busco_out_list
            print("\nbusco_out_list: <list of busco dirs>")
            print("\tsingle column file of busco output dirs")
            print("\torder should match the list of fasta file")
            print("\tspecified in fasta_files parameter")
            # fasta_files_list explanation
            print("\nfasta_files_list: <list of fasta files>")
            print("\tsingle column file of fasta files")
            print("\torder should match the list of busco dirs")
            print("\tspecified in busco_out parameter")
            print("\tNOTE: fasta file headers per fa file")
            print("\tshould be named in the following convention...")
            print("\t>ID|1")
            print("\t>ID|2")
            print("\t>ID|3")
            print("\t> ...")
            print("\t>ID|N")
            ## parameters 
            # occupancy explanation
            print("\noccupancy: <taxon occupancy>")
            print("\ttaxon occupancy for USCOs")
            print("\tvalue should be between 0 and 1\n")
            sys.exit()
        elif opt == '-c':
            if os.path.isfile(arg):
                config_file = arg
            else:
                # error message
                print("\n\nThe specified configuration file does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-t':
            print("\n[programs]")
            print("mafft: pathway to mafft")
            print("trimAl: pathway to trimAl")
            print("\n[input_files]")
            print("busco_out_list: single column file with names of busco output dirs")
            print("fasta_files_list: single column file with names of fasta files")
            print("\n[output_files]")
            print("concat_fasta: output concatenation fasta file name")
            print("partition_file: RAxML partition file")
            print("\n[parameters]")
            print("occupancy: taxon occupancy per USCO, a value between 0 and 1\n")
            sys.exit()

    # pass to read_config parameter
    read_config(
        config_file
        )

if __name__ == '__main__':
    main(sys.argv[1:])