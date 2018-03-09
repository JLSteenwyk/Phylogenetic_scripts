# Scripts to enable phylogenetic/phylogenomic analyses

This repository houses numerous scripts to facilitate phylogenetic/phylogenomic sequence analyses

## Scripts

A short description of each script

### Calculate_distance_between_two_taxa.py
Calculates phylogenetic distance between two taxa in a newick tree file.
Taxa names are provided as arguments. <br />
For detailed information use the -h argument <br />
```
python3
|- sys
|- getopt
|- os.path
|- re
|- Bio
   |- Phylo
      |- BaseTree
         |- TreeMixin
```
Basic usage: python Calculate_distance_between_two_taxa.py -i newick_tree_file -o taxa1 -t taxa2 <br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### Calculate_pairwise_distances_among_taxa.py
Calculates all pairwise phylogenetic distances between two taxa in a newick tree file.
Taxa names are provided as a secondary file where taxa names are a single column. <br />
For detailed information use the -h argument <br />
```
python3
|- sys
|- getopt
|- os.path
|- re
|- Bio
   |- Phylo
      |- BaseTree
         |- TreeMixin
|- itertools
```
Basic usage: python Calculate_pairwise_distances_among_taxa.py -i newick_tree_file -l target_taxa_list <br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### calculate_clade1_clade2_branch_len.bash
Calculate the internode branch length of the branch that leads up to clade 1 and clade 2.
Input is a file of newick trees, which is the same as the input for ASTRAL coalescence based
phylogenetic inference. <br />
Variables Clade1, Clade2, and All_Clade are hardcoded and should be changed for each use.<br />
```
newick utilities v1.6
|- nw_clade
|- nw_labels
|- nw_distance
awk v3.1.7
```
Basic usage: bash calculate_clade1_clade2_branch_len.bash file_of_newick_trees <br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### busco2alignment.py
Takes the full table output of busco runs and creates a concatenated fasta file for phylogenetic inference using the concatenation method. Only argument is a configuration file specified with -c.
Configuration file should specify the path to [programs] mafft, trimAl, [input_files] a single column file with the names of the busco output directories, a single column file of the fasta files, [output_files] concatenation fasta file name, partition file output name, and [parameters] taxon occupancy with a value between 0 and 1. Example format is the following: <br /> 
[programs] <br />
mafft: pathway to mafft <br />
trimAl: pathway to trimAl <br /> 
<br />
[input_files] <br />
busco_out_list: busco_dirs.list <br />
fasta_files_list: fasta_files.list <br />
<br /> 
[output_files] <br />
concat_fasta: concat.fa <br />
partition_file: USCO_partition.txt <br /> 
<br />
[parameters] <br />
occupancy: value between 0 and 1 <br />
<br />
Exemplary template file can be printed out using the -t option
NOTE: busco_out and fasta_files should have files in the same order <br />
and the fasta file header names should be formatted in the following manner:<br /> 
\>indivID|1<br /> 
...<br /> 
\>indivID|2<br /> 
...<br /> 
where '|' is used in every gene name to denote specific genes but text prior to denotes a unique identifier for all genes from that fasta file.
For detailed information about usage, use the -h argument <br />
```
python3
|- sys
|- Bio
   |- SeqIO
   |- SeqRecord
      |- SeqRecord
   |- Alphabet
      |- IUPAC
   |- Seq
      |- Seq
|- getopt
|- os.path
|- os
|- re
|- subprocess
|- re
|- configparser
```
Basic usage: python busco2alignment.py -c config.busco2alignment<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

This script is also broken up into two parts busco2uscofa.py and trimal2concat.py. <br />
Input and usage are the same but alignment and trimming must be done without these scripts. For trimal2concat.py, trimal files are assumed to have the following naming scheme buscoID.fa.mafft.trimal.

### busco_occupancy.py
Takes the full table output of busco runs and creates an occupancy matrix where columns are taxa and rows are busco IDs. The input file is a list of busco output dirs which is the same as "busco_out_list" in busco2alignment.py. Specifically, it is a single column file with busco output dir file names with one 
```
python3
|- sys
|- getopt
|- os.path
|- os
|- re
```
Basic usage: python busco_occupancy.py -b busco_output_dirs -o output_file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### threadNuclProt.py
Takes a protein sequence alignment and the corresponding nucleotide fasta file and threads the corresponding codons on top of the protein alignment. The order of genes in the nucleotide and protein fasta file must be the same. To facilitate use of this script for downstream programs, a boolean for whether the stop codon should be maintained or replaced with gaps is specified with the '-s' parameter representing a true or false for if the stop codon should be kept. Additionally, ambiguous amino acids 'X' are replaced with gaps. 
```
python3
|- sys
|- getopt
|- os.path
|- os
|- Bio
   |- SeqIO
```
Basic usage: python threadNuclProt.py -p protein.MSA.fasta -n nucleotide.fasta -s T/F<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### genomic_BUSCOs2uscofa.py
Creates BUSCO fasta files across given taxa using BUSCO outputs created using the '-m genome' parameter - that is, BUSCOs/genome completeness was determined from the genome fasta file and NOT the proteome or cds file. Requires a single column file with the names of the busco output directories and a list of the taxa names. It is required that the order of the busco output directories and the taxa names are in the same order. Taxon occupancy cut-offs for each BUSCO can be specified using a decimal value between 0 and 1. This script can create BUSCO fasta files for either nucleotide or protein fasta files.
```
python3
|- sys
|- getopt
|- os.path
|- os
|- Bio
   |- SeqIO
|- datetime
   |- datetime
```
Basic usage: python genomic_BUSCOs2uscofa.py -b busco_dirs.list -o 0-1 -t taxon_names.list -c nucl/prot<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### create_concat_matrix.py
Create a concatenated matrix for phylogenomic inference. The input files and parameters include a list of alignment files to concatenate, a list of taxa to include, whether the sequences are proteins or nucleotides, and a prefix for output files. Output files include a fasta file of concatenated sequence with '.fa' appended to the end, a RAxML style partition file, and a file that summarizes the occupancy of each gene from each alignment.
```
import sys
import getopt
import os.path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from datetime import datetime
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- SeqIO
   |- SeqRecord
      |- SeqRecord
   |- Alphabet
      |- IUPAC
   |- Seq
      |- Seq
|- datetime
   |- datetime
```
Basic usage: python create_concat_matrix.py -a alignment.list -t taxa.list -p prefix -c nucl/prot<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### variable_sites_of_alignment.py
Determines the number and percentage of variable sites in a fasta alignment.
```
python3
|- sys
|- getopt
|- os.path
|- re
|- Bio
   |- SeqIO
   |- Seq
      |- Seq
```
Basic usage: python genomic_BUSCOs2uscofa.py -b busco_dirs.list -o 0-1 -t taxon_names.list -c nucl/prot<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### Identify_terminal_branches_xx_longer_than_the_median_branch_length.py
Determines the branch length of each branch (internal and terminal) in a newick tree file. Then it determines the median branch length and reports any terminal branches that are equal to or greater than 'xx' the median terminal branch length. 'xx' is defined by the -n parameter. This script was designed specifically to identify potentially spurious sequences and subsequently remove them from a single gene tree.
```
python3
|- sys
|- getopt
|- os.path
|- statistics
|- Bio
   |- Phylo
```
Basic usage: python Identify_terminal_branches_xx_longer_than_the_median_branch_length.py -t tree.file -n 20<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### bootstrap_statistics.py
Prints out various statistics regarding the bootstrap support of a phylogenetic tree. Statistics printed include the mean, median, standard deviation, variance, and the 25th and 75th percentile.
```
python3
|- sys
|- getopt
|- os.path
|- statistics
|- numpy
|- Bio
   |- Phylo
```
Basic usage: python bootstrap_statistics.py -t tree.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

## Authors

* **Jacob Steenwyk** - [Github page](https://jlsteenwyk.github.io/)
* The online community of bioinformaticians. Other/original authors are listed per script.

## Acknowledgments

* [Rokas lab personnel](https://as.vanderbilt.edu/rokaslab/people/)
* [ACCRE](http://www.accre.vanderbilt.edu/)
* John Soghigian - [Github page](http://www.vector-eco-evo.com/)
* Xingxing Shen - [Github page](https://xingxingshen.github.io/)