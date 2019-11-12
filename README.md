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
Basic usage: python variable_sites_of_alignment.py -i alignment.fasta -c nucl/prot<br />
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

### shuffle_each_site_in_mfasta.py
Shuffles each site in an aligned fasta file such that each column is randomly shuffled. The proportion of each nucleotide or protein is maintained for each column.
```
python3
|- sys
|- getopt
|- os.path
|- re
|- random
|- Bio
   |- Phylo
   |- SeqRecord
      |- SeqRecord
   |- Seq
      |- Seq
```
Basic usage: python shuffle_each_site_in_mfasta.py -i fasta.file -t prot/nucl -s seed<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### shuffle_taxa_in_a_tree.py
Shuffles taxa in a newick tree. This script also strips out bootstrap values.
```
python3
|- sys
|- getopt
|- os.path
|- random
|- Bio
   |- Phylo
```
Basic usage: python shuffle_each_site_in_mfasta.py -t tree.file -s seed<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### print_tree.py
Prints a phylogenetic newick tree in ASCII format.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- Phylo
```
Basic usage: python print_tree.py -i tree.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### create_labels_for_internodes_as_support.py
Creates internode labels from 1 to n where n is the number of internodes in the tree. Internode labels are reported as confidence scores that can be easily viewed in FigTree of other similar software. The only argument is -t to specify a newick tree file. The output will be a tree with the same file name and ".internodeLabels.tree" appended to the end.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- Phylo
      |- BaseTree
         |- Clade
```
Basic usage: python create_labels_for_internodes_as_support.py -t tree.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### report_ic_as_internode_support.py
Internode certainty scores as outputted from RAxML are modified to be reported branch support labels. The input newick tree is specified with the -i parameter. Internode certainty values can then be easily viewed in FigTree of other similar software.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- Phylo
      |- BaseTree
         |- Clade
|- re
```
Basic usage: python report_ic_as_internode_support.py -i tree.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### retrieve_internode_certainty_values.py
Extracting internode certainty values may be of interest to plot the distribution of internode certainty values for various trees. This script extract the internode certainty values of a single tree produced from a RAxML output. The phylogenetic tree in newick format is specified with the -t parameter.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- Phylo
|- re
```
Basic usage: python retrieve_internode_certainty_values.py -t tree.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### root_on_clade.py
Roots tree on specified clade of individuals. The individuals to root the tree on are specified in a separate file that contains a single column file list of taxa. This file is specified using the -o parameter. The input tree is specified with the -i parameter. The output is a newick tree file with the same name as the original tree file with .rooted.tre appended to the end of the file name.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- Phylo
      |- BaseTree
         |- Clade
```
Basic usage: python root_on_clade.py -t tree.file -o individuals.list<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### prune_list_of_taxa.py
Prunes list of a taxa from a newick tree file. Arguments include -i to specific a newick tree and -l to specify a file with a single column list of taxa to prune. Output is a tree with the same name as the input tree with .pruned appended to the end of the name.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- Phylo
      |- BaseTree
         |- TreeMixin
```
Basic usage: python prune_list_of_taxa.py -i tree.file -l individuals.list<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### remove_distances_from_tree.py
Removes distances from a tree in newick format. This can be used if you are only interested in the topology of the tree or want to reestimate branch lengths after making a consensus tree. The input newick tree file is specified with the -t parameter.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- Phylo
      |- BaseTree
         |- Clade
```
Basic usage: python remove_distances_from_tree.py -t tree.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### DVMC_degree_of_violation_of_a_molecular_clock.py
Calculates the degree of violation of a molecular clock (DVMC) for a tree. The tree, specified with the -t parameter, is first rooted on the outgroup taxa and then those outgroup taxa are pruned from the tree - outgroup taxa are specified in a single column separate file using the -o parameter. After that, DVMC is calculated. The outgroup taxa should be the taxa that are not used in the timetree but present in the phylogenetic tree. This is calculated according to the following paper: http://www.pnas.org/content/114/35/E7282
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- Phylo
      |- BaseTree
         |- TreeMixin
|- math
```
Basic usage: python DVMC_degree_of_violation_of_a_molecular_clock.py -t tree.file -o outgroup.list <br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### create_taxon_bootstrap_alignments.py
This script is designed to make subsets of a given fasta file. Subsets are created based on taxa and not the alignment such that a fraction of the total taxa represented in the input fasta file are used to create the resulting fasta files. This script samples without replacement. The input fasta file is specified using the -i argument, the number of taxa to subset is specified using the -l argument, the number of replicates to perform is specified using the -r argument, the seed is set using the -s argument, and the prefix of the output files is specified using the -p argument.
```
python3
|- sys
   |- stdout
|- getopt
|- os.path
|- Bio
   |- SeqIO
|- random
|- itertools
```
Basic usage: python create_taxon_bootstrap_alignments.py -i input.fa -l 50 -r 100 -s 1 -p subset <br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### RCV.py
Calculate relative nucleotide composition variability (RCV) as defined by Phillips and Penny 2003, specifically in the manuscript titled 'The root of the mammalian tree inferred from whole mitochondrial genomes. Mol Phylogenet Evol. 28:171–185. https://www.ncbi.nlm.nih.gov/pubmed/12878457. The only file is an input fasta file which, is specified using the -i argument.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- SeqIO
|- Bio.Seq
   |- Seq
```
Basic usage: python RCV.py -i fasta.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### compare_lengths_of_two_alignments.py
Compare the alignment length between two alignments. This script was designed for phylogenomics projects. More specifically, it is good to filter out genes who's trimmed and aligned length is less than 50 percent the aligned gene's length. The output is the length of the -t alignment, the length of the -m alignment followed by the percentage of the -t alignment divided by the -m alignment. The logic of the -t and -m denomination is the -m represents a mafft alignment and -t represents a trimAl alignment.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- SeqIO
```
Basic usage: python compare_lengths_of_two_alignments.py -m alignment.file -t alignment.file<br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### total_tree_length.py
Calculates tree length of a phylogenetic tree in newick format.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- Phylo
```
Basic usage: python compare_lengths_of_two_alignments.py -t newick.tre <br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### compare2trees.py
This script compares two trees and reports the differences between them.
Input trees must be rooted and contain the same set of taxa. There is no
test within the script for these parameters so ensure the trees have been
properly processed before using this script. Note, output files will be
rewritten from previous runs of the script.
```
python3
|- sys
|- getopt
|- os.path
|- Bio
   |- Phylo
|- datetime
   |- datetime
```
Basic usage: python compare2trees.py -a tree1 -b tree2 <br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)

### prepare_tags_for_BayesTraits.py
Creates tags for BayesTraits run/parameters file by looping through all internal nodes and identifying children tips from there. Nodes are labeled according to how they are encoutnered in the tree file. The tree file should be a nexus file and be the same one as is used as input for BayesTraits.
```
python3
|- sys
|- getopt
|- os.path
|- Bio.Phylo.BaseTree
   |- Clade
|- pprint
   |- pprint
```
Basic usage: python prepare_tags_for_BayesTraits.py -t tree <br />
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)


## Authors

* **Jacob Steenwyk** - [Github page](https://jlsteenwyk.github.io/)
* The online community of bioinformaticians. Other/original authors are listed per script.

## Acknowledgments

* [Rokas lab personnel](https://as.vanderbilt.edu/rokaslab/people/)
* [ACCRE](http://www.accre.vanderbilt.edu/)
* John Soghigian - [Github page](http://www.vector-eco-evo.com/)
* Xingxing Shen - [Github page](https://xingxingshen.github.io/)