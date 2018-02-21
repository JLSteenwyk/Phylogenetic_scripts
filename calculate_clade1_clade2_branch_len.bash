#!/bin/bash
# usage: bash myscript.bash file_of_newick_trees
#
# this script will determine if taxa from Clade1 and Clade2 are monophyletic
# if monophyletic, it will calculate the distance of the branches leading up to Clade1 and Clade2
# given a file that has one tree per line in newick format (ie., a file with n number of lines 
# where each line is a newick tree). This is the same input file that is used for ASTRAL
#
# With each usage change variables
# Clade1, Clade2, and All_Clade 
# and replace with taxa of interest
#
# output is:
# col1: distance for Clade1 internode branch
# col2: distance for Clade2 internode branch
#

INfile=$1
Clade1=$(echo "Hjakobsenii_35692-160519 Hsingularis_35696-160519 Khatyaiensis_35684-160519 Hclermontiae_34963-160519 Hmeyeri_34960-160519 Huvarum_DSM2768 Huvarum_34956-170307 Huvarum_AWRI3580 Huvarum_34-9 Hthailandica_35697-160519 Hnectarophila_35693-160613 Hguilliermondii Hguilliermondii_34962-170307 Hlachancei_34961-170307 Hpseudoguilliermondii_35695-160519 Hopuntiae Hvalbyensis Hlibkind_35685-160519")
Clade2=$(echo "Hoccidentalis_var_citrica_35694-160519 Hoccidentalis_var_occidentalis_34959-160519 Han_19 Hosmophila Hosmophila_34957-160519 Hvinae Hvineae_26147-160519")
All_Clade=$(echo "Hoccidentalis_var_citrica_35694-160519 Hoccidentalis_var_occidentalis_34959-160519 Han_19 Hosmophila Hosmophila_34957-160519 Hvinae Hvineae_26147-160519 Hjakobsenii_35692-160519 Hsingularis_35696-160519 Khatyaiensis_35684-160519 Hclermontiae_34963-160519 Hmeyeri_34960-160519 Huvarum_DSM2768 Huvarum_34956-170307 Huvarum_AWRI3580 Huvarum_34-9 Hthailandica_35697-160519 Hnectarophila_35693-160613 Hguilliermondii Hguilliermondii_34962-170307 Hlachancei_34961-170307 Hpseudoguilliermondii_35695-160519 Hopuntiae Hvalbyensis Hlibkind_35685-160519")

while read line
	do 
	
	# determine monophyly of clade 1 (FE) 
	FE_LABELS_Num=$(echo -n "$line " | \
		nw_clade -m - $Clade1 | \
		nw_labels - | wc -l)
	
	# determine monophyly of clade 2 (SE)	
	SE_LABELS_Num=$(echo -n "$line " | \
		nw_clade -m - $Clade2 | \
		nw_labels - | wc -l)

	# test for monophyly
	if [ $FE_LABELS_Num == "18" ] && [ $SE_LABELS_Num == "7" ]; then

		### CALCULATE LENGTH OF BRANCH LEADING TO CLADE 1 (FE)
		# get name of shortest branch to LCA among clade 1 (FE)
		SHORT_LEAF_2_LCA_LABEL=$(echo -n "$line" | nw_clade -m - \
			$Clade1 | nw_distance -n - | sort -nk2,2 | head -n 1 | awk '{print $1}')

		# get length of shortest branch to LCA among 1 and 2 clades (FE and SE)
		SHORT_LEAF_2_LCA_DIST=$(echo -n "$line" | nw_clade -m - \
			$All_Clade | nw_distance -n - | grep -w "$SHORT_LEAF_2_LCA_LABEL" | awk '{print $2}')

		# get length of shortest branch to LCA among clade 1 (FE)
		SHORT_LEAF_2_PARENT_DIST=$(echo -n "$line" | nw_clade -m - \
			$Clade1 | nw_distance -n - | grep -w $SHORT_LEAF_2_LCA_LABEL | awk '{print $2}')

		# calculate difference between LCA of clade 1 and 2 (FE and SE) and parent of shortest branch in clade 1 (FE)
		# this therefore calculates the branch length of the internode branch leading up to clade 1
		CLADE1_LEN=$(echo -e "$SHORT_LEAF_2_LCA_DIST\t$SHORT_LEAF_2_PARENT_DIST" | awk '{print $1-$2}')
		echo -n "$CLADE1_LEN"
		###

		### CALCULATE LENGTH OF BRANCH LEADING TO CLADE 2 (SE)
		SHORT_LEAF_2_LCA_LABEL=$(echo -n "$line" | nw_clade -m - \
			$Clade2 | nw_distance -n - | sort -nk2,2 | head -n 1 | awk '{print $1}')						
		
		# get length of shortest branch to LCA among 1 and 2 clades (FE and SE)
		SHORT_LEAF_2_LCA_DIST=$(echo -n "$line" | nw_clade -m - \
			$All_Clade | nw_distance -n - | grep -w "$SHORT_LEAF_2_LCA_LABEL" | awk '{print $2}')

		# get length of shortest branch to LCA among clade 2 (SE)		
		SHORT_LEAF_2_PARENT_DIST=$(echo -n "$line" | nw_clade -m - \
			$Clade2 | nw_distance -n - | grep -w $SHORT_LEAF_2_LCA_LABEL | awk '{print $2}')		
		
		# calculate difference between LCA of clade 1 and 2 (FE and SE) and parent of shortest branch in clade 2 (SE)
		# this therefore calculates the	branch length of the internode branch leading up to clade 2
		echo -e "$SHORT_LEAF_2_LCA_DIST\t$SHORT_LEAF_2_PARENT_DIST" | awk '{print "\t"$1-$2}'
		###		
	else
		continue
	fi

done <$INfile
