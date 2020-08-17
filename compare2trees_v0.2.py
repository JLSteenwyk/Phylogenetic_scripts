#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import Phylo
from datetime import datetime

def compare_trees(
	tree1,
	tree2,
	tree1Name,
	tree2Name,
	startTime
	):
    """
    compares the two trees and identifies differences

    Parameters
    ----------
    argv: tree1
        reference tree
    argv: tree2
        second tree
    argv: tree1Name
        reference tree name
    argv: tree2Name
        second tree name
    argv: startTime
        time script was started
    """

    # initialize conflict counter
    conflictCNT = 0

    # initialize list to keep track of conflicting topologies
    taxa_conflicting_between_two_trees_tree1 = []
    taxa_conflicting_between_two_trees_tree2 = []

    # Initialize script to create plots of all conflicting topologies
    fileNameALL='create_plots_of_conflicting_topologies.R'
    with open(fileNameALL, "w+") as f:
        # tested using R version 3.5.2
        # tested using ggtree v1.14.6
        # tested using cowplot v0.9.4
        # tested using emojifont v0.5.2
        f.write("library(ggtree)\n")
        f.write("library(cowplot)\n")
        f.write("library(emojifont)\n")
        f.write("library(ggplot2)\n")

    # Initialize text file that describes all conflicts
    fileNameALL='Table_summarizing_conflicts.txt'
    with open(fileNameALL, "w+") as f:
        f.write("Conflict #\tTaxa in "+tree1Name+"\tTaxa in "+tree2Name+"\tTaxa in collapsed node\n")

    ## loop through tree1 and find similar clade in tree2
    for clade1 in tree1.get_nonterminals():
        # initialize and populate a list of tip names in tree1
        tipNames1 = []
        for leaf in clade1.get_terminals():
            tipNames1.append(leaf.name)
        # get common ancestor of tree1 tip names in tree2
        clade2 = tree2.common_ancestor(tipNames1)
        # initialize and populate a list of tip names in tree2
        tipNames2 = []
        for leaf in clade2.get_terminals():
            tipNames2.append(leaf.name)

        # compare the list of tip names
        if set(tipNames1) == set(tipNames2):
            continue
        else:
            try:
                # get parent clades
                parent1=tree1.get_path(clade1)[-2]
                parent2=tree2.get_path(clade2)[-2]
            except IndexError:
                print("There is an error at the root")
                print("Are both phylogenies rooted on the same taxa?")
                print("Exiting now...")
                sys.exit()

            # determine monophyletic children clades
            for child1 in parent1.get_nonterminals()[0]:
                subtreeLeaves1 = []
                subtreeLeaves2 = []
                for leaf1 in child1.get_terminals():
                    subtreeLeaves1.append(leaf1.name)
                child2 = tree2.common_ancestor(subtreeLeaves1)
                for leaf2 in child2.get_terminals():
                    subtreeLeaves2.append(leaf2.name)
                #Phylo.draw_ascii(child1)
                if (set(subtreeLeaves1)==set(subtreeLeaves2)) and (len(subtreeLeaves1) >= 5):
                    break
                    
            
            # # create file name to with ascii trees of conflicting topologies
            # fileName='conflict_'+str(conflictCNT)+'.txt'
            # with open(fileName, "w+") as f:
            #     f.write("-------------------\n")
            #     f.write("| Tree 1 Topology |\n")
            #     f.write("-------------------\n")
            #     Phylo.draw_ascii(parent1, file=f)
            #     f.write("\n\n")
            #     f.write("-------------------\n")
            #     f.write("| Tree 2 Topology |\n")
            #     f.write("-------------------\n")
            #     Phylo.draw_ascii(parent2, file=f)
            #     f.close()
            
            # # create file with newick trees of conflicting topologies
            # fileNameT1='conflict_'+str(conflictCNT)+'_'+tree1Name
            # with open(fileNameT1, "w+") as f:
            #     Phylo.write(parent1, f, 'newick')
            # fileNameT2='conflict_'+str(conflictCNT)+'_'+tree2Name
            # with open(fileNameT2, "w+") as f:
            #     Phylo.write(parent2, f, 'newick')
            
            # # create an R script to plot all conflicting topologies
            # fileNameALL='create_plots_of_conflicting_topologies.R'
            # with open(fileNameALL, "a+") as f:

            #     # write lines for t1
            #     f.write("\n\n#code for t1\n")
            #     line='tre1<-read.tree("'+fileNameT1+'")\n'
            #     f.write(line)
            #     line="trePlot1<-ggtree(tre1, branch.length='none') + geom_tiplab() + geom_point2(aes(subset=(node==MRCA(tre1, tip=c("+str(tipNames1).replace("[","").replace("]","")+")))), size=4, color='firebrick2')"
            #     f.write(line)
            #     line="+ labs(title='"+tree1Name+"', caption='Red internodes show conflict')\n"
            #     f.write(line)
            #     if len(subtreeLeaves1) >= 5:
            #         line="internodeLABEL<-MRCA(tre1, tip=c("+str(subtreeLeaves1).replace("[","").replace("]","")+"))\n"
            #         f.write(line)
            #         line="trePlot1<-collapse(trePlot1, internodeLABEL, mode= 'none') + geom_cladelabel(node=internodeLABEL, 'Collapsed "+str(len(subtreeLeaves1))+" taxa', geom='label', offset=.25)"
            #         f.write(line)
            #         line="+ geom_cladelabel(node=internodeLABEL, intToUtf8(9664), family = 'OpenSansEmoji', fontsize=10, offset=-.08)\n"
            #         f.write(line)

            #     # write lines for t2
            #     f.write("\n\n#code for t2\n")
            #     line='tre2<-read.tree("'+fileNameT2+'")\n'
            #     f.write(line)
            #     line="trePlot2<-ggtree(tre2, branch.length='none') + geom_tiplab() + geom_point2(aes(subset=(node==MRCA(tre2, tip=c("+str(tipNames1).replace("[","").replace("]","")+")))), size=4, color='firebrick2')"
            #     f.write(line)
            #     line="+ labs(title='"+tree2Name+"', caption='Red internodes are parent to the conflict in t1')\n"
            #     f.write(line)
            #     if len(subtreeLeaves2) >= 5:
            #         line="internodeLABEL<-MRCA(tre2, tip=c("+str(subtreeLeaves2).replace("[","").replace("]","")+"))\n"
            #         f.write(line)
            #         line="trePlot2<-collapse(trePlot2, internodeLABEL, mode= 'none') + geom_cladelabel(node=internodeLABEL, 'Collapsed "+str(len(subtreeLeaves1))+" taxa', geom='label', offset=.25)"
            #         f.write(line)
            #         line="+ geom_cladelabel(node=internodeLABEL, intToUtf8(9664), family = 'OpenSansEmoji', fontsize=10, offset=-.08)\n"
            #         f.write(line)

            #     # create plot w/ trees side by side
            #     f.write("\n\n#create plot w/ trees side by side\n")
            #     line='pdf("conflict_'+str(conflictCNT)+'.pdf")\n'
            #     f.write(line)
            #     f.write("plot_grid(trePlot1, trePlot2, ncol=2)\n")
            #     f.write("dev.off()\n")

            # write to summary table
            fileName='Table_summarizing_conflicts.txt'
            with open(fileName, "a+") as f:
                #	Conflict #		Taxa in tree1Name		Taxa in tree2Name		Taxa in collapsed node
                if (len(subtreeLeaves1) >= 5):
                    f.write("conflict_"+str(conflictCNT)+"\t"+str(tipNames1).replace("[","").replace("]","").replace("'","")+"\t")
                    f.write(str(tipNames2).replace("[","").replace("]","").replace("'","")+"\t")
                    f.write(str(subtreeLeaves2).replace("[","").replace("]","").replace("'","")+"\n")
                else:
                    f.write("conflict_"+str(conflictCNT)+"\t"+str(tipNames1).replace("[","").replace("]","").replace("'","")+"\t")
                    f.write(str(tipNames2).replace("[","").replace("]","").replace("'","")+"\tNA\n")


        # append to taxa_conflicting_between_two_tree1/2
        conflicting2LeafLabels1=[]
        conflicting2LeafLabels2=[]

        for leaf in clade1.get_terminals():
            conflicting2LeafLabels1.append(leaf.name)
        for leaf in clade2.get_terminals():
            conflicting2LeafLabels2.append(leaf.name)
        taxa_conflicting_between_two_trees_tree1.append(conflicting2LeafLabels1)
        taxa_conflicting_between_two_trees_tree2.append(conflicting2LeafLabels2)

        # add one to conflict counter
        conflictCNT+=1
    
    # create plot for full tree side by side
    fileNameALL='create_summary_plot_comparing_two_topologies.R'
    with open(fileNameALL, "w+") as f:
        # tested using R version 3.5.2
        # tested using ggtree v1.14.6
        # tested using cowplot v0.9.4
        # tested using emojifont v0.5.2
        f.write("library(phytools)\n")
        f.write("library(cowplot)\n")
        f.write("library(ape)\n")
        
        # read in trees
        f.write("\n\n# Read in trees\n")
        line='tre1<-read.tree("'+tree1Name+'")\n'
        f.write(line)
        line='tre2<-read.tree("'+tree2Name+'")\n'
        f.write(line)

        # create function to resize fonts automatically
        f.write("foo<-function(tree,...){\n")
        f.write("fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)\n")
        f.write("plotTree(tree,fsize=fsize,lwd=1,...)\n")
        f.write("}\n")

        # create plot for t1 with geom_points
        f.write("\n\n# Code for plotting tree1\n")
        f.write("pdf('tree1_conflicts.pdf')\n")
        for ii in taxa_conflicting_between_two_trees_tree1:
            if taxa_conflicting_between_two_trees_tree1.index(ii) == 0:
                line="foo(tre1, edge.width=2)" + "\n" + "nodelabels(node=findMRCA(tre1, tips=c("+str(ii).replace("[","").replace("]","")+"), type = 'node'), frame='circle', cex=0.25, bg='red', text='')" + "\n"
                f.write(line)
            elif taxa_conflicting_between_two_trees_tree1.index(ii)+1 < len(taxa_conflicting_between_two_trees_tree1):
                line="nodelabels(node=findMRCA(tre1, tips=c("+str(ii).replace("[","").replace("]","")+"), type = 'node'), frame='circle', cex=0.25, bg='red', text='')" + "\n"
                f.write(line)
            elif taxa_conflicting_between_two_trees_tree1.index(ii)+1 == len(taxa_conflicting_between_two_trees_tree1):
                line="nodelabels(node=findMRCA(tre1, tips=c("+str(ii).replace("[","").replace("]","")+"), type = 'node'), frame='circle', cex=0.25, bg='red', text='')"
                f.write(line)
        f.write("\ndev.off()")

        # create plot for t1 with geom_points and no branch lengths
        f.write("\n\n# Code for plotting tree1\n")
        f.write("pdf('tree1_conflicts_no_bls.pdf')\n")
        f.write("tre1$edge.length<-NULL\n")
        for ii in taxa_conflicting_between_two_trees_tree1:
            if taxa_conflicting_between_two_trees_tree1.index(ii) == 0:
                line="foo(tre1, edge.width=2)" + "\n" + "nodelabels(node=findMRCA(tre1, tips=c("+str(ii).replace("[","").replace("]","")+"), type = 'node'), frame='circle', cex=0.25, bg='red', text='')" + "\n"
                f.write(line)
            elif taxa_conflicting_between_two_trees_tree1.index(ii)+1 < len(taxa_conflicting_between_two_trees_tree1):
                line="nodelabels(node=findMRCA(tre1, tips=c("+str(ii).replace("[","").replace("]","")+"), type = 'node'), frame='circle', cex=0.25, bg='red', text='')" + "\n"
                f.write(line)
            elif taxa_conflicting_between_two_trees_tree1.index(ii)+1 == len(taxa_conflicting_between_two_trees_tree1):
                line="nodelabels(node=findMRCA(tre1, tips=c("+str(ii).replace("[","").replace("]","")+"), type = 'node'), frame='circle', cex=0.25, bg='red', text='')"
                f.write(line)
        f.write("\ndev.off()")

        # # create plot for t1 with geom_points
        # f.write("\n\n# Code for plotting tree1\n")
        # f.write("pdf('tree2_conflicts.pdf')\n")
        # for ii in taxa_conflicting_between_two_trees_tree1:
        #     if taxa_conflicting_between_two_trees_tree2.index(ii) == 0:
        #         line="foo(tre2, edge.width=2)" + "\n" + "nodelabels(node=findMRCA(tre2, tips=c("+str(ii).replace("[","").replace("]","")+"), type = 'node'), frame='circle', cex=0.25, bg='red', text='')" + "\n"
        #         f.write(line)
        #     elif taxa_conflicting_between_two_trees_tree1.index(ii)+1 < len(taxa_conflicting_between_two_trees_tree1):
        #         line="nodelabels(node=findMRCA(tre2, tips=c("+str(ii).replace("[","").replace("]","")+"), type = 'node'), frame='circle', cex=0.25, bg='red', text='')" + "\n"
        #         f.write(line)
        #     elif taxa_conflicting_between_two_trees_tree1.index(ii)+1 == len(taxa_conflicting_between_two_trees_tree1):
        #         line="nodelabels(node=findMRCA(tre2, tips=c("+str(ii).replace("[","").replace("]","")+"), type = 'node'), frame='circle', cex=0.25, bg='red', text='')"
        #         f.write(line)
        # f.write("\ndev.off()")


    print("\n-----------------")
    print("| Compare2Trees |")
    print("-----------------\n")
    print("author: Jacob L. Steenwyk")
    print("citation: NA\n\n")
    print("Number of conflicting topologies identified:", str(conflictCNT),"\n")
    print("Plot individual conflicting topologies w/ create_plots_of_conflicting_topologies.R")
    print("Plot all conflicting topologies w/ create_summary_plot_comparing_two_topologies.R")
    print("Get a summary of each topolgoy w/ Table_summarizing_conflicts.txt\n")
    print("Complete!")
    print("Total time:", datetime.now() - startTime)
    print("")


def read_inputs(
    tree1,
    tree2,
    startTime
    ):
    """
    Read in input tree files

    Parameters
    ----------
    argv: tree1
        reference tree
    argv: tree2
        second tree
    argv: startTime
        time the script was started
    """

    tree1Name = tree1
    tree2Name = tree2
    tree1 = Phylo.read(tree1, 'newick')
    tree2 = Phylo.read(tree2, 'newick')

    # pass to compare_trees function
    compare_trees(
        tree1, 
        tree2,
        tree1Name,
        tree2Name,
        startTime
        )


def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree1 = ''
    tree2 = ''

    try:
        opts, args = getopt.getopt(argv, "ha:b:")
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
            print("\n-----------------")
            print("| Compare2Trees |")
            print("-----------------\n")
            print("author: Jacob L. Steenwyk")
            print("citation: NA\n")
            print("This script compares two trees and reports the differences between them.")
            print("Input trees must be rooted and contain the same set of taxa. There is no")
            print("test within the script for these parameters so ensure the trees have been")
            print("properly processed before using this script. Note, output files will be")
            print("rewritten from previous runs of the script.\n\n")
            print("Usage:")
            print("python compare2trees.py -a tree1 -b tree2\n")
            print("• tree1 is a newick tree file. This tree is the 'reference' tree")
            print("• tree2 is a newick tree file that will be compared to tree1\n\n")
            print("The output will be a series of files.")
            print("(1) A set of newick tree files with the topology in tree1 and tree2 for each")
            print("\tconflicting internode. They will be named conflict_x_n.tre")
            print("\twhere x is the numbered conflict starting from zero and")
            print("\tn is the name of the tree that the topology is observed in.\n")
            print("(2) A set of text files with the ascii tree with the topology observed in")
            print("\ttree1 and tree2. This is to get a quick idea of the topologies.\n")
            print("(3) An Rscript titled create_plots_of_conflicting_topologies.R which,")
            print("\tcreates individual plots of each conflicting topology.")
            print("\tThe topologies are displayed side by side. There will be")
            print("\tred dots at the conflicting internode of tree1 and a red")
            print("\tdot at the parent of the conflicting internode of tree2.")
            print("\tThis script relies on the following packages and has been")
            print("\ttested on the following versions of packages:")
            print("\t• R version 3.5.2")
            print("\t• ggtree v1.14.6")
            print("\t• cowplot v0.9.4")
            print("\t• emojifont v0.5.2\n")
            print("(4) An Rscript titled create_summary_plot_comparing_two_topologies.R which,")
            print("\tcreates a plot of the full topology from tree1 and tree2")
            print("\twith red dots at conflicting internodes of tree1 and the")
            print("\tparent of the conflicting internode for tree2.\n")
            print("(5) Lastly, a file titled Table_summarizing_conflicts.txt which, summarizes")
            print("\teach conflicting internode. Col 1 is the number of the conflict,")
            print("\tis the taxa in tree1, Col 3 is the taxa in tree2, and Col 4")
            print("\tis the taxa in the collapsed internode. If there is no collapsed")
            print("\tinternode, there will be an 'NA.'")

            print("\n")
            sys.exit()
        elif opt == '-a':
            if arg:
                tree1 = arg
            else:
                # error message
                print("\n\nThe specified file of arg -a does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-b':
            if arg:
                tree2 = arg
                startTime = datetime.now()
            else:
                # error message
                print("\n\nThe specified input of arg -b must be an integer.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()


    # pass to read inputs function
    read_inputs(
        tree1,
        tree2,
        startTime
        )

if __name__ == '__main__':
    main(sys.argv[1:])
