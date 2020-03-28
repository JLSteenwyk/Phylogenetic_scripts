#!/usr/bin/env Rscript

# The Robinson-Foulds distance is the number of splits that are present 
# in one tree but not in the other one, and vice versa. Since unrooted 
# n-taxa trees have a maximum of n - 3 inner branches, the maximal 
# Robinson-Foulds distance is 2(n - 3). Normalized Robinson-Foulds 
# distance is the RF divided by 2(n - 3). This yields a value between 
# 0% and 100%, which can be interpreted as the percentage of false or 
# missing splits in the inferred tree compared to the true tree. 

args <- commandArgs(trailingOnly=TRUE)

tre1<-ape::read.tree(args[1])
tre2<-ape::read.tree(args[2])

dist<-phangorn::RF.dist(tre1, tre2, normalize = TRUE)

cat(dist,"\n")