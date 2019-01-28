#!/usr/bin/env python

import sys

"""
concatenates sequences from a list of fasta files 
assuming 100% taxon occupancy 
"""

fastaList  = [line.rstrip('\n') for line in open(sys.argv[1])]
fasta      = {}
partition  = []

for file in fastaList:
    tempPartition = []
    tempPartition.append((file.split('.',1)[0]))
    with open(file) as fastaFile:
        for line in fastaFile:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                length = int(line.split(' ')[1])
                line = line.split(' ', 1)[0]
                active_sequence_name = line[1:]
                if active_sequence_name not in fasta:
                    fasta[active_sequence_name] = []
                continue
            sequence = line
            fasta[active_sequence_name].append(sequence)
    tempPartition.append(length)
    partition.append(tempPartition)

open('partition.file', "w")

with open('partition.file', 'a') as f:
    start=1
    for ele in partition:
        end=start+ele[1]-1
        entry="LG+G4, "+ele[0]+" = "+str(start)+"-"+str(end)+"\n"
        f.write(entry)
        start+=ele[1]

open('concatenation.fasta', "w")

with open('concatenation.fasta', 'a') as f:

    fasta_dict = dict()

    for k,v in fasta.items():
        if k not in fasta_dict:fasta_dict[k] = []
        fasta_dict[k].append(v)

    for k, v in fasta_dict.items():
        entry0='>'+str(k)+"\n"
        f.write(entry0)
        for elel in v:
           entry1=''.join(elel)+"\n"
           f.write(entry1)
