# -*- coding: utf-8 -*-
"""
get_clustered.py

Created on Fri May 13 14:06:31 2022

@author: Elena Aramendía

This script was written to parse SILVAngs output files; this software clusters 
the query sequences into OTUs and outputs a file with the sequences included 
in each OTU. This script reads that cluster mapping file to get all the 
sequences included in each OTU. 

INPUT:
    - Cluster mapping file, includes the header of the reference sequence for
    the OTU and of the sequences clustered in each of them.
    - File with header of reference OTU sequences of interest, the output will contain 
    this sequences and all of the sequences included (clustered) in them.

OUTPUT
    - Text file with all the headers of the reference sequences provided and
    the sequences clustered together with that reference. 
    Name of the output file can be specified, if not the the default is
    output_clustered.txt

EXAMPLE USAGE:
    python get_clustered.py cluster_map.clstr seqs.txt [output_clustered.txt]
"""

#%% Input and output files

# Modules
import sys

## Read input files

# Clusters
cluster_map = sys.argv[1]
# Brev ref SRAs
seqs = sys.argv[2]

## Ouput
all_seqs = sys.argv[3]
# Default output file
if not all_seqs:
    all_seqs = "output_clustered.txt"
    
## Errors
if not cluster_map and seqs:
    raise(Exception("Please specify two input files, first the cluster map and second a file including names of the sequences of interest."))



#%% Code
# Dictionary with seqs in each cluster
# Key is the ref sequence for the OTU
otus = {}

for line in cluster_map:
    if line.startswith(">"):
        sra1 = line.split("x")[0].strip(">")
        if sra1 not in otus:
            otus[sra1] = list()
    else:
        line = line.split(",")[1]
        sra2 = line.split("x")[0].strip()
        sra2 = sra2.strip(">")
        
        if sra1 != sra2:
            otus[sra1].append(sra2)
            
# Get SRAs of interests -> Brev sequences
for line in seqs:
    line = line.strip()
    if line in otus:
        print(line, file=all_seqs)
        for i in otus[line]:
            print(i, file=all_seqs)