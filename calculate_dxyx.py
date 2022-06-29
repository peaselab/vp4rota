#!/usr/bin/env python3
import sys
from itertools import combinations
from time import time
time0 = time()
# Calculates dxy for multiple aligned fasta files.
# Each fasta file contains data for a single locus (aligned sequences for both ecotypes from multiple locations).
# For each locus, loops through all locations and calculates dxy between ecotypes.
# Version 0.2

def parse_fasta(infile):
    name, seq = None, []
    for line in infile: # for loop
        line = line.rstrip() # strip the white spaces on the right
        if line.startswith(">"): 
            if name: yield (name, ''.join(seq)) # generator function, append None, seq 
            name, seq = line, [] # name = line, add sequence into list
        else:
            seq.append(line) # if there is no name then you append the seq to the seq list. 
    if name: yield (name, ''.join(seq))


def dx(pop1, mode="nucleotide"):
    valid_chars = 'ATGCU'
    if mode == 'amino':
        valid_chars = 'ACDEFGHIKLMNPQRSTVWY'
    sum_of_pairwise_distances = 0.0 # running total of pairwise distances
    total_comparisons = 0 # number of sequences
    for i, seq1 in enumerate(pop1):
        seq1 = seq1.upper()
        seq1_sites = [i for i in range(len(seq1)) if seq1[i] in valid_chars]
        for seq2 in pop1[i + 1:]:
            seq2 = seq2.upper()
            # valid_sites = [i for i in seq1_sites if seq2[i] in valid_chars]
            n_diff = [seq1[i] != seq2[i] for i in seq1_sites if seq2[i] in valid_chars]
            sum_of_pairwise_distances += float(sum(n_diff) / len(n_diff)) if n_diff else 0
            total_comparisons += 1 if n_diff else 0
    if total_comparisons == 0:
        return "na", 'na'          
  
    return sum_of_pairwise_distances, total_comparisons

# function to calculate dxy, adjusted from script by Simon Martin
def dxy(pop1,pop2, mode='nucleotide'):
    valid_chars = 'ATGCU'
    if mode == 'amino':
        valid_chars = 'ACDEFGHIKLMNPQRSTVWY'
    sum_of_pairwise_distances = 0.0 # running total of pairwise distances
    total_comparisons = 0 # number of sequences
    for seq1 in pop1:
        seq1 = seq1.upper()
        seq1_sites = [i for i in range(len(seq1)) if seq1[i] in valid_chars]
        for seq2 in pop2:
            seq2 = seq2.upper()
            # valid_sites = [i for i in seq1_sites if seq2[i] in valid_chars]
            n_diff = [seq1[i] != seq2[i] for i in seq1_sites if seq2[i] in valid_chars]
            sum_of_pairwise_distances += float(sum(n_diff) / len(n_diff)) if n_diff else 0
            total_comparisons += 1 if n_diff else 0
    if total_comparisons == 0:
        return "na", 'na'         
    return sum_of_pairwise_distances, total_comparisons

# generate results file
outFileName = sys.argv[2]
outFile = open(outFileName, "w")
outFile.write("\t".join([
    "measure",
    "geno1",
    "geno2",
    "dxy",
    "dx1", 
    "dxyx1", 
    "dx2", 
    "dxyx2", 
    "n1", 
    "n2",
    'ncomp'])+
     "\n")

# make list of populations
# genotypes = ["P4", "P6", "P7", "P8", "P13", "P15"]
genotypes =[]
dx_values = {}
fafile = open(sys.argv[1], "r")
for label, seq in parse_fasta(fafile): 
    genotype = label.split("/")[6]
    genotypes.append(genotype)
fafile.close()
genotypes = list(set(genotypes))
print(genotypes)

for i, geno1 in enumerate(genotypes):
    pop1 = []
    fafile = open(sys.argv[1], "r")
    for label, seq in parse_fasta(fafile): 
        genotype = label.split("/")[6]
        if genotype == geno1:
            pop1.append(seq)
    if pop1: 
        dx_total, dx_comparisons = dx(pop1)
        print('dx', geno1, len(pop1), dx_comparisons, 
              dx_total / dx_comparisons if dx_comparisons !='na' else 'na')
        if dx_total != 'na':
            dx_values[geno1] = dx_total / dx_comparisons

for i, geno1 in enumerate(genotypes):
    for geno2 in genotypes[i+1:]:
        pop1 = []
        pop2 = []
        fafile = open(sys.argv[1], "r")
        for label, seq in parse_fasta(fafile): 
            genotype = label.split("/")[6]
            if genotype == geno1:
                pop1.append(seq)
            if genotype == geno2:
                pop2.append(seq)
        if pop1 and pop2: 
            dxy_total, dxy_comparisons = dxy(pop1, pop2)
            dxy_value = 'na' if dxy_total == 'na' else dxy_total / dxy_comparisons
            print('dxy', geno1, geno2, len(pop1), len(pop2), dxy_comparisons, 
                  dxy_value)
            outFile.write("\t".join([
                "dxyx", 
                geno1, 
                geno2, 
                str(dxy_value), 
                str(dx_values.get(geno1, 'na')),
                str(dxy_value - dx_values[geno1]) if geno1 in dx_values else 'na', 
                str(dx_values.get(geno2, 'na')),
                str(dxy_value - dx_values[geno2]) if geno2 in dx_values else 'na',
                str(len(pop1)), 
                str(len(pop2)),
                str(dxy_comparisons)
                ]) + 
                "\n")
        fafile.close()
    
outFile.close()
print("Completed in: ", time() - time0, "seconds")