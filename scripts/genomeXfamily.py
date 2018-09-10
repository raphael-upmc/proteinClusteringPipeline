#! /usr/bin/env python


import os,sys,re
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='creating a matrix with genomes as rows and family as columns')
parser.add_argument('orf2family_filename', help='the path of the ORF2FAMILY_FILENAME that is a tab-separated file with col1: orf name col2: family. WARNING: the first line is skipped (header)')
parser.add_argument('orf2genome_filename', help='the path of the ORF2GENOME_FILENAME that is a tab-separated file with col1: orf name col2: genome. WARNING: the first line is skipped (header)')
parser.add_argument('matrix_filename', help='the path of the MATRIX_FILENAME')
parser.add_argument('--min-size',type=int,default=3,help='minimal number of distinct genomes where the protein families is present to be reported in the matrix (default: 3)')
parser.add_argument('--count',action='store_true',default=False,help='create a count matrix (default: create a binary matrix)')
parser.add_argument('--genomeList',help='the path to the GENOME_LIST_FILENAME to consider')
args = parser.parse_args()


# checking arguments
if not os.path.exists(args.orf2family_filename) :
    sys.exit(args.orf2family_filename+' does not exist, exit')
else:
    orf2family_filename = os.path.abspath(args.orf2family_filename)
    
if not os.path.exists(args.orf2genome_filename) :
    sys.exit(args.orf2genome_filename+' does not exist, exit')
else:
    orf2genome_filename = os.path.abspath(args.orf2genome_filename)

genomeList = 0    
if args.genomeList != None and os.path.exists(args.genomeList) :    
    genomeList = 1
    genomeSet = set()
    file = open(args.genomeList,'r')
    for line in file :
        genome = line.rstrip()
        genomeSet.add(genome)
    file.close()

    
family_threshold = args.min_size
        
orf2genome = dict()
file = open(orf2genome_filename,'r')
next(file)
for line in file :
    line = line.rstrip()
    orf,genome = line.split('\t')
    orf2genome[ orf ] = genome
file.close()

genome2family2count = dict()
family2genomes = defaultdict(set)
genome2family = defaultdict(set)
file = open(orf2family_filename,'r')
next(file)
for line in file :
    line = line.rstrip()
    orf,family = line.split('\t')
    genome = orf2genome[ orf ]
    if genomeList == 1 and genome not in genomeSet :        
            continue
        
    genome2family[ genome ].add(family)
    family2genomes[ family ].add(genome)
    if genome not in genome2family2count :
        genome2family2count[genome] = defaultdict(int)
        genome2family2count[genome][family] += 1
    else:
        genome2family2count[genome][family] += 1
file.close()

familySet = set()
for family,genomeSet in family2genomes.items() :
    if len(genomeSet) < family_threshold :
        continue
    else:
        familySet.add(family)

family2profile = defaultdict(list)
output = open(args.matrix_filename,'w')
output.write('\t'+'\t'.join(list(familySet))+'\n')
for genome in genome2family :
    output.write(genome)
    for family in familySet :
        if family  in genome2family[ genome ]  :
            if not args.count :
                output.write('\t'+'1')
            else:
                output.write('\t'+str( genome2family2count[genome][family] ) )
        else:
            output.write('\t'+'0')
    output.write('\n')
output.close()


print('Genome X family matrix created ('+str(len(genome2family))+' genomes) ('+str(len(familySet))+' families)')
