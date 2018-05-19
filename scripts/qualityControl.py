#! /usr/bin/python

import os,sys,re
from collections import defaultdict
import argparse


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Run protein clustering quality control')
    parser.add_argument('annotation_filename', help='the path of the ANNOTATION_FILENAME')
    parser.add_argument('orf2family_filename', help='the path of the ORF2FAMILY_FILENAME')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILENAME')


    args = parser.parse_args()


    # checking arguments
    if not os.path.exists(args.annotation_filename) :
        sys.exit(args.fasta_filename+' does not exist, exit')
    else:
        annotation_filename = os.path.abspath(args.annotation_filename)

    if not os.path.exists(args.orf2family_filename) :
        sys.exit(args.orf2family_filename+' does not exist, exit')
    else:
        orf2family_filename = os.path.abspath(args.orf2family_filename)
        

    # reading annotation filename
    orf2annotation = dict()
    annotation2count = defaultdict(int)
    file = open(annotation_filename,'r')
    for line in file :
        line = line.rstrip()
        orf,annotation = line.split('\t')
        orf2annotation[orf] = annotation
        annotation2count[annotation] += 1
    file.close()

    # reading orf2family filename
    annotation2family = defaultdict(list)
    family2count = defaultdict(int)
    file = open(orf2family_filename,'r')
    for line in file :
        line = line.rstrip()
        orf,family = line.split('\t')
        family2count[family] += 1
        if orf in orf2annotation :
            annotation = orf2annotation[ orf ]
            annotation2family[ annotation ].append(family)
        else:
            continue
    file.close()
    
    # writting the results
    output = open(args.output_filename,'w')
    for annotation,nb in annotation2count.items() :
        output.write(annotation+'\t'+str(nb)+'\n')
        family2annot = defaultdict(int)
        cpt = 0
        for family in annotation2family[ annotation ] :
            family2annot[family] += 1
            cpt += 1
        
        for family,annot in family2annot.items() :
            output.write('\t\t'+str(annot)+'\t'+family+' ('+str(family2count[family])+')'+'\n')

        output.write('\t\t'+'#orfans (# unclustered proteins)'+'\t'+str(annotation2count[annotation] - cpt)+'\n')
    output.close()
