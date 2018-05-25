#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import argparse

def checkingFasta(fasta_filename) :
    seqIdSet = set()
    badSequenceList = list()
    duplicatedSeqIdSet = set()
    
    noFunkyCharacter = r'^[\w]+\*{0,1}$'
    for record in SeqIO.parse(fasta_filename,'fasta') :

        # duplicated seqId
        if record.id in seqIdSet :
            duplicatedSeqIdSet.add(record.id)
            seqIdSet.add(record.id)
        else:
            seqIdSet.add(record.id)

        # funky characters in the aa sequence
        if re.search(noFunkyCharacter,str(record.seq)) :
            continue
        else:
            badSequenceList.append(record.id)
    return badSequenceList,duplicatedSeqIdSet


def writingCleanedFasta(fasta_filename,badSequenceList,duplicatedSeqIdSet,output_filename):

    output = open(output_filename,'w')

    seqSet = set()
    for record in SeqIO.parse(fasta_filename,'fasta') :
        if record.id in set(badSequenceList) :
            continue

        if record.id in duplicatedSeqIdSet :
            if record.id in seqSet :
                continue
            else:
                seqSet.add(seqSet)
                SeqIO.write(record,output,'fasta') # write only one of the duplicated sequence
                continue
        
        SeqIO.write(record,output,'fasta')
        
        
if __name__ == "__main__":

    cwd = os.path.abspath(os.getcwd())
    
    parser = argparse.ArgumentParser(description='a script that looks for funky characters in the aa sequences and for duplicated seqId')
    parser.add_argument('fasta_filename', help='the path of the FASTA_FILENAME that contains the proteins sequences to cluster')
    parser.add_argument('--cleaned-fasta-filename',help='the new fasta_filename cleaned and ready to use in the protein clustering pipeline')

    args = parser.parse_args()


    # checking arguments
    if not os.path.exists(args.fasta_filename) :
        sys.exit(args.fasta_filename+' does not exist, exit')
    else:
        fasta_filename = os.path.abspath(args.fasta_filename)


    # looking for funky characters and duplicated seqId in the fasta file
    print('reading the fasta file, looking for errors')
    badSequenceList,duplicatedSeqIdSet = checkingFasta(fasta_filename)
    print('done')
    if len(badSequenceList) > 0 or len(duplicatedSeqIdSet) > 0 :
        print('Sequences with funky characters: '+str(len(badSequenceList)))
        print('Sequences with duplicated seqId: '+str(len(duplicatedSeqIdSet))+'\n')
        if args.cleaned_fasta_filename != None :
            writingCleanedFasta(fasta_filename,badSequenceList,duplicatedSeqIdSet,os.path.abspath(args.cleaned_fasta_filename))
            print('Creation of a cleaned fasta file: '+os.path.abspath(args.cleaned_fasta_filename))
            sys.exit()
        else:
            sys.exit()
    else:
        print('Sequences with funky characters: '+str(len(badSequenceList)))
        print('Sequences with duplicated seqId: '+str(len(duplicatedSeqIdSet))+'\n')        
        print('Fasta file looks okay: '+fasta_filename)
        if args.cleaned_fasta_filename != None :
            print('No need to create a new fasta file: '+os.path.abspath(args.cleaned_fasta_filename))
            sys.exit()
        else:
            sys.exit()

