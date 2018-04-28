#! /usr/bin/env python

import sys, os, re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import shutil

cluster_filename = sys.argv[1]
fasta_filename = sys.argv[2]

seq2subfam_filename = sys.argv[3]
config_filename = sys.argv[4]

output_directory = seq2subfam_filename.rsplit('/',1)[0]

cluster2sequences = defaultdict(set)
sequence2cluster = dict()
file = open(cluster_filename,'r')
for line in file :
    line = line.rstrip()
    represent,sequence = line.split('\t')
    sequence2cluster[ sequence ] = represent
    cluster2sequences[ represent ].add(sequence)
file.close()


# renaming cluster name by subfamXXX
cluster2subfam = dict()
subfam = 0
for cluster,sequenceList in cluster2sequences.items() :
    if len(sequenceList) == 1 : # removing orphans
        continue
    else :
        subfam += 1
        cluster2subfam[cluster] = subfam

l = len(str(subfam))

print('there are '+str(subfam)+' subfams')
#print(l)


cluster2seqList = defaultdict(list)
for record in SeqIO.parse(fasta_filename,'fasta') :
    if record.id in sequence2cluster :
        cluster2seqList[ sequence2cluster[ record.id ] ].append( SeqRecord(seq=record.seq,id=record.id,description="") )
    else :
        continue


if os.path.exists(output_directory+'/'+'fasta') :
    shutil.rmtree(output_directory+'/'+'fasta')

os.mkdir(output_directory+'/'+'fasta')

output = open(seq2subfam_filename,'w')        
output.write('orf\tsubfamily\n')
for cluster,subfam in cluster2subfam.items() :
    subfamName = 'subfam'
    nb = l - len(str(subfam))
    for i in range(nb) :
        subfamName += '0'
    subfamName += str(subfam)
    
    for sequence in cluster2sequences[ cluster ] :
        output.write(sequence+'\t'+subfamName+'\n')

    fasta_output_filename = output_directory+'/'+'fasta'+'/'+str(subfamName)+'.fa'
    SeqIO.write(cluster2seqList[cluster],fasta_output_filename,'fasta')
        
output.close()

output = open(config_filename,'w')

liste = list()
output.write('{'+'\n')
output.write('\t\"clusters\":{\n')
for (path, dirs, files) in os.walk(output_directory+'/'+'fasta'):
    for filename in files :
        liste.append('\t\t\"'+filename.replace('.fa','\"')+':\"'+path+'/'+filename+'\"')
output.write(',\n'.join(liste)+'\n')
output.write('\t'+'}'+'\n')
output.write('}'+'\n')
output.close()
