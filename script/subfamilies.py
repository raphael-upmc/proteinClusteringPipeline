#! /usr/bin/env python


import os,sys,re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


fasta_filename = sys.argv[1]
directory = sys.argv[2]

cwd = os.path.abspath(os.getcwd())

if not os.path.exists(fasta_filename) :
    sys.exit(fasta_filename+' does not exist, exit')
fasta_filename = os.path.abspath(fasta_filename)
print(fasta_filename)


# checking directory #
if not os.path.exists(directory) :
    try :
        directory = os.path.realpath(directory)
        os.mkdir(directory)                
    except:
        print(directory+' can not be created, creating the outputs in: '+cwd+'/'+os.path.basename(fasta_filename)+'_subfamilies')        
        directory = cwd+'/'+os.path.basename(fasta_filename)+'_subfamilies'
        if not os.path.exists(directory) :
            os.mkdir(directory)
        else:            
            sys.exit(os.path.realpath(directory)+' already exists, remove it first! exit')
else:
    sys.exit(os.path.realpath(directory)+' already exists, remove it first! exit')

print(directory)    

mmseqs_directory = directory+'/mmseqs'
os.mkdir(directory+'/mmseqs')



db_filename = mmseqs_directory+'/'+os.path.basename(fasta_filename)+'.mmseqsDB'

mmseqs_db_cmd = 'mmseqs createdb '+fasta_filename+' '+db_filename
print(mmseqs_db_cmd)
status = os.system(mmseqs_db_cmd)

tmp_directory = mmseqs_directory+'/tmp'
if os.path.exists(db_filename) :
    os.mkdir(tmp_directory)
else :
    sys.exit(db_filename+' is absent! exit')


cluster_filename = mmseqs_directory+'/'+os.path.basename(fasta_filename)+'.mmseqsDB_clu'
mmseqs_cluster_cmd = 'mmseqs cluster '+db_filename+' '+cluster_filename+' '+tmp_directory+' '+'--threads 1 -s 6 -c 0.7 --cov-mode 0 --max-seqs 5000 -e 0.001 --cluster-mode 0'
print(mmseqs_cluster_cmd)
os.system(mmseqs_cluster_cmd)

tsv_filename = mmseqs_directory+'/'+os.path.basename(fasta_filename)+'.mmseqsDB_clu.tsv'
mmseqs_createtsv_cmd = 'mmseqs createtsv '+db_filename+' '+db_filename+' '+cluster_filename+' '+tsv_filename
print(mmseqs_createtsv_cmd)
os.system(mmseqs_createtsv_cmd)



#################################
### writting the final output ###
#################################


cluster2sequences = defaultdict(set)
sequence2cluster = dict()
file = open(tsv_filename,'r')
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

print('the '+str(len(sequence2cluster))+' orfs were clustered in '+str(subfam)+' subfams')




cluster2seqList = defaultdict(list)
for record in SeqIO.parse(fasta_filename,'fasta') :
    if record.id in sequence2cluster :
        cluster2seqList[ sequence2cluster[ record.id ] ].append( SeqRecord(seq=record.seq,id=record.id,description="") )
    else :
        continue

os.mkdir(directory+'/'+'subfamiliesFasta')


seq2subfam_filename = directory+'/orf2subfamily.tsv'
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

    fasta_output_filename = directory+'/'+'subfamiliesFasta'+'/'+str(subfamName)+'.fa'
    SeqIO.write(cluster2seqList[cluster],fasta_output_filename,'fasta')

        
output.close()



