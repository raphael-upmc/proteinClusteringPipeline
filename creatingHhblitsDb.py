#! /usr/bin/env python

import os,sys,re
from collections import defaultdict


def addingNameToA3m(input_filename,output_filename,name) :
    output = open(output_filename,'w')
    output.write('#'+name+'\n')
    
    file=open(input_filename,'r')
    for line in file :
        output.write(line)
    file.close()
    output.close()

cluster_filename = sys.argv[1]

output_directory = cluster_filename.split('/')[0]+'/'+'database'+'/'+'hhblits'
os.mkdir(output_directory+'/'+'aln')
os.mkdir(output_directory+'/'+'a3m')
os.mkdir(output_directory+'/'+'cs219')
os.mkdir(output_directory+'/'+'hhm')

fasta_directory = cluster_filename.split('/')[0]+'/'+'output'+'/'+'subfamily'+'/'+'fasta'

for (path, dirs, files) in os.walk(fasta_directory):
    for filename in files :
        fasta_filename = path+'/'+filename
        mafft_filename = output_directory+'/'+'aln'+'/'+filename.replace('.fa','.mafft')
        cmd = 'mafft --auto '+fasta_filename+' > '+mafft_filename
        print(cmd)
        status = os.system(cmd)
        if status != 0 :
            continue
            print(cmd)


        a3m_filename_noss = output_directory+'/'+'aln'+'/'+filename.replace('.fa','.a3m')
        cmd = 'reformat.pl fas a3m '+mafft_filename+' '+a3m_filename_noss
        print(cmd)
        status = os.system(cmd)
        if status != 0 :
            print(cmd)
        a3m_filename_noss_renamed = output_directory+'/'+'aln'+'/'+filename.replace('.fa','_renamed.a3m')
        name = filename.replace('.fa','')
        addingNameToA3m(a3m_filename_noss,a3m_filename_noss_renamed,name)
            
        a3m_filename = output_directory+'/'+'a3m'+'/'+filename.replace('.fa','.a3m')
        cmd = 'addss.pl '+a3m_filename_noss_renamed+' '+a3m_filename+' -a3m'
        print(cmd)
        status = os.system(cmd)
        if status != 0 :
            print(cmd)
       
        hhm_filename = output_directory+'/'+'hhm'+'/'+filename.replace('.fa','.hhm')
        cmd = 'hhmake -i '+a3m_filename+' -o '+hhm_filename+' -name '+filename.replace('.fa','')
        print(cmd)
        status = os.system(cmd)
        if status != 0 :
            print(cmd)

            

cmd = 'cd '+output_directory+'/'+'a3m'+' && '+'ffindex_build -as '+'..'+'/'+'db_a3m.ffdata '+'..'+'/'+'db_a3m.ffindex '+'*.a3m'
print(cmd)
status = os.system(cmd)
print(status)
cmd = 'cd '+output_directory+'/'+'hhm'+' && '+'ffindex_build -as '+'..'+'/'+'db_hhm.ffdata '+'..'+'/'+'db_hhm.ffindex '+'*.hhm'
print(cmd)
status = os.system(cmd)
print(status)


cmd = 'cstranslate -A $HHLIB/data/cs219.lib -D $HHLIB/data/context_data.lib -I a3m -x 0.3 -c 4'+' -f -i '+output_directory+'/'+'db_a3m'+' -o '+output_directory+'/'+'db_cs219'+' -b'
print(cmd)
status = os.system(cmd)
print(status)
