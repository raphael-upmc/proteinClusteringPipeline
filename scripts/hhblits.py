#! /usr/bin/env python

import json
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys,os,re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import logging
from datetime import date, time, datetime

def creatingFiles4HhblitsDb(a3m_filename,hhm_directory) :
    basename = os.path.basename(a3m_filename).split('.')[0]

    cpt = 0
    
    # reformat.pl
    cmd = '/home/meheurap/programs/hhsuite-3.0-beta.3-Linux/scripts/reformat.pl a3m a3m '+a3m_filename+' '+a3m_filename+' -M 50 -r >/dev/null 2>&1'
    status = os.system(cmd)
    if status :
        cpt += 1

    # addss.pl
    cmd = '/home/meheurap/programs/hhsuite-3.0-beta.3-Linux/scripts/addss.pl '+a3m_filename+' '+a3m_filename+' -a3m >/dev/null 2>&1'
    status = os.system(cmd)
    if status :
        cpt += 1


    # hhmake
    hhm_filename = hhm_directory+'/'+basename+'.hhm'
    cmd = '/home/meheurap/programs/hhsuite-3.0-beta.3-Linux/bin/hhmake -add_cons -M 50 -i '+a3m_filename+' -o '+hhm_filename+' >/dev/null 2>&1'
    status = os.system(cmd)
    if status :
        cpt += 1

    if cpt == 0 :
        return True
    else :
        return False

def creatingHhblitsDb(output_directory) :
    cmd = 'cd '+output_directory+'/'+'a3m'+' && '+'ffindex_build -as '+'..'+'/'+'db_a3m.ffdata '+'..'+'/'+'db_a3m.ffindex '+'*.a3m'+' >/dev/null 2>&1'
    status = os.system(cmd)
    if status != 0 :
        return False

    cmd = 'cd '+output_directory+'/'+'hhm'+' && '+'ffindex_build -as '+'..'+'/'+'db_hhm.ffdata '+'..'+'/'+'db_hhm.ffindex '+'*.hhm'+' >/dev/null 2>&1'
    status = os.system(cmd)
    if status != 0 :
        return False

    cmd = 'cstranslate -A $HHLIB/data/cs219.lib -D $HHLIB/data/context_data.lib -I a3m -x 0.3 -c 4'+' -f -i '+output_directory+'/'+'db_a3m'+' -o '+output_directory+'/'+'db_cs219'+' -b'+' >/dev/null 2>&1'
    status = os.system(cmd)
    if status != 0 :
        return False

    return True



def orf2familyFunction(tsv_filename) :
    orf2family = dict()
    file = open(tsv_filename,'r')
    next(file)
    for line in file :
        line = line.rstrip()
        orf,subfamily = line.split('\t')
        orf2family[ orf ] = subfamily
    file.close()
    return orf2family


def creatingA3m(a3m_filename,subfamily,defline2seq) :
    output = open(a3m_filename,'w')
    output.write('#'+subfamily+'\n')
    for defline,record in defline2seq.items() :
        SeqIO.write(record,output,'fasta')

def readingMsa(msa_filename,orf2family) :
    subfamily2seqList = defaultdict(list)
    tag = 0
    file = open(msa_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.search('>',line):
            if tag == 1 :
                if subfamily in subfamily2seqList :
                    subfamily2seqList[ subfamily ][ defline.replace('>','') ] =  SeqRecord(seq=Seq(seq,IUPAC.protein),id=defline,description="")
                else:
                    subfamily2seqList[ subfamily ] = defaultdict(list)
                    subfamily2seqList[ subfamily ][ defline.replace('>','') ] =  SeqRecord(seq=Seq(seq,IUPAC.protein),id=defline,description="")
        
            defline = line.split('>')[1].split()[0]
            if defline in orf2subfamily :
                tag = 1
                subfamily = orf2subfamily[ defline ]
                seq = ''
            else:
                tag = 0
        else:
            if tag == 1 :
                seq += line
            else:
                continue
    return subfamily2seqList
        
def checkMSA(fasta_filename,seqId2seq,nb) :
    cpt = 0
    for record in SeqIO.parse(fasta_filename,'fasta') :
        cpt += 1
        if record.id not in seqId2seq :
            print(record.id)
            return False
        else:
            continue
            msa = str(seqId2seq[ record.id ].seq).replace('-','')
            if not re.search(msa,str(record.seq)) :
                print(record.id)

    if nb == cpt and nb == len(seqId2seq) :
        return True
    else:
        return False

if __name__ == "__main__":
    t1 = datetime.now()
    
    config_filename = sys.argv[1]
    if not os.path.exists(config_filename) :
        sys.exit(config_filename+' does not exist, exit')
    
    with open(config_filename) as f:
        data = json.load(f)

    subfamily2nb = data['clusters']
    msa_filename = data['msa_filename']
    cwd = data['directory']

    subfam_directory = cwd+'/'+'subfamiliesFasta'
    subfam_directory = os.path.abspath(subfam_directory)

    tsv_filename = os.path.abspath(cwd+'/'+'orf2subfamily.tsv')
    orf2subfamily = orf2familyFunction(tsv_filename)

    subfamily2seqList = readingMsa(msa_filename,orf2subfamily)


    # creating a log file
    logging_filename = cwd+"/logs/"+'hhblits.log'
    logging.basicConfig(filename=logging_filename, level=logging.INFO,filemode="w")
    logging.getLogger(__name__)
    logging.info('creating log file '+logging_filename+"\n")
    logging.info('datetime start: '+str(t1)+'\n')


    logging.info("command line: "+' '.join(sys.argv))
    logging.info('config_filename: '+config_filename)    

    

    ##########################
    # creating the directory #
    ##########################   
    hhblits_directory = os.path.abspath(cwd+'/'+'hhblits')
    a3m_directory = os.path.abspath(hhblits_directory+'/'+'a3m')
    hhm_directory = os.path.abspath(hhblits_directory+'/'+'hhm')

    logging.info('creating directories')
    logging.info('\t'+hhblits_directory)
    logging.info('\t'+a3m_directory)
    logging.info('\t'+hhm_directory)
    
    if os.path.exists(hhblits_directory) :
        logging.error(hhblits_directory+' already exists, remove it')
        sys.exit(hhblits_directory+' already exists, remove it')
    else:
        os.makedirs(hhm_directory)
        os.makedirs(a3m_directory)    


    ###################################################
    # checking msa file and creating individual files #
    ###################################################
    print('checking '+msa_filename+'...')
    logging.info('checking '+msa_filename+'...')
    for root, dirs, files in os.walk(subfam_directory):
        for filename in files :
            subfamily = filename.split('.')[0]

            fasta_filename = root+'/'+filename
            nb = int(subfamily2nb[ subfamily ])
            if not checkMSA(fasta_filename,subfamily2seqList[ subfamily ],nb) :
                print(subfamily+' ==> ERROR')
            else:
                print(subfamily+' ==> okay')

            a3m_filename = os.path.abspath(a3m_directory+'/'+subfamily+'.a3m')
            creatingA3m(a3m_filename,subfamily,subfamily2seqList[ subfamily ])
    logging.info('done')

    
    #################################
    # creating the hhblits database #
    #################################
    print('creating the hhblits database...')
    logging.info('creating the hhblits database...')
    # to parallelize
    
    for subfamily,nb in subfamily2nb.items() :
        a3m_filename = os.path.abspath(a3m_directory+'/'+subfamily+'.a3m')
        cpt = creatingFiles4HhblitsDb(a3m_filename,hhm_directory)
        if cpt :
            print('\t'+subfamily+'\t'+str(nb)+'\t'+'==>'+'\t'+'Okay' )
            logging.info('\t'+subfamily+'\t'+str(nb)+'\t'+'==>'+'\t'+'Okay' )
        else:
            print('\t'+subfamily+'\t'+str(nb)+'\t'+'==>'+'\t'+'Error' )
            logging.error('\t'+subfamily+'\t'+str(nb)+'\t'+'==>'+'\t'+'Okay' )

            
    if creatingHhblitsDb(hhblits_directory) :
        logging.info('hhblits database created!')
        print('hhblits database created!')
    else:
        print('something went wrong during the hhblits database creation!')
        logging.error('something went wrong during the hhblits database creation!')

        
    ###################
    # running hhblits #
    ###################
    # to parallelize
    
