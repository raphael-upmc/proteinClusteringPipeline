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
import multiprocessing as mp


def runningHhblits(hhm_filename,hhblits_database,hhr_filename) :
    basename = os.path.basename(hhm_filename).split('.')[0]
    cmd = '/home/meheurap/programs/hhsuite-3.0-beta.3-Linux/bin/hhblits -i '+hhm_filename+' -o '+hhr_filename+' -d '+hhblits_database+'  -v 0 -p 50 -E 0.001 -z 1 -Z 32000 -B 0 -b 0 -n 2 -cpu 1'
    status = os.system(cmd)
    if status == 0 :
        return basename,True
    else:
        return basename,False

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
    cmd = '/home/meheurap/programs/hhsuite-3.0-beta.3-Linux/bin/hhmake -add_cons -M 50 -diff 100 -i '+a3m_filename+' -o '+hhm_filename+' >/dev/null 2>&1'
    status = os.system(cmd)
    if status :
        cpt += 1

    if cpt == 0 :
        return basename,True
    else :
        return basename,False

def creatingHhblitsDb(output_directory) :

    output = open(output_directory+'/'+'a3m.list','w')
    for root, dirs, files in os.walk(output_directory+'/'+'a3m'):
        for filename in files :
            output.write(filename+'\n')
    output.close()


    cmd = 'cd '+output_directory+'/'+'a3m'+' && '+'ffindex_build -as '+'..'+'/'+'db_a3m.ffdata '+'..'+'/'+'db_a3m.ffindex '+'-f ../a3m.list'+' >/dev/null 2>&1'
    status = os.system(cmd)
    if status != 0 :
        return False


    output = open(output_directory+'/'+'hhm.list','w')
    for root, dirs, files in os.walk(output_directory+'/'+'hhm'):
        for filename in files :
            output.write(filename+'\n')
    output.close()

    
    cmd = 'cd '+output_directory+'/'+'hhm'+' && '+'ffindex_build -as '+'..'+'/'+'db_hhm.ffdata '+'..'+'/'+'db_hhm.ffindex '+'-f ../hhm.list'+' >/dev/null 2>&1'
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
    subfamily2defline2seqList = dict()
    defline = None
    file = open(msa_filename,'r')
    for line in file :
        line = line.strip('\t\r\n\0')
        if re.search(r'^>',line):
            if defline != None and defline in orf2family : # store the sequence
                subfamily = orf2family[ defline ]
                if subfamily in subfamily2defline2seqList :
                    subfamily2defline2seqList[ subfamily ][ defline ] =  SeqRecord(seq=Seq(seq,IUPAC.protein),id=defline,description="")
                else:
                    subfamily2defline2seqList[ subfamily ] = defaultdict(list)
                    subfamily2defline2seqList[ subfamily ][ defline ] =  SeqRecord(seq=Seq(seq,IUPAC.protein),id=defline,description="")

            # new sequence to create
            defline = line.replace('>','').split()[0]
            seq = ''
        else:
            seq += line


    # the last sequence of the MSA file ==> sneaky sequence !
    if defline != None and defline in orf2family : # store the sequence
        subfamily = orf2family[ defline ]
        if subfamily in subfamily2defline2seqList :
            subfamily2defline2seqList[ subfamily ][ defline ] =  SeqRecord(seq=Seq(seq,IUPAC.protein),id=defline,description="")
        else:
            subfamily2defline2seqList[ subfamily ] = defaultdict(list)
            subfamily2defline2seqList[ subfamily ][ defline ] =  SeqRecord(seq=Seq(seq,IUPAC.protein),id=defline,description="")

    return subfamily2defline2seqList
        
def checkingMSA(fasta_filename,seqId2seq,nb) :
    msaSet = set()
    for defline in seqId2seq :
        msaSet.add(defline)
    
    
    fastaSet = set()
    cpt = 0
    for record in SeqIO.parse(fasta_filename,'fasta') :
        cpt += 1
        fastaSet.add(record.id)
        if record.id not in seqId2seq :
            print(record.id)
            return False
        else:
            continue
            msa = str(seqId2seq[ record.id ].seq).replace('-','')
            if not re.search(msa,str(record.seq)) :
                print(record.id)

    if msaSet.issubset(fastaSet) :
        return True
    else:
        return False

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run all-vs-all hhblits on the subfamilies')
    parser.add_argument('config_filename', help='the path of the CONFIG_FILENAME created by the subfamilies.py script')
    parser.add_argument('--cpu',type=int,default=1,help='number of CPUs used by mmseqs (default: 1)')
    parser.add_argument('--min-size',type=int,default=2,help='minimal size of the protein families to consider (default: 2)')

    args = parser.parse_args()

    
    t1 = datetime.now()

    if not os.path.exists(args.config_filename) :
        sys.exit(args.config_filename+' does not exist, exit')
    else:
        config_filename = os.path.abspath(args.config_filename)
        
    with open(config_filename) as f:
        data = json.load(f)

    subfamily2nb = data['clusters']
    msa_filename = data['msa_filename']
    cwd = data['directory']

    subfam_directory = cwd+'/'+'subfamiliesFasta'
    subfam_directory = os.path.abspath(subfam_directory)

    tsv_filename = os.path.abspath(cwd+'/'+'orf2subfamily.tsv')
    orf2subfamily = orf2familyFunction(tsv_filename)


    # creating a log file
    logging_filename = cwd+"/logs/"+'hhblits.log'
    logging.basicConfig(filename=logging_filename, level=logging.INFO,filemode="w")
    logging.getLogger(__name__)
    logging.info('creating log file '+logging_filename+"\n")
    logging.info('datetime start: '+str(t1)+'\n')

    logging.info("command line: "+' '.join(sys.argv))
    logging.info('config_filename: '+config_filename)
    logging.info('cpu: '+str(args.cpu))
    logging.info('min-size: '+str(args.min_size)+'\n')   


    

    ##########################
    # creating the directory #
    ##########################   
    hhblits_directory = os.path.abspath(cwd+'/'+'hhblits')
    a3m_directory = os.path.abspath(hhblits_directory+'/'+'a3m')
    hhm_directory = os.path.abspath(hhblits_directory+'/'+'hhm')
    hhr_directory = os.path.abspath(hhblits_directory+'/'+'hhr')

    logging.info('creating directories')
    logging.info('\t'+hhblits_directory)
    logging.info('\t'+a3m_directory)
    logging.info('\t'+hhm_directory)
    logging.info('\t'+hhr_directory+'\n')    
    
    if os.path.exists(hhblits_directory) :
        logging.error(hhblits_directory+' already exists, remove it')
        sys.exit(hhblits_directory+' already exists, remove it')
    else:
        os.makedirs(hhm_directory)
        os.makedirs(a3m_directory)
        os.makedirs(hhr_directory)    


    ###################################################
    # checking msa file and creating individual files #
    ###################################################

    problematicSubfamiliesSet = set()
    
    print('checking '+msa_filename+'...')    
    logging.info('checking '+msa_filename+'...')

    subfamily2defline2seqList = readingMsa(msa_filename,orf2subfamily)
    error = 0
    for subfamily,nb in subfamily2nb.items() :
        nb = int(nb)

        if nb < args.min_size :
            continue

        fasta_filename = subfam_directory+'/'+subfamily+'.fa'
        if not checkingMSA(fasta_filename,subfamily2defline2seqList[ subfamily ],nb) :
            logging.error(subfamily+'\t'+str(nb)+' ==> ERROR')
            problematicSubfamiliesSet.add(subfamily)
            error += 1

        a3m_filename = os.path.abspath(a3m_directory+'/'+subfamily+'.a3m')
        creatingA3m(a3m_filename,subfamily,subfamily2defline2seqList[ subfamily ])

    if error != 0 :
        logging.info('\n'+str(error)+' subfamilies failed to checking\n')


    logging.info('done\n')
    print('done')
    
    #################################
    # creating the hhblits database #
    #################################
    t2 = datetime.now()
    print('\ncreating the hhblits database...')
    logging.info('creating the hhblits database... ('+str(t2)+')')

    # parallelizing
    results = list()
    pool = mp.Pool(processes=args.cpu,maxtasksperchild=1) # start 20 worker processes and 1 maxtasksperchild in order to release memory
    
    for subfamily,nb in subfamily2nb.items() :

        nb = int(nb)
        if nb < args.min_size :
            continue

        a3m_filename = os.path.abspath(a3m_directory+'/'+subfamily+'.a3m')
        results.append( pool.apply_async( creatingFiles4HhblitsDb, args= (a3m_filename,hhm_directory,) ))      
    pool.close() # Prevents any more tasks from being submitted to the pool
    pool.join() # Wait for the worker processes to exit


    error = 0
    for elt in results :
        subfamily,result = elt.get()
        if not  result :
            logging.error('\t'+subfamily+' '+'==>'+' '+'Error' )
            problematicSubfamiliesSet.add(subfamily)
            error += 1

    if error != 0 :
        logging.info('\n'+str(error)+' subfamilies failed to be transformed into hhm\n')


    # removing problematic subfamilies to avoid error during the hhblits db creation
    if len(problematicSubfamiliesSet) != 0 :
        logging.info('removing problematic subfamilies to avoid error during the hhblits db creation. Those files will not be considered for the Hmm-Hmm comparison')
        for subfamily in problematicSubfamiliesSet :
            a3m_filename = os.path.abspath(a3m_directory+'/'+subfamily+'.a3m')
            if os.path.exists(a3m_filename) :                
                os.remove(a3m_filename)
                logging.info('\t'+a3m_filename)

            hhm_filename = os.path.abspath(hhm_directory+'/'+subfamily+'.hhm')
            if os.path.exists(hhm_filename) :
                os.remove(hhm_filename)
                logging.info('\t'+hhm_filename)
        logging.info('removing files done\n')
    else:
        logging.info('No problematic subfamilies detected!')

    if creatingHhblitsDb(hhblits_directory) :
        logging.info('hhblits database creation done\n')
    else:
        logging.error('something went wrong during the hhblits database creation! exit')
        sys.exit('something went wrong during the hhblits database creation! exit')
    print('done\n')
    



    ###################
    # running hhblits #
    ###################
    t3 = datetime.now()
    logging.info('running the hhblits... ('+str(t3)+')')
    print('running the hhblits...')

    # parallelizing
    results = list()
    pool = mp.Pool(processes=args.cpu,maxtasksperchild=1) # start 20 worker processes and 1 maxtasksperchild in order to release memory
    
    for subfamily,nb in subfamily2nb.items() :

        nb = int(nb)
        if nb < args.min_size :
            continue

        if subfamily in problematicSubfamiliesSet :
            continue

        hhr_filename = os.path.abspath(hhr_directory+'/'+subfamily+'.hhr')
        hhm_filename = os.path.abspath(hhm_directory+'/'+subfamily+'.hhm')
        hhblits_database = hhblits_directory+'/'+'db'

        results.append( pool.apply_async( runningHhblits, args= (hhm_filename,hhblits_database,hhr_filename,) ))      
    pool.close() # Prevents any more tasks from being submitted to the pool
    pool.join() # Wait for the worker processes to exit
        
    error = 0
    for elt in results :
        subfamily,result = elt.get()
        if not  result :
            logging.error('\t'+subfamily+' '+'==>'+' '+'Error' )
            error += 1

    if error != 0 :
        logging.info('\n'+str(error)+' hhm failed to run hhblits\n')


    logging.info('done\n')
    print('done')

    
    t4 = datetime.now()
    logging.info('Script completed at '+str(t4))

    sys.exit(0)
