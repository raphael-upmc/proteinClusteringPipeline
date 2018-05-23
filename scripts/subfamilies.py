#! /usr/bin/env python


import os,sys,re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import logging
from datetime import date, time, datetime


if __name__ == "__main__":
    t1 = datetime.now()
    cwd = os.path.abspath(os.getcwd())
    
    parser = argparse.ArgumentParser(description='Run protein sequences clustering')
    parser.add_argument('fasta_filename', help='the path of the FASTA_FILENAME that contains the proteins sequences to cluster')
    parser.add_argument('--output-directory',help='the output directory where the results will be store (default: ./FASTA_FILENAME_proteinClutering)')
    parser.add_argument('--cpu',type=int,default=1,help='number of CPUs used by mmseqs (default: 1)')

    args = parser.parse_args()


    # checking arguments
    if not os.path.exists(args.fasta_filename) :
        sys.exit(args.fasta_filename+' does not exist, exit')
    else:
        fasta_filename = os.path.abspath(args.fasta_filename)


    if args.output_directory == None :
        directory = cwd+'/'+os.path.basename(fasta_filename)+'_proteinClustering'
    else:    
        directory = args.output_directory


    
    # checking and creating the output directory
    if not os.path.exists(directory) :
        try :
            directory = os.path.realpath(directory)
            os.mkdir(directory)                
        except:
            print(directory+' cannot be created, creating the output in: '+cwd+'/'+os.path.basename(fasta_filename)+'_proteinClustering')        
            directory = cwd+'/'+os.path.basename(fasta_filename)+'_proteinClustering'
            if not os.path.exists(directory) :
                os.mkdir(directory)
            else:            
                sys.exit(os.path.realpath(directory)+' already exists, remove it first! exit')
    else:
        sys.exit(os.path.realpath(directory)+' already exists, remove it first! exit')


    # creating a log file
    os.mkdir(directory+"/"+"logs")
    logging_filename = directory+"/logs/"+'subfamilies.log'
    logging.basicConfig(filename=logging_filename, level=logging.INFO,filemode="w")
    logging.getLogger(__name__)
    logging.info('creating '+directory+"/"+"logs and log file "+logging_filename+"\n")
    logging.info('datetime start: '+str(t1)+'\n')


    logging.info("command line: "+' '.join(sys.argv))
    logging.info('fasta_filename: '+fasta_filename)    
    logging.info('output directory: '+directory)
    logging.info('cpu: '+str(args.cpu)+'\n')




    mmseqs_directory = directory+'/mmseqs'
    os.mkdir(directory+'/mmseqs')
    logging.info('creating '+mmseqs_directory)


    db_filename = mmseqs_directory+'/'+os.path.basename(fasta_filename)+'.mmseqsDB'
    log_filename = directory+"/logs/"+'mmseqs_createdb.log'
    mmseqs_db_cmd = '/home/meheurap/programs/mmseqs2/bin/mmseqs createdb '+fasta_filename+' '+db_filename+' >'+log_filename
    logging.info('running '+mmseqs_db_cmd)
    status = os.system(mmseqs_db_cmd)


    tmp_directory = mmseqs_directory+'/tmp'
    logging.info('creating '+tmp_directory)
    if os.path.exists(db_filename) :
        os.mkdir(tmp_directory)
    else :
        logging.error(db_filename+' is absent! exit')
        sys.exit(db_filename+' is absent! exit')


    cluster_filename = mmseqs_directory+'/'+os.path.basename(fasta_filename)+'.mmseqsDB_clu'
    log_filename = directory+'/logs/'+'mmseqs_cluster.log'
    mmseqs_cluster_cmd = '/home/meheurap/programs/mmseqs2/bin/mmseqs cluster '+db_filename+' '+cluster_filename+' '+tmp_directory+' '+'--threads '+str(args.cpu)+' -s 7.5 -c 0.5 --cov-mode 0 --max-seqs 5000 -e 0.001 --cluster-mode 0'+' >'+log_filename
    logging.info('running '+mmseqs_cluster_cmd)
    status = os.system(mmseqs_cluster_cmd)


    tsv_filename = mmseqs_directory+'/'+os.path.basename(fasta_filename)+'.mmseqsDB_clu.tsv'
    log_filename = directory+'/logs/'+'mmseqs_createtsv.log'
    mmseqs_createtsv_cmd = '/home/meheurap/programs/mmseqs2/bin/mmseqs createtsv '+db_filename+' '+db_filename+' '+cluster_filename+' '+tsv_filename+' >'+log_filename
    logging.info('running '+mmseqs_createtsv_cmd)
    status = os.system(mmseqs_createtsv_cmd)



    #################################
    ### writting the final output ###
    #################################

    logging.info('reading the output and renaming the subfamilies...')
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
    clusteredSeq = 0
    cluster2subfam = dict()
    subfam = 0
    for cluster,sequenceList in cluster2sequences.items() :
        if len(sequenceList) == 1 : # removing orphans
            continue
        else :
            subfam += 1
            cluster2subfam[cluster] = subfam
            clusteredSeq += len(sequenceList)
    l = len(str(subfam))

    logging.info('On the '+str(len(sequence2cluster))+' protein sequences, '+str(clusteredSeq)+' were clustered into '+str(subfam)+' subfams')




    cluster2seqList = defaultdict(list)
    for record in SeqIO.parse(fasta_filename,'fasta') :
        if record.id in sequence2cluster :
            cluster2seqList[ sequence2cluster[ record.id ] ].append( SeqRecord(seq=record.seq,id=record.id,description="") )
        else :
            continue

    logging.info('creating one fasta file per subfamily in '+directory+'/'+'subfamiliesFasta')
    os.mkdir(directory+'/'+'subfamiliesFasta')

    subfam2nb = defaultdict(int)
    seq2subfam_filename = directory+'/orf2subfamily.tsv'
    logging.info('creating tsv file '+seq2subfam_filename)
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
            subfam2nb[ subfamName ] += 1
            
        fasta_output_filename = directory+'/'+'subfamiliesFasta'+'/'+str(subfamName)+'.fa'
        SeqIO.write(cluster2seqList[cluster],fasta_output_filename,'fasta')
        
    output.close()



    config_filename = directory+'/config.json'
    logging.info('creating the config file for HMM-HMM comparison '+config_filename+'\n')
    output = open(config_filename,'w')
    liste = list()
    output.write('{'+'\n')
    output.write('\t\"directory\":\"'+directory+'\",\n')
    output.write('\t\"clusters\":{\n')
    for (path, dirs, files) in os.walk(directory+'/'+'subfamiliesFasta') :
        for filename in files :
#            liste.append('\t\t\"'+filename.replace('.fa','\"')+':\"'+path+'/'+filename+'\"')
            liste.append('\t\t\"'+filename.replace('.fa','\"')+':\"'+str(subfam2nb[ filename.replace('.fa','') ])+'\"')
    output.write(',\n'.join(liste)+'\n')
    output.write('\t'+'}'+'\n')
    output.write('}'+'\n')
    output.close()

    t2 = datetime.now()
    t3 = t2 - t1

    logging.info("done in "+str(t3.seconds)+" sec. ("+str(t3.seconds/3600)+" hours)")
    logging.info('datetime end: '+str(t2))
    sys.exit()