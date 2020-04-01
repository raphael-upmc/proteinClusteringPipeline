#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import argparse
import logging
from datetime import date, time, datetime
import shutil
import json
from Bio import SeqIO

def orf2familyFunction(cwd , subfamily2family) :
    subfamilySet = set()
    orf2subfamily = dict()
    orf2subfamily_filename = cwd+'/'+'orf2subfamily.tsv'
    file = open(orf2subfamily_filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        orf,subfamily = line.split('\t')
        orf2subfamily[ orf ] = subfamily
        subfamilySet.add(subfamily)
    file.close()

    orf2subfamily2family_filename = cwd+'/'+'orf2subfamily2family.tsv'
    output2 = open(orf2subfamily2family_filename,'w')
    output2.write('orf'+'\t'+'subfamily'+'\t'+'family'+'\n')
    
    orf2family_filename = cwd+'/'+'orf2family.tsv'
    output = open(orf2family_filename,'w')
    output.write('orf'+'\t'+'family'+'\n')
    familySet = set()
    for orf,subfamily in orf2subfamily.items() :
        family = subfamily2family[ subfamily ]
        output.write(orf+'\t'+family+'\n')
        output2.write(orf+'\t'+subfamily+'\t'+family+'\n')
        familySet.add(family)
    output.close()
    output2.close()

    cptOrf = len(orf2subfamily)
    cptSubfam = len(subfamilySet)
    cptFam = len(familySet)
    return cptOrf,cptSubfam,cptFam

def creatingMclNetworkFile(hhrHitsList,network_filename) :
    edge2weight = dict()
    for hit in sorted(hhrHitsList,key= lambda x:x[0]) :
        query = hit[0]
        subject = hit[1]
        probs = float(hit[2])/100.0
        qcover = float(hit[3])
        scover = float(hit[4])
        weight = probs * max([qcover,scover])
        edge = '\t'.join( sorted( [query,subject] ) )
        if edge in edge2weight :
            if weight > edge2weight[edge] :
                edge2weight[ edge ] = weight
        else:
            edge2weight[ edge ] = weight

    output = open(network_filename,'w')
    for edge,weight in edge2weight.items() :
        output.write(edge+'\t'+str(weight)+'\n')
    output.close()


    
def parsingMclFile(mcl_filename,subfam2nb) :
    subfam2fam = dict()
    cpt = 0 # number of family
    file = open(mcl_filename,'r')
    for line in file :
        cpt += 1 
        line = line.rstrip() # one line ==> one cluster
        liste = line.split()
        for subfam in liste :
            subfam2fam[ subfam ] = cpt
    file.close()

    for subfam,nb in subfam2nb.items() :
        if subfam not in subfam2fam :
            cpt += 1
            subfam2fam[ subfam ] = cpt
        else:
            continue

    l1 = len(str(cpt))
    for subfam,fam in subfam2fam.items() :
        l2 = len(str(fam))
        family = 'fam'+''.join( ['0' for s in range(l1 - l2)] ) + str(fam)
        subfam2fam[ subfam ] = family

    print(len(subfam2fam))
    print(cpt)
    return subfam2fam

def readingHhrFile(hhr_filename,coverage_threshold,probs_threshold) :
#    print(hhr_filename)
    result = list()
    tag = 0
    file = open(hhr_filename,'r')
    for line in file :
        line = line.strip()
        if re.match(r'Query',line) :
            query = line.split()[1]
            continue
        
        if re.match(r'Match_columns',line) :
            qlen = float( line.split()[1] )
            continue

        if re.match(r'No Hit',line) :
            tag = 1
            continue
        
        if tag == 1 and line != '':
            liste = line.split()
            if len(liste) == 11 :
                target = liste[1]
                probs = float(liste[2])
                qcoord = liste[8]
                qstart = float(qcoord.split('-')[0])
                qend = float(qcoord.split('-')[1])            
                qcover = ( qend - qstart + 1.0 ) / qlen            
                scoord = liste[9]
                sstart = float(liste[9].split('-')[0])
                send = float(liste[9].split('-')[1])            
                slen = float(liste[10].replace('(','').replace(')',''))
                scover = ( send - sstart + 1 ) / slen
            elif len(liste) == 10: # special case when liste[9] = 120-198(200)
                target = liste[1]
                probs = float(liste[2])
                qcoord = liste[8]
                qstart = float(qcoord.split('-')[0])
                qend = float(qcoord.split('-')[1])            
                qcover = ( qend - qstart + 1.0 ) / qlen            
                scoord,slen = liste[9].split('(')
                sstart = float(scoord.split('-')[0])
                send = float(scoord.split('-')[1])            
                slen = float(slen.replace('(','').replace(')',''))
                scover = ( send - sstart + 1 ) / slen
            else:
                sys.exit('error')
                continue
            if probs > probs_threshold and ( scover > coverage_threshold or qcover > coverage_threshold ) :
                if query == target : # removing self loop
                    continue

                if query < target :
                    result.append([query,target,probs,qcover,scover])
                else:
                    result.append([target,query,probs,scover,qcover])

    file.close()
    return result



def checkingDirectory(directory,cluster2nb):
    subfamSet = set()
    for root, dirs, files in os.walk(directory):
        for filename in files :
            subfam = filename.split('.')[0]
            subfamSet.add(subfam)
            if subfam not in cluster2nb :
                print(subfam)
                return False
            else:
                continue
    
    if len(subfamSet) != len(cluster2nb) :
        return False
    else:
        return True



if __name__ == "__main__":
    t1 = datetime.now()
    
    parser = argparse.ArgumentParser(description='Run protein sequences clustering')
    parser.add_argument('config_filename', help='the path of the CONFIG_FILE')
    parser.add_argument('--coverage',type=float,default=0.75,help='coverage threshold (default: 0.75)')
    parser.add_argument('-I',type=float,default=2,help='mcl inflation parameter (default: 2.0)')
    parser.add_argument('--probs',type=float,default=0.95,help='probability threshold (default: 0.95)')
    parser.add_argument('--min-size',type=int,default=5,help='minimal size of the protein families to keep (default: 5)')
    parser.add_argument('--force',action='store_true',default=False,help='force MCL clustering (default: False)')
    parser.add_argument('--fasta',action='store_true',default=False,help='creating a folder with fasta file for each family (default: False)')
    parser.add_argument('--cpu',type=int,default=1,help='number of CPUs used by mcl (default: 1)')
    
    args = parser.parse_args()


    # checking arguments
    if not os.path.exists(args.config_filename) :
        sys.exit(args.config_filename+' does not exist, exit')
    else:
        config_filename = os.path.abspath(args.config_filename)


    with open(config_filename) as f:
        data = json.load(f)

    cwd = data['directory']
    subfam2nb = data['clusters']

    if args.fasta :
        if os.path.exists(cwd+'/'+'familiesFasta') :
            shutil.rmtree(cwd+'/'+'familiesFasta', ignore_errors=True)

            

    #######################
    # creating a log file #
    #######################
    
    logging_filename = cwd+"/logs/"+'runningMclClustering.log'
    logging.basicConfig(filename=logging_filename, level=logging.INFO,filemode="w")
    logging.getLogger(__name__)
    logging.info('creating log file '+logging_filename+"\n")
    logging.info('datetime start: '+str(t1)+'\n')

    logging.info("command line: "+' '.join(sys.argv))
    logging.info('config_filename: '+config_filename)
    logging.info('coverage threshold: '+str(args.coverage))
    logging.info('probs threshold: '+str(args.probs))   
    logging.info('min-size: '+str(args.min_size)+'\n')
    logging.info('force: '+str(args.force))
    logging.info('fasta: '+str(args.fasta)+'\n')   
    logging.info('cpu: '+str(args.cpu))
    logging.info('inflation parameter: '+str(args.I))
    
    ###################################
    # checking if hhblits runned fine #
    ###################################

    logging.info('checking the hhm and hhr files...')
    fasta_dir = cwd+'/'+'subfamiliesFasta'
    hhm_dir = cwd+'/'+'hhblits'+'/'+'hhm'
    hhr_dir = cwd+'/'+'hhblits'+'/'+'hhr'

    if checkingDirectory(hhm_dir,subfam2nb) :
        logging.info(hhm_dir+' looks okay')
    else:
        if not args.force :
            logging.error(hhm_dir+' does not look okay, exit')
            sys.exit(hhm_dir+' does not look okay, exit')
        else:
            logging.info(hhm_dir+' does not look okay but continue')


    if checkingDirectory(hhr_dir,subfam2nb) :
        logging.info(hhr_dir+' looks okay')
    else:
        if not args.force :
            logging.error(hhr_dir+' does not look okay, exit')
            sys.exit(hhr_dir+' does not look okay, exit')
        else:
            logging.info(hhr_dir+' does not look okay but continue')
    logging.info('done\n')

    ###########################################        
    # collecting the results in each hhr file #
    ###########################################

    logging.info('Collecting the results in each hhr file...')
    hhrHitsList = list()
    for root, dirs, files in os.walk(hhr_dir):
        for filename in files :
            hhr_filename = root+'/'+filename
            hhrHitsList.extend( readingHhrFile(hhr_filename,args.coverage,args.probs) )


    mcl_network_filename = cwd+'/'+'hhblits'+'/'+'hhr.network'
    creatingMclNetworkFile(hhrHitsList,mcl_network_filename)
    logging.info('done\n')
    
    ###############
    # running mcl #
    ###############
    logging.info('running MCL clustering...')
    mcl_log_filename = cwd+'/'+'logs'+'/'+'mcl.log'
    mcl_filename = cwd+'/'+'hhblits'+'/'+'hhr.network.mcl'

    cmd = "/shared/software/bin/mcl "+mcl_network_filename+" --abc -I "+str(args.I)+" -te "+str(args.cpu)+" -o "+mcl_filename+' >'+mcl_log_filename+' 2>&1'
    logging.info(cmd)
    status = os.system(cmd)
    
    if status != 0 :
        logging.error("/shared/software/bin/mcl returned a non 0 status! exit!")
        sys.exit("/shared/software/bin/mcl returned a non 0 status! exit!")

    logging.info('done\n')

    #######################
    # parsing mcl results #
    #######################
    logging.info('parsing MCL clustering results and creating the final output files...')
    subfamily2family = parsingMclFile(mcl_filename,subfam2nb)    



    # network_attr_filename = cwd+'/'+'hhblits'+'/'+'hhr.network.attr'
    # output = open(network_attr_filename,'w')
    # for subfam,fam in subfam2fam.items() :
    #     output.write(subfam+'\t'+'fam'+str(fam)+'\n')
    # output.close()

    
    orf,subfamily,family = orf2familyFunction(cwd , subfamily2family)    
    logging.info(str(orf)+' ORFs were clustered into '+str(subfamily)+' subfamilies and '+str(family)+' families')
    logging.info('done\n')


    if args.fasta :
        os.mkdir(cwd+'/'+'familiesFasta')
        family2subfamilies = defaultdict(list)
        for subfamily,family in subfamily2family.items() :
             family2subfamilies[ family ].append(subfamily)

        for family,liste in family2subfamilies.items() :
            family_filename = cwd+'/'+'familiesFasta'+'/'+family+'.fa'
            output = open(family_filename,'w')
            for subfamily in liste :
                subfamily_filename = cwd+'/'+'subfamiliesFasta'+'/'+subfamily+'.fa'
                for record in SeqIO.parse(subfamily_filename,'fasta') :
                    SeqIO.write(record,output,'fasta')
            output.close()
        
    t2 = datetime.now()
    logging.info('script ended at '+str(t2))
    sys.exit(0)
