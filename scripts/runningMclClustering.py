#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import argparse
import logging
from datetime import date, time, datetime
import shutil
import json

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

    return subfam2fam,cpt

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
            if probs > probs_threshold and ( scover > coverage_threshold and qcover > coverage_threshold ) :
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
    parser.add_argument('config_filename', help='the path of the FASTA_FILENAME that contains the proteins sequences to cluster')
    parser.add_argument('--coverage',type=float,default=0.75,help='coverage threshold (default: 0.75)')
    parser.add_argument('--probs',type=float,default=0.95,help='probability threshold (default: 0.95)')
    parser.add_argument('--force',action='store_true',default=False,help='force MCL clustering (default: False)')

    args = parser.parse_args()


    # checking arguments
    if not os.path.exists(args.config_filename) :
        sys.exit(args.config_filename+' does not exist, exit')
    else:
        config_filename = os.path.abspath(args.config_filename)


    with open(config_filename) as f:
        data = json.load(f)


    # checking if hhblits runned fine 
    cwd = data['directory']
    print(cwd)
    subfam2nb = data['clusters']
    print(len(subfam2nb))

    fasta_dir = cwd+'/'+'subfamiliesFasta'
    hhm_dir = cwd+'/'+'hhblits'+'/'+'hhm'
    hhr_dir = cwd+'/'+'hhblits'+'/'+'hhr'

    if checkingDirectory(hhm_dir,subfam2nb) :
        print(hhm_dir+' looks okay')
    else:
        if not args.force :
            sys.exit(hhm_dir+' does not look okay')
        else:
            print(hhm_dir+' does not look okay but continue')


    if checkingDirectory(hhr_dir,subfam2nb) :
        print(hhr_dir+' looks okay')
    else:
        if not args.force :
            sys.exit(hhr_dir+' does not look okay')
        else:
            print(hhr_dir+' does not look okay but continue')

    # collecting the results in each hhr file
    hhrHitsList = list()
    for root, dirs, files in os.walk(hhr_dir):
        for filename in files :
            hhr_filename = root+'/'+filename
            hhrHitsList.extend( readingHhrFile(hhr_filename,args.coverage,args.probs) )

    output_network_filename = cwd+'/'+'hhblits'+'/'+'hhr.network'
    output = open(output_network_filename,'w')
    for hit in sorted(hhrHitsList,key= lambda x:x[0]) :
        output.write('\t'.join(str(elt) for elt in hit[0:3])+'\n')
    output.close()


    # running mcl
    mcl_log_filename = cwd+'/'+'logs'+'/'+'mcl.log'
    mcl_filename = cwd+'/'+'hhblits'+'/'+'hhr.network.mcl'
    cmd = "/usr/bin/mcl "+output_network_filename+" --abc -I 2.0 -o "+mcl_filename+' >'+mcl_log_filename+' 2>&1'
    print(cmd)
    status = os.system(cmd)
    if status != 0 :
        sys.exit("/usr/bin/mcl has a non 0 status! exit!")


    # parsing mcl results
    subfam2fam,cpt = parsingMclFile(mcl_filename,subfam2nb)    
    print(len(subfam2fam))
    print(cpt)

    network_attr_filename = cwd+'/'+'hhblits'+'/'+'hhr.network.attr'
    output = open(network_attr_filename,'w')
    for subfam,fam in subfam2fam.items() :
        output.write(subfam+'\t'+'fam'+str(fam)+'\n')
    output.close()
