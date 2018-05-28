#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import argparse
import logging
from datetime import date, time, datetime
import shutil
import json


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
                # if query < target :
                #     result.append([query,target,probs,qcover,scover])
                # else:
                #     result.append([target,query,probs,qcover,scover])
                result.append([query,target,probs,qcover,scover])
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
#    parser.add_argument('output_filename', help='the path of the OUTPUT_FILENAME that contains the result of the MCL clustering')
    parser.add_argument('--coverage',type=float,default=0.75,help='number of CPUs used by mmseqs (default: 1)')
    parser.add_argument('--probs',type=float,default=0.95,help='minimal size of the protein families to retain (default: 2)')

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
        sys.exit(hhm_dir+' does not look okay')


    if checkingDirectory(hhr_dir,subfam2nb) :
        print(hhr_dir+' looks okay')
    else:
        sys.exit(hhr_dir+' does not look okay')

    # collecting the results in each hhr file
    hhrHitsList = list()
    for root, dirs, files in os.walk(hhr_dir):
        for filename in files :
            hhr_filename = root+'/'+filename
            hhrHitsList.extend( readingHhrFile(hhr_filename,args.coverage,args.probs) )

    output_network_filename = cwd+'/'+'hhblits'+'/'+'hhr.network'
    output = open(output_network_filename,'w')
    for hit in hhrHitsList :
        output.write('\t'.join(str(elt) for elt in hit[0:3])+'\n')
    output.close()


    # running mcl
    output_mcl_filename = cwd+'/'+'hhblits'+'/'+'hhr.network.mcl'
    cmd = "mcl "+output_network_filename+" --abc -I 2.0 -o "+output_mcl_filename
    print(cmd)
    status = os.system(cmd)
