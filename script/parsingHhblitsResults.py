#! /usr/bin/env python

import os,sys,re
from collections import defaultdict


def readingOutput(hhr_filename) :
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

            target = liste[1]
            prob = liste[2]
            qcoord = liste[8]
            qstart = float(qcoord.split('-')[0])
            qend = float(qcoord.split('-')[1])            
            qcover = ( qend - qstart + 1.0 ) / qlen
            
            scoord = liste[9]
            sstart = float(liste[9].split('-')[0])
            send = float(liste[9].split('-')[1])            
            slen = float(liste[10].replace('(','').replace(')',''))
            scover = ( send - sstart + 1 ) / slen
            result.append(query+'\t'+target+'\t'+prob+'\t'+str(qcover)+'\t'+str(scover))
    file.close()
    return result




output_filename = sys.argv[1]

output_directory = sys.argv[1].rsplit('/',1)[0]

output = open(output_filename,'w')
output.write('query'+'\t'+'subject'+'\t'+'prob'+'\t'+'qcover'+'\t'+'scover'+'\n')

for (path, dirs, files) in os.walk(output_directory+'/'+'hhr'):
    for filename in files :
        hhr_filename = path+'/'+filename
        print(filename)
        for line in readingOutput(hhr_filename):
            output.write(line+'\n')
            
output.close()
