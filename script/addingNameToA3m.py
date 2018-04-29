#! /usr/bin/env python

import sys,os

input_filename = sys.argv[1]
output_filename = sys.argv[2]

name = os.path.basename(input_filename).replace('.a3m','')

output = open(output_filename,'w')
output.write('#'+name+'\n')
    
file=open(input_filename,'r')
for line in file :
    output.write(line)
file.close()
output.close()
