#! /usr/bin/env python

import os,sys,re
from collections import defaultdict

output_directory = sys.argv[1].rsplit('/',1)[0]

cmd = 'cd '+output_directory+'/'+'a3m'+' && '+'ffindex_build -as '+'..'+'/'+'db_a3m.ffdata '+'..'+'/'+'db_a3m.ffindex '+'*.a3m'
print(cmd)
status = os.system(cmd)
if status != 0 :
    sys.exit(status)
cmd = 'cd '+output_directory+'/'+'hhm'+' && '+'ffindex_build -as '+'..'+'/'+'db_hhm.ffdata '+'..'+'/'+'db_hhm.ffindex '+'*.hhm'
print(cmd)
status = os.system(cmd)
if status != 0 :
    sys.exit(status)

cmd = 'cstranslate -A $HHLIB/data/cs219.lib -D $HHLIB/data/context_data.lib -I a3m -x 0.3 -c 4'+' -f -i '+output_directory+'/'+'db_a3m'+' -o '+output_directory+'/'+'db_cs219'+' -b'
print(cmd)
status = os.system(cmd)
if status != 0 :
    sys.exit(status)

