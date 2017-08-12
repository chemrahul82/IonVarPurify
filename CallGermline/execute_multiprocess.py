"""
Execute bamreadcount in multiprocessing mode
"""

from __future__ import division

__Author__ =    "Rahul K. Das"
__Version__ = "2.0"
__LastUpdated__ = "Jan 31, 2017"


import os
import math
import multiprocessing

def split_file(nchunk,regionFile):
    with open(regionFile,'r') as f:
        nline = len([line for line in f])

    nline_chunk = int(math.floor(nline/nchunk))

    #clean up from priviously unfinished/crashed jobs
    if os.path.exists('out.xaa'):
        os.system('rm out.xa*')
    if os.path.exists('xaa'):
        os.system('rm xa*')

    command = 'split -l %s %s' %(nline_chunk,regionFile)
    os.system(command)


def bamrc(regionFile,bam,refFasta):

    bname = os.path.basename(regionFile)
    command = 'bam-readcount -q 30 -b 20 -f %s -d 500 -l %s %s > out.%s 2>/dev/null' %(refFasta,regionFile,bam,bname)
    os.system(command)

    
def execute_multiprocess(regionFile,bam,refFasta,ncpu):
    
    split_file(ncpu,regionFile)
    
    chunks = sorted(i for i in os.listdir(os.getcwd()) if 'xa' in i)
    jobs = []
    for i in chunks:
        p = multiprocessing.Process(target=bamrc,args=(i,bam,refFasta))
        jobs.append(p)
        p.start()
        
    for j in jobs:
        j.join() 

    os.system('cat out.xa* > bamrc.out')
    os.system('rm out.xa* xa*')


