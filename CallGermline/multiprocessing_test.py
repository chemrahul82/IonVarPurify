from __future__ import division
import os
import math
import multiprocessing


def split_file(nchunk,regionFile):
    #regionFile = '/rawdata/Rahul/Development/Somatic_Calling_Pipeline/Dev/version2.0/CallSomatic/testing1/Somatic_Variants_UHQ/_tmp_/bamrc.region' 
    with open(regionFile,'r') as f:
        nline = len([line for line in f])

    nline_chunk = int(math.floor(nline/nchunk))

    command = 'split -l %s %s' %(nline_chunk,startFile)
    os.system(command)

def bamrc(regionFile,bam,refFasta):
    bam = "/mnt/Orcus/projects/LungBio/pair_035/Tumor/Merged/PTRIM.bam"
    refFasta = "/results/referenceLibrary/tmap-f3/hg19/hg19.fasta"

    bname = os.path.basename(regionFile)
    command = 'bam-readcount -q 30 -b 20 -f %s -d 2000 -l %s %s > out.%s 2>/dev/null' %(refFasta,regionFile,bam,bname)
    os.system(command)

    
def main(regionFile,ncpu,bam,refFasta):
    
	split_file(ncpu,regionFile)
    
    chunks = sorted(i for i in os.listdir(os.getcwd()) if 'xa' in i)
    jobs = []
    for i in chunks:
        p = multiprocessing.Process(target=bamrc,args=(i,))
        jobs.append(p)
        p.start()


if __name__ == "__main__":
	regionFile = '/rawdata/Rahul/Development/Somatic_Calling_Pipeline/Dev/version2.0/CallSomatic/testing1/Somatic_Variants_UHQ/_tmp_/bamrc.region' 
	bam = "/mnt/Orcus/projects/LungBio/pair_035/Tumor/Merged/PTRIM.bam"
	fasta = "/results/referenceLibrary/tmap-f3/hg19/hg19.fasta"
	ncpu=4
	
	main(regionFile,ncpu,bam,refFasta)

