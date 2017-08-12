"""
Extract sequence context/motif around a given position
"""

__Author__ =    "Rahul K. Das"
__LastModified__ = "Nov 7, 2016"

import os,sys
import subprocess

def getSequence(varFile,refFasta,lenMotif):
    """
    @params varFile: file with list of variants with contig name & position: chr1,1000000 (+1 positions & tab-delimited)
    @params refFasta: reference fasta file
    @params lenMotif: Integer::length of desired motif in one side of the position
    @returns motifList: list with sequence motif
    """
    
    command = '''awk '{print $1":"($2-%s)"-"($2+%s)}' %s | xargs samtools faidx %s | grep -v "chr"''' %(
            lenMotif,
            lenMotif,
            varFile,
            refFasta)

    motifSeq = subprocess.check_output(command, shell=True).strip()

    motifList = motifSeq.split('\n')
    return motifList

#if __name__ == "__main__":
#
#    varFile = sys.argv[1]
#    refFasta = sys.argv[2]
#    lenMotif = 5
#    getSequence(varFile,refFasta,lenMotif)


