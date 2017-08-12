"""

calculate estimated readlengths based on amplicon lengths
spanning a variant

"""

from __future__ import division

__Author__ =    "Rahul K. Das"
__Date__ = "May 23, 2016"
__Version_ = "2.0"
__LastModified_ = "Feb 14, 2017"

import os, sys, time
import numpy as np
from itertools import groupby
from operator import itemgetter
from collections import Counter


def expected_readlen(bedfile, bedHeader, contig, start):
    
    """Calculates the expected length of the reads supporting a variant from the length of amplicon(s)
    supporting it
    @param: the project bed file, 
    @param: header present in project bed? 
    @param: the annovar list for variants
    @return: the expected read length for a variant based on supporting amplicon(s)
    """

    print contig,start
    #chromosomes from project bed
    with open(bedfile, 'r') as f:
        if bedHeader=='yes':
            chroms = [key for key in dict(Counter([line.split('\t')[0] for il, line in enumerate(f) if il>0])).keys()]
        else:
            chroms = [key for key in dict(Counter([line.split('\t')[0] for il, line in enumerate(f)])).keys()]

    #dictionaries: keys are chrom and values are lists with start/end positions
    ampStarts = {}; ampEnds = {}
    for chrom in chroms:
        ampStarts[chrom] = []
        ampEnds[chrom] = []

    with open(bedfile, 'r') as f:
        if bedHeader=='yes':
            next(f)
        for line in f:
            ampStarts[line.split('\t')[0]].append(int(line.split('\t')[1]))

    with open(bedfile, 'r') as f:
        if bedHeader=='yes':
            next(f)    
        for line in f:
            ampEnds[line.split('\t')[0]].append(int(line.split('\t')[2]))

    
    #exp_readlength = []
    #of = open('expected_readlen.txt','w')
        
    nVarAmp = 0
    amplength = []
    #no. of amplicons spanning the chromosome
    namp = len(ampStarts.get(contig))
    for ia in range(0,namp):
        #print int(ampStarts.get(chrom)[ia]),int(ampEnds.get(chrom)[ia])
        #no. of amplicons that span the variant
        if (int(start) >= int(ampStarts.get(contig)[ia])) and (int(start) <= int(ampEnds.get(contig)[ia])):
            nVarAmp += 1
            amplength.append(int(ampEnds.get(contig)[ia]) - int(ampStarts.get(contig)[ia]) + 1)
        # if amplicon start is greater than variant pos, go to next variant
        if int(start) < int(ampStarts.get(contig)[ia]):
            break
        
    exp_readlength = np.mean(amplength)
    #exp_readlength.append(np.mean(amplength))

    #of.write(str(np.mean(amplength)) + '\n')
    
    #of.close()
   
    return exp_readlength


    ## save the chrom#, begin, end of amplicons from the project bed file
    #with open(bedfile,'r') as bf:
    #    bedpos = [line.split('\t')[:3] for line in bf if 'chr' in line.split('\t')[0]]
    #
    ## iterate through the variants
    #exp_readlen = []
    #for idx, var in enumerate(annoList):
    #    varchrom = var.split('\t')[0]
    #    varpos = var.split('\t')[1]
    #    namp = 0
    #    amplen = []

    #    # iterate through the amplicons
    #    for amp in bedpos:
    #        ampchrom = amp[0]
    #        ampstart = amp[1]
    #        ampend = amp[2]

    #        # find the amplicons that belong to the chrom and pos of the variant
    #        if (varchrom == ampchrom) and \
    #            (int(varpos) >= int(ampstart)  and int(varpos) <= int(ampend)):
    #            amplen.append(int(ampend) - int(ampstart) + 1)
    #            namp += 1
    #        # if amplicon start is greater than variant pos in same chrom, go to next variant
    #        if (varchrom == ampchrom) and (int(varpos) < int(ampstart)):
    #            break
    #        # if amplicon/s found and chrom is different, go to next variant
    #        if (varchrom != ampchrom) and (namp > 0):
    #            break
    #        
    #    #print varchrom, varpos, namp 
    #    exp_readlen.append(np.mean(amplen))


