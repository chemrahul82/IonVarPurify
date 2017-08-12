"""

Applies the following filters after the variants are filtered
based on Fisher's Exact test & known germline entries

Concept of some of these filters are borrowed from VarScan: 
    Table 1: Koboldt et al., Genome Research 2012
Read position, Strandedness, Distance to 3', Homopolymer,
Map quality difference, Read length difference, MMQS difference

Strandedness & homopolymer length are extracted from TVC generated VCF file
& are extracted in a separate module

"""

from __future__ import division

__Author__ =    "Rahul K. Das"
__Date__ = "May 23, 2016"
__Version_ = "2.0"
__LastModified_ = "Jan 24, 2017"

import os, sys
import numpy as np
import time
from itertools import groupby
from operator import itemgetter
from collections import Counter
from execute_multiprocess import *


def processor2(annoList, bamreadcount, ref_fasta, normal_bam, ncpu):
    
    """ Run bam-readcount on variants & get the QC metrics

    Input: 1.list with fields from annovar annotations table, 2. Tumor bam
    
    Returns: Map quality difference, MMQS difference,
            Read length difference,Read position, Distance to 3'
    
    """

    #create the region file for the first set of filtered variants
    # use annovar's list because it already provides adjusted first and last positions of indels
    tfile = open('bamrc.region', 'w')
    for el in annoList:
        tfile.write('\t'.join(el.split('\t')[:5])+'\n')

    tfile.close()

    # run bam-readcount on the tumor bam file using the region file
    oldbamrc_outfile = '_tmp_/bamrc.out'

    if not os.path.exists(oldbamrc_outfile):
        execute_multiprocess('bamrc.region',normal_bam,ref_fasta,ncpu)
        #command = '%s -q 30 -b 20 -d 500 -l bamrc.region -f %s %s\
        #    > bamrc.out 2>/dev/null' %(bamreadcount,ref_fasta,normal_bam)
        #os.system(command)
        bamrc_outfile = 'bamrc.out'
    else:
        print " bamrc has already been run, delete the bamrc.out file to rerun"
        bamrc_outfile = oldbamrc_outfile


    ref = []; alt = []; fields = []
    diff_q = []; diff_mmqs = []; diff_readlen = []; pos = []; dist3p = []
    alt_depth = []; depth = []; alt_maf = []

    bases = ['A', 'C', 'G', 'T']

    with open('bamrc.region', 'r') as region:
        ref = [line.split('\t')[3] for line in region]
    with open('bamrc.region', 'r') as region:
        alt = [line.split('\t')[4].strip('\n')for line in region]
    with open('bamrc.region', 'r') as region:
        start = [line.split('\t')[1] for line in region]

    
    # get the info from bam-readcount output % store as list
    with open(bamrc_outfile, 'r') as bamread:
        for line in bamread:
            fields.append(line)
    
    #store the positions for all lines in bamrc out
    bamrcpos = [int(line.split('\t')[1]) for line in fields]

    # if a variant start position is very close to the amplicon edge, bam-readcount
    # will not return any values (for snps or short indels) or 
    #return values for positions that are far from the amplicon edge (for long indels), 
    #identify these variants; these should be potential FPs anyways

    missidx= []
    with open('bamrc.region', 'r') as region:
        idlist1 = ['_'.join(line.split('\t')[:2]) for line in region]
    with open('bamrc.region', 'r') as region:
        idlist2 = ['_'.join([line.split('\t')[0],line.split('\t')[2]]) for line in region]
    with open(bamrc_outfile, 'r') as bamread:
        idlist3 = ['_'.join(line.split('\t')[:2]) for line in bamread]

    #expected_vars = set(idlist1) | set(idlist3)
    expected_vars = set(idlist1)
    observed_vars = set(idlist3) #based on start position in bamrc output

    missingvars = list(expected_vars - observed_vars)
    missidx = [idlist1.index(var) for var in missingvars] # based on the start position in bamrc output

    #------------ the main loop -----------

    # extract the relevant metrics from bam-readcount file and
    # append those into the existing info for the stage-1 filtered variants 

    #------------ the main loop -----------

    linecount = 0 # counter for keeping track of variant's start position in bam-readcount output

    for idx in range(len(ref)):
        
        #--------- if amplicon-edge variants------------ 
        # set the metrics values to ridiculously high/low values so these 
        # will be filtered out later

        #both start and end positions don't exist in bamrc output
        if idx in missidx and (idlist2[idx] not in idlist3):
            diff_q.append('1000')
            diff_mmqs.append('1000')
            diff_readlen.append('1000')
            pos.append('-1')
            dist3p.append('-1')
            alt_depth.append('0')
            depth.append('1')
            alt_maf.append('0')
        
        #------first deal with snps-----------------
        elif len(ref[idx]) == 1 and len(alt[idx]) == 1 and\
            ref[idx] != '-' and alt[idx] != '-':
            depth.append(int(fields[linecount].split('\t')[3]))

            for ib, base in enumerate(bases):
                if ref[idx] == base:
                    stats = fields[linecount].split('\t')[ib+5] # 6-9 columns are for A,C,G,T
                    ref_q = float(stats.split(':')[2])
                    ref_mmqs = float(stats.split(':')[9])
                    ref_readlen = float(stats.split(':')[-2])
                if alt[idx] == base:
                    stats = fields[linecount].split('\t')[ib+5]
                    alt_q = float(stats.split(':')[2])
                    pos.append(float(stats.split(':')[7]))
                    alt_mmqs = float(stats.split(':')[9])
                    alt_readlen = float(stats.split(':')[-2])
                    dist3p.append(float(stats.split(':')[-1]))
                    alt_depth.append(int(stats.split(':')[1]))
        
            diff_q.append(ref_q - alt_q)
            diff_mmqs.append(alt_mmqs - ref_mmqs)
            diff_readlen.append(ref_readlen - alt_readlen)
        
            if int(depth[idx])>0:
                alt_maf.append('%3.2f' %(alt_depth[idx]/depth[idx]))
            else:
                alt_maf.append('0.0')
            
            
            #update the line no. to match the start of next variant in the bam-readcount output 
            linecount += 1
            

        #--------------Multi Nucleotide Substitutions---------------
    
        elif len(ref[idx]) > 1 and (len(ref[idx]) == len(alt[idx])):
            refb = list(ref[idx]) # ref. bases spanning
            altb = list(alt[idx]) # alt. bases spanning
            runlen = len(ref[idx]) # length of MNPs

            ##check if bamrc outputs QC for all positions, if not, correct for that
            ## this happens for really long positions where TVC reports MNPs
            #templen = []
            #for k, g in groupby(enumerate(bamrcpos[linecount:linecount+runlen]), lambda (i,x):i-x):
            #    group = map(itemgetter(1), g)
            #    templen.append(len(group))

            ##length of bps where bamrc computes metrics
            #bamrc_runlen = templen[0]

            #if runlen == bamrc_runlen:
            #    refb = list(ref[idx]) # reference bases spanning these positions
            #elif runlen > bamrc_runlen:
            #    refb = list(ref[idx])[:bamrc_runlen] # reference bases spanning these positions
            
            tot = len(fields[linecount].split('\t')) #total tab-seprated fields
        
            depth_i = []; alt_depth_i = []
            ref_q_i = []
            ref_mmqs_i = []
            ref_readlen_i = []
            alt_q_i = []
            alt_mmqs_i = []
            alt_readlen_i = []
            pos_i = []
            dist3p_i = []
            npos = 0

            for jj in range(linecount,linecount+runlen):
                npos += 1
                depth_i.append(int(fields[jj].split('\t')[3]))
                for ib, base in enumerate(bases):
                    if refb[npos-1] == base:
                        stats = fields[jj].split('\t')[ib+5]
                        ref_q_i.append(float(stats.split(':')[2]))
                        ref_mmqs_i.append(float(stats.split(':')[9]))
                        ref_readlen_i.append(float(stats.split(':')[-2]))
                    if altb[npos-1] == base: 
                        stats = fields[jj].split('\t')[ib+5]
                        alt_q_i.append(float(stats.split(':')[2]))
                        pos_i.append(float(stats.split(':')[7]))
                        alt_mmqs_i.append(float(stats.split(':')[9]))
                        alt_readlen_i.append(float(stats.split(':')[-2]))
                        dist3p_i.append(float(stats.split(':')[-1]))
                        alt_depth_i.append(int(stats.split(':')[1]))

            
            depth.append(int(np.mean(depth_i)))
            alt_depth.append(int(np.mean(alt_depth_i)))
            diff_q.append(np.mean(ref_q_i)-np.mean(alt_q_i))
            diff_mmqs.append(np.mean(alt_mmqs_i)-np.mean(ref_mmqs_i))
            diff_readlen.append(np.mean(ref_readlen_i) - np.mean(alt_readlen_i))
            pos.append(np.mean(pos_i))
            dist3p.append(np.mean(dist3p_i))

            if int(depth[idx])>0:
                alt_maf.append('%3.2f' %(alt_depth[idx]/depth[idx]))
            else:
                alt_maf.append('0.0')
            
            linecount += runlen 
        

        #------------------ now deal with indels------------------------

        # Insertions and substitution------
        
        elif ref[idx] == '-' or \
            (len(alt[idx])-len(ref[idx])>0 and ref[idx] != '-'):
            runlen = len(ref[idx]) # length of the actual region 
            
            # get the reference base/bases spanning these positions
            if ref[idx] != '-':
                refb = list(ref[idx]) 
            else:
                refb = [fields[linecount].split('\t')[2]]
            
            ref_q_i = []
            ref_mmqs_i = []
            ref_readlen_i = []
            dist3p_i = []
            npos = 0

            for jj in range(linecount,linecount+runlen):
                npos += 1
            
                for ib, base in enumerate(bases):
                    if refb[npos-1] == base:
                        stats = fields[jj].split('\t')[ib+5].strip('\n')
                        ref_q_i.append(float(stats.split(':')[2]))
                        ref_mmqs_i.append(float(stats.split(':')[9]))
                        ref_readlen_i.append(float(stats.split(':')[-2]))
            
            # stats for insertion reads
            depth.append(int(fields[linecount].split('\t')[3]))
            tot = len(fields[linecount].split('\t')) #total tab-seprated fields

            # get readcounts for all insertions at this position
            # the insertions will be 11th & subsequent fields
            counts = [fields[linecount].split('\t')[j].split(':')[1] for j in range(10,tot)]
            
            if len(counts) > 0:
                # index of the most abundant insertion
                ib = counts.index(max(counts)) + 10
        
                isnstats = fields[linecount].split('\t')[ib].strip('\n')
                alt_q = float(isnstats.split(':')[2])
                pos.append(float(isnstats.split(':')[7]))
                alt_mmqs = float(isnstats.split(':')[9])
                alt_readlen = float(isnstats.split(':')[-2])
                dist3p.append(float(isnstats.split(':')[-1]))
                alt_depth.append(int(isnstats.split(':')[1]))
                
                diff_q.append(np.mean(ref_q_i) - alt_q)
                diff_mmqs.append(alt_mmqs - np.mean(ref_mmqs_i))
                diff_readlen.append(np.mean(ref_readlen_i) - alt_readlen)
            else:
                ###determine the most abundant variant allele when no insertion is detected by bamrc
                allelefields = fields[linecount].strip().split('\t')[5:]
                #index of the reference base at the 1st insertion position
                refidx = [el.split(':')[0] for el in allelefields].index(refb[0])

                ref_q = float(allelefields[refidx].split(':')[2])
                ref_mmqs = float(allelefields[refidx].split(':')[9])
                ref_readlen = float(allelefields[refidx].split(':')[-2])
                #ref_pos = float(allelefields[refidx].split(':')[7])

                #remove ref. fields
                allelefields.pop(refidx)
                #max depth of non-reference alleles
                maxdepth =  max([int(el.split(':')[1]) for el in allelefields])
                if int(maxdepth) > 0:
                    #variant allele (SNP) detected at this position
                    maxalleleidx = [int(el.split(':')[1]) for el in allelefields].index(maxdepth)

                    alt_depth.append(int(allelefields[maxalleleidx].split(':')[1]))
                    alt_q = float(allelefields[maxalleleidx].split(':')[2])
                    alt_mmqs = float(allelefields[maxalleleidx].split(':')[9])
                    alt_readlen = float(allelefields[maxalleleidx].split(':')[-2])
                    pos.append(float(allelefields[maxalleleidx].split(':')[7]))
                    dist3p.append(float(allelefields[maxalleleidx].split(':')[-1]))

                    diff_q.append('%3.2f' %(ref_q - alt_q))
                    diff_mmqs.append('%3.2f' %(alt_mmqs - ref_mmqs))
                    diff_readlen.append('%3.2f' %(ref_readlen - alt_readlen))
                else:
                    #no variant allele (SNP) detected at this position; set the metrics to '.'
                    diff_q.append('.')
                    diff_mmqs.append('.')
                    diff_readlen.append('.')
                    pos.append('.')
                    dist3p.append('.')
                    alt_depth.append('0')

        
            if int(depth[idx])>0:
                alt_maf.append('%4.3f' %(int(alt_depth[idx])/depth[idx]))
            else:
                alt_maf.append('0.0')
            
            linecount += runlen

        
        # Deletions---------    
        elif (alt[idx] == '-') or \
                (len(ref[idx]) > len(alt[idx]) and alt[idx] != '-'):
            runlen = len(ref[idx]) # length of the actual region
            
            #check if bamrc outputs QC for all positions, if not, correct for that
            # this happens for really long positions where TVC reports deletion
            templen = []
            for k, g in groupby(enumerate(bamrcpos[linecount:linecount+runlen]), lambda (i,x):i-x):
                group = map(itemgetter(1), g)
                templen.append(len(group))

            #length of bps where bamrc computes metrics
            bamrc_runlen = templen[0]

            if runlen == bamrc_runlen:
                refb = list(ref[idx]) # reference bases spanning these positions
            elif runlen > bamrc_runlen:
                refb = list(ref[idx])[:bamrc_runlen] # reference bases spanning these positions
           
            runlen = bamrc_runlen
            
            ref_q_i = []
            ref_mmqs_i = []
            ref_readlen_i = []
            dist3p_i = []
            npos = 0

            for jj in range(linecount,linecount+runlen):
                npos += 1
            
                for ib, base in enumerate(bases):
                    if refb[npos-1] == base:
                        stats = fields[jj].split('\t')[ib+5].strip('\n')
                        ref_q_i.append(float(stats.split(':')[2]))
                        ref_mmqs_i.append(float(stats.split(':')[9]))
                        ref_readlen_i.append(float(stats.split(':')[-2]))
            
            # stats for deletion reads
            # get readcounts for all deletions at the first position
            depth.append(int(fields[linecount].split('\t')[3]))
            tot = len(fields[linecount].split('\t')) #total tab-seprated fields
            counts = [fields[linecount].split('\t')[j].split(':')[1] for j in range(10,tot)]
            
            if len(counts) >0:
                # index of the most abundant deletion
                ib = counts.index(max(counts)) + 10
                    
                delstats = fields[linecount].split('\t')[ib].strip('\n')
                alt_q = float(delstats.split(':')[2])
                pos.append(float(delstats.split(':')[7]))
                alt_mmqs = float(delstats.split(':')[9])
                alt_readlen = float(delstats.split(':')[-2])
                dist3p.append(float(delstats.split(':')[-1]))
                alt_depth.append(int(delstats.split(':')[1]))
                
                diff_q.append(np.mean(ref_q_i) - alt_q)
                diff_mmqs.append(alt_mmqs - np.mean(ref_mmqs_i))
                diff_readlen.append(np.mean(ref_readlen_i) - alt_readlen)
            else:
                ###determine the most abundant variant allele when no deletion is detected by bamrc
                allelefields = fields[linecount].strip().split('\t')[5:]
                #index of the reference base at the 1st deletion position
                refidx = [el.split(':')[0] for el in allelefields].index(refb[0])

                ref_q = float(allelefields[refidx].split(':')[2])
                ref_mmqs = float(allelefields[refidx].split(':')[9])
                ref_readlen = float(allelefields[refidx].split(':')[-2])
                #ref_pos = float(allelefields[refidx].split(':')[7])
                

                #remove ref. fields
                allelefields.pop(refidx)
                #max depth of non-reference alleles
                maxdepth =  max([int(el.split(':')[1]) for el in allelefields])
                if int(maxdepth) > 0:
                    #variant allele (SNP) detected at this position
                    maxalleleidx = [int(el.split(':')[1]) for el in allelefields].index(maxdepth)

                    alt_depth.append(int(allelefields[maxalleleidx].split(':')[1]))
                    alt_q = float(allelefields[maxalleleidx].split(':')[2])
                    alt_mmqs = float(allelefields[maxalleleidx].split(':')[9])
                    alt_readlen = float(allelefields[maxalleleidx].split(':')[-2])
                    pos.append(float(allelefields[maxalleleidx].split(':')[7]))
                    dist3p.append(float(allelefields[maxalleleidx].split(':')[-1]))

                    diff_q.append('%3.2f' %(ref_q - alt_q))
                    diff_mmqs.append('%3.2f' %(alt_mmqs - ref_mmqs))
                    diff_readlen.append('%3.2f' %(ref_readlen - alt_readlen))
                else:
                    #no variant allele (SNP) detected at this position; set the metrics to '.'
                    diff_q.append('.')
                    diff_mmqs.append('.')
                    diff_readlen.append('.')
                    pos.append('.')
                    dist3p.append('.')
                    alt_depth.append('0')
            
            
            if int(depth[idx])>0:
                alt_maf.append('%4.3f' %(int(alt_depth[idx])/depth[idx]))
            else:
                alt_maf.append('0.0')
            
            linecount += runlen

    #clean up temp files
    #os.remove('temp.out')
    #os.remove('temp.region')

    return diff_q, diff_mmqs, diff_readlen, pos, dist3p, alt_maf, alt_depth, depth
