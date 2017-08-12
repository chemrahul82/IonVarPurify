"""

Applies the following filters after the variants are filtered
based on Fisher's Exact test & known germline entries

Concept of some of these filters are borrowed from VarScan: 
    Table 1: Koboldt et al., Genome Research 2012
Read position, Strandedness, Distance to 3', Homopolymer,
Map quality difference, Read length difference, MMQS difference
Also applying:
1. Difference in length between the expected value and the obtained value for reference supporting reads
2. Re/Alt Difference in read position

All of the above are obtained by runnin bam-readcount from WashU

"""

from __future__ import division

__Author__ =    "Rahul K. Das"
__Date__ = "May 23, 2016"
__Version_ = "2.0"
__LastModified_ = "Feb 22, 2017"

import os, sys, time
import numpy as np
import time
from itertools import groupby
from operator import itemgetter
from collections import Counter
from execute_multiprocess import *



def processor2(readq, baseq, bamrcPath, ref_fasta, tumor_bam, bedfile, bedHeader, varList, annoList, ncpu):
    
    """ Run bam-readcount on variants & get the QC metrics

    @param: read quality threshold
    @param: base quality threshold
    @param: bam-readcount executable path
    @param: path for the reference fasta file
    @param: path for the tumor bam file
    @param: path for the project bed file
    @param: header present in project bed? 
    @param: list with fields from vcf file
    @param: list with fields from annovar annotations table
    @param: number of CPU to run bamreadcount
    
    Returns: 1. Ref/Alt Map quality difference
             2. Alt/Ref MMQS difference
             3. Ref/Alt Read length difference
             4. Read position
             5. Ref/Alt read position difference
             6. Distance to 3'
             7. Alt read length
             8. Expected/Ref read length difference
             9. read-Q threshold adjusted AF of the variant
             10. read-Q threshold adjusted depth of the variant
             11. read-Q threshold adjusted total depth 
    
    """

    #create the region file for the first set of filtered variants
    # use annovar's list because it already provides adjusted first and last positions of indels
    
    tfile = open('bamrc.region', 'w')
    tvc_maf = []; tvc_depth = []
    for el in annoList:
        tfile.write('\t'.join(el.split('\t')[:5])+'\n')
        
    
    # store the TVC-reported MAF and depth for the variant in the T
    # this will be used for those TVC variants that don't match with bam-readcount
    for el in varList:
        tvc_maf.append(float(el.split('\t')[9])) 
        tvc_depth.append(int(el.split('\t')[10]))
    
    tfile.close()


    # get the expected length of reads supporting the variant
    old_exprlenFile = '_tmp_/expected_readlen.txt'
    if not os.path.exists(old_exprlenFile):
        print " Calculating estimated read length around the variant positions..."
        exp_readlen = expected_readlen(bedfile, bedHeader, annoList)
    else:
        print " Expected read lengths already have been calculated, delete the old file _tmp_/expected_readlen.txt to recompute"
        with open(old_exprlenFile,'r') as rf:
            exp_readlen = [float(line.strip()) for line in rf] 


    # run bam-readcount on the tumor bam file using the region file
    oldbamrc_outfile = '_tmp_/bamrc.out'
    
    if not os.path.exists(oldbamrc_outfile):
        execute_multiprocess(bamrcPath,'bamrc.region',tumor_bam,ref_fasta,ncpu)
        
        #command = '%s -q %s -b %s -d 2000 -l bamrc.region -f %s %s\
        #        > bamrc.out 2>/dev/null' %(bamreadcount, readq, baseq, ref_fasta, tumor_bam)
        #os.system(command)
        
        bamrc_outfile = 'bamrc.out'
    else:
        print " bamrc has already been run, delete the bamrc.out file to rerun"
        bamrc_outfile = oldbamrc_outfile


    ref = []; alt = []; fields = []; diff_q = []; diff_mmqs = []; diff_readlen = []; diffpos = []; altpos = []; dist3p = []; 
    alt_readlen_ = []; diff_expec_readlen = []; alt_depth = []; depth = []; alt_maf = []

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


    # if a variant position is very close to the amplicon edge, bam-readcount
    # will not return any values, identify these variants; these should be potential
    # FPs anyways

    missidx= []
    with open('bamrc.region', 'r') as region:
        idlist1 = ['_'.join(line.split('\t')[:2]) for line in region]
    with open(bamrc_outfile, 'r') as bamread:
        idlist2 = ['_'.join(line.split('\t')[:2]) for line in bamread]

    missingvars = list(set(idlist1) - set(idlist2))
    missidx = [idlist1.index(var) for var in missingvars]

    #------------ the main loop -----------

    # extract the relevant metrics from bam-readcount file and
    # append those into the existing info for the stage-1 filtered variants 

    #------------ the main loop -----------

    linecount = 0 # counter for keeping track of variant's start position in bam-readcount output

    for idx in range(len(ref)):
        #print idx,linecount 
        #--------- if amplicon-edge variants------------ 
        # set the metrics values to ridiculously high/low values so these 
        # will be filtered out later

        if idx in missidx:
            diff_q.append('1000')
            diff_mmqs.append('1000')
            diff_readlen.append('1000')
            altpos.append('0')
            diffpos.append('1')
            dist3p.append('1')
            alt_readlen_.append('0')
            alt_depth.append('0')
            depth.append('1')
            alt_maf.append('0')
            diff_expec_readlen.append('100')            

        
        #-------------------------------------------------------
        #------------------first deal with snps-----------------
        #-------------------------------------------------------

        elif len(ref[idx]) == 1 and len(alt[idx]) == 1 and\
            ref[idx] != '-' and alt[idx] != '-':
            depth.append(int(fields[linecount].split('\t')[3]))

            for ib, base in enumerate(bases):
                if ref[idx] == base:
                    stats = fields[linecount].split('\t')[ib+5] # 6-9 columns are for A,C,G,T
                    ref_q = float(stats.split(':')[2])
                    ref_mmqs = float(stats.split(':')[9])
                    ref_readlen = float(stats.split(':')[-2])
                    refpos = float(stats.split(':')[7])
                if alt[idx] == base:
                    stats = fields[linecount].split('\t')[ib+5]
                    alt_q = float(stats.split(':')[2])
                    altpos.append(float(stats.split(':')[7]))
                    alt_mmqs = float(stats.split(':')[9])
                    alt_readlen = float(stats.split(':')[-2])
                    dist3p.append(float(stats.split(':')[-1]))
                    alt_depth.append(int(stats.split(':')[1]))
        
            diff_q.append(ref_q - alt_q)
            diff_mmqs.append(alt_mmqs - ref_mmqs)
            diff_readlen.append(ref_readlen - alt_readlen)
            diffpos.append(refpos - altpos[idx])
            alt_readlen_.append(alt_readlen)
            diff_expec_readlen.append('%6.2f' %(float(exp_readlen[idx])-alt_readlen))
            #diff_expec_readlen.append('%3.2f' %(float(exp_readlen[idx])-ref_readlen))
       
            if int(depth[idx])>0:
                alt_maf.append('%5.4f' %(alt_depth[idx]/depth[idx]))
            else:
                alt_maf.append('0.0')


            #update the line no. to match the start of next variant in the bam-readcount output 
            linecount += 1
        

        #-------------------------------------------------------
        #--------------Multi Nucleotide Substitutions-----------
        #-------------------------------------------------------
    
        elif len(ref[idx]) > 1 and (len(ref[idx]) == len(alt[idx])):
            refb = list(ref[idx]) # ref. bases spanning
            altb = list(alt[idx]) # alt. bases spanning
            runlen = len(ref[idx]) # length of MNPs

            tot = len(fields[linecount].split('\t')) #total tab-seprated fields
        
            depth_i = []; ref_q_i = []; ref_mmqs_i = []; ref_readlen_i = []; refpos_i = []
            alt_q_i = []; alt_mmqs_i = []; alt_readlen_i = []; altpos_i = []; dist3p_i = []; alt_depth_i = []
            npos = 0

            for jj in range(linecount,linecount+runlen):
                npos += 1
                #only get the stats for positions that are mismatches
                if refb[npos-1] != altb[npos-1]:
                    depth_i.append(int(fields[jj].split('\t')[3]))
                    for ib, base in enumerate(bases):
                        if refb[npos-1] == base:
                            stats = fields[jj].split('\t')[ib+5]
                            ref_q_i.append(float(stats.split(':')[2]))
                            ref_mmqs_i.append(float(stats.split(':')[9]))
                            ref_readlen_i.append(float(stats.split(':')[-2]))
                            refpos_i.append(float(stats.split(':')[7]))
                        if altb[npos-1] == base: 
                            stats = fields[jj].split('\t')[ib+5]
                            alt_q_i.append(float(stats.split(':')[2]))
                            altpos_i.append(float(stats.split(':')[7]))
                            alt_mmqs_i.append(float(stats.split(':')[9]))
                            alt_readlen_i.append(float(stats.split(':')[-2]))
                            dist3p_i.append(float(stats.split(':')[-1]))
                            alt_depth_i.append(int(stats.split(':')[1]))
        
            depth.append(int(np.mean(depth_i)))
            diff_q.append('%6.2f' %(np.mean(ref_q_i)-np.mean(alt_q_i)))
            diff_mmqs.append('%6.2f' %(np.mean(alt_mmqs_i)-np.mean(ref_mmqs_i)))
            diff_readlen.append('%6.2f' %(np.mean(ref_readlen_i) - np.mean(alt_readlen_i)))
            altpos.append('%3.2f' %np.mean(altpos_i))
            dist3p.append('%3.2f' %np.mean(dist3p_i))
            diffpos.append('%3.2f' %(np.mean(refpos_i) - np.mean(altpos_i)))
            alt_readlen_.append(np.mean(alt_readlen_i))
            diff_expec_readlen.append('%6.2f' %(float(exp_readlen[idx])-np.mean(alt_readlen_i)))
            #diff_expec_readlen.append('%3.2f' %(float(exp_readlen[idx])-np.mean(ref_readlen_i)))
            alt_depth.append(int(np.mean(alt_depth_i)))

            if int(depth[idx])>0:
                alt_maf.append('%5.4f' %(alt_depth[idx]/depth[idx]))
            else:
                alt_maf.append('0.0')

            linecount += runlen 
        


        #-------------------------------------------------------
        #--------------- Insertions and substitution------------
        #-------------------------------------------------------
        
        elif ref[idx] == '-' or \
            (len(alt[idx])-len(ref[idx])>0 and ref[idx] != '-'):

            runlen = len(ref[idx]) # length of the actual region 
            
            # get the reference base/bases spanning these positions
            if ref[idx] != '-':
                refb = list(ref[idx]) 
            else:
                refb = [fields[linecount].split('\t')[2]]
            
            ref_q_i = []; ref_mmqs_i = []; ref_readlen_i = []; refpos_i = []
            npos = 0

            for jj in range(linecount,linecount+runlen):
                npos += 1
            
                for ib, base in enumerate(bases):
                    if refb[npos-1] == base:
                        stats = fields[jj].split('\t')[ib+5].strip('\n')
                        ref_q_i.append(float(stats.split(':')[2]))
                        ref_mmqs_i.append(float(stats.split(':')[9]))
                        ref_readlen_i.append(float(stats.split(':')[-2]))
                        refpos_i.append(float(stats.split(':')[7]))
            
            # stats for insertion reads
            depth.append(int(fields[linecount].split('\t')[3]))
            tot = len(fields[linecount].split('\t')) #total tab-seprated fields

            # get readcounts for all insertions at this position
            # the insertions will be 11th & subsequent fields
            counts = [int(fields[linecount].split('\t')[j].split(':')[1]) for j in range(10,tot)]
            
            if len(counts) > 0:
                # index of the most abundant insertion
                ib = counts.index(max(counts)) + 10
        
                isnstats = fields[linecount].split('\t')[ib].strip('\n')
                alt_q = float(isnstats.split(':')[2])
                altpos.append(float(isnstats.split(':')[7]))
                alt_mmqs = float(isnstats.split(':')[9])
                alt_readlen = float(isnstats.split(':')[-2])
                dist3p.append(float(isnstats.split(':')[-1]))
                alt_depth.append(int(isnstats.split(':')[1]))
                
                diff_q.append('%6.2f' %(np.mean(ref_q_i) - alt_q))
                diff_mmqs.append('%6.2f' %(alt_mmqs - np.mean(ref_mmqs_i)))
                diffpos.append('%3.2f' %(np.mean(refpos_i) - altpos[idx]))
                diff_readlen.append('%6.2f' %(np.mean(ref_readlen_i) - alt_readlen))
                alt_readlen_.append(alt_readlen)
                diff_expec_readlen.append('%6.2f' %(float(exp_readlen[idx]) - alt_readlen))
                #diff_expec_readlen.append('%3.2f' %(float(exp_readlen[idx])-np.mean(ref_readlen_i)))

                #for insertion total depth will be total depth + insertion depth
                depth[idx] = depth[idx] + alt_depth[idx]
                if int(depth[idx])>0:
                    alt_maf.append('%5.4f' %(alt_depth[idx]/depth[idx]))
                else:
                    alt_maf.append('0.0')

            else:
                ###determine the most abundant variant allele when no insertion is detected by bamrc
                allelefields = fields[linecount].strip().split('\t')[5:]
                #index of the reference base at the 1st insertion position
                refidx = [el.split(':')[0] for el in allelefields].index(refb[0])

                ref_q = float(allelefields[refidx].split(':')[2])
                ref_mmqs = float(allelefields[refidx].split(':')[9])
                ref_readlen = float(allelefields[refidx].split(':')[-2])
                ref_pos = float(allelefields[refidx].split(':')[7])

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
                    altpos.append(float(allelefields[maxalleleidx].split(':')[7]))
                    dist3p.append(float(allelefields[maxalleleidx].split(':')[-1]))

                    #alt_q_all.append(alt_q)
                    diff_q.append('%6.2f' %(ref_q - alt_q))
                    diff_mmqs.append('%6.2f' %(alt_mmqs - ref_mmqs))
                    diffpos.append('%3.2f' %(ref_pos - altpos[idx]))
                    diff_readlen.append('%6.2f' %(ref_readlen - alt_readlen))
                    diff_expec_readlen.append('%6.2f' %(float(exp_readlen[idx])-alt_readlen))
                    #diff_expec_readlen.append('%3.2f' %(float(exp_readlen[idx])-ref_readlen))
                else:
                    #no variant allele (SNP) detected at this position; set the metrics to '.'
                    altpos.append('.')
                    alt_depth.append('0')
                    diff_q.append('.')
                    diff_mmqs.append('.')
                    diff_readlen.append('.')
                    diffpos.append('.')
                    dist3p.append('.')
                    diff_expec_readlen.append('.')
                    #diff_expec_readlen.append('%3.2f' %(float(exp_readlen[idx])-np.mean(ref_readlen_i)))
                
                if int(depth[idx])>0:
                    alt_maf.append('%5.4f' %(int(alt_depth[idx])/depth[idx]))
                else: 
                    alt_maf.append('0.0')
        
            linecount += runlen
        
        
        #-------------------------------------------------------
        #-----------------------Deletions------------------------
        #--------------------------------------------------------
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

            ref_q_i = []; ref_mmqs_i = []; ref_readlen_i = []; refpos_i = []
            npos = 0

            for jj in range(linecount,linecount+runlen):
                npos += 1
            
                for ib, base in enumerate(bases):
                    if refb[npos-1] == base:
                        stats = fields[jj].split('\t')[ib+5].strip('\n')
                        ref_q_i.append(float(stats.split(':')[2]))
                        ref_mmqs_i.append(float(stats.split(':')[9]))
                        ref_readlen_i.append(float(stats.split(':')[-2]))
                        refpos_i.append(float(stats.split(':')[7]))
            
            # stats for deletion reads
            depth.append(int(fields[linecount].split('\t')[3]))
            tot = len(fields[linecount].split('\t')) #total tab-seprated fields
            
            # get readcounts for all deletions at the first position
            counts = [int(fields[linecount].split('\t')[j].split(':')[1]) for j in range(10,tot)]
           
            if len(counts) >0:
                # index of the most abundant deletion
                ib = counts.index(max(counts)) + 10
                    
                delstats = fields[linecount].split('\t')[ib].strip('\n')
                alt_q = float(delstats.split(':')[2])
                altpos.append(float(delstats.split(':')[7]))
                alt_mmqs = float(delstats.split(':')[9])
                alt_readlen = float(delstats.split(':')[-2])
                dist3p.append(float(delstats.split(':')[-1]))
                alt_depth.append(int(delstats.split(':')[1]))

                diff_q.append('%6.2f' %(np.mean(ref_q_i) - alt_q))
                diff_mmqs.append('%6.2f' %(alt_mmqs - np.mean(ref_mmqs_i)))
                diffpos.append('%3.2f' %(np.mean(refpos_i) - altpos[idx]))
                diff_readlen.append('%6.2f' %(np.mean(ref_readlen_i) - alt_readlen))
                alt_readlen_.append(alt_readlen)
                diff_expec_readlen.append('%6.2f' %(float(exp_readlen[idx])-alt_readlen))
                #diff_expec_readlen.append('%3.2f' %(float(exp_readlen[idx])-np.mean(ref_readlen_i)))
                
                if int(depth[idx])>0:
                    alt_maf.append('%5.4f' %(alt_depth[idx]/depth[idx]))
                else:
                    alt_maf.append('0.0')

            else:
                ###determine the most abundant variant allele when no deletion is detected by bamrc
                allelefields = fields[linecount].strip().split('\t')[5:]
                #index of the reference base at the 1st deletion position
                refidx = [el.split(':')[0] for el in allelefields].index(refb[0])

                ref_q = float(allelefields[refidx].split(':')[2])
                ref_mmqs = float(allelefields[refidx].split(':')[9])
                ref_readlen = float(allelefields[refidx].split(':')[-2])
                ref_pos = float(allelefields[refidx].split(':')[7])

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
                    altpos.append(float(allelefields[maxalleleidx].split(':')[7]))
                    dist3p.append(float(allelefields[maxalleleidx].split(':')[-1]))

                    #alt_q_all.append(alt_q)
                    diff_q.append('%6.2f' %(ref_q - alt_q))
                    diff_mmqs.append('%6.2f' %(alt_mmqs - ref_mmqs))
                    diffpos.append('%3.2f' %(ref_pos - altpos[idx]))
                    diff_readlen.append('%6.2f' %(ref_readlen - alt_readlen))
                    diff_expec_readlen.append('%6.2f' %(float(exp_readlen[idx])-alt_readlen))
                    #diff_expec_readlen.append('%3.2f' %(float(exp_readlen[idx])-ref_readlen))

                else:
                    #no variant allele (SNP) detected at this position; set the metrics to '.'
                    altpos.append('.')
                    alt_depth.append('0')
                    diff_q.append('.')
                    diff_mmqs.append('.')
                    diff_readlen.append('.')
                    diffpos.append('.')
                    dist3p.append('.')
                    diff_expec_readlen.append('.')
                    #diff_expec_readlen.append('%3.2f' %(float(exp_readlen[idx])-np.mean(ref_readlen_i)))
                
                if int(depth[idx])>0:
                    alt_maf.append('%5.4f' %(int(alt_depth[idx])/depth[idx]))
                else: 
                    alt_maf.append('0.0')
            
            linecount += runlen
        
    #clean up temp files
    #os.remove('temp.out')
    #os.remove('temp.region')
        
    return diff_q, diff_mmqs, diff_readlen, altpos, diffpos, dist3p, alt_readlen_, diff_expec_readlen ,alt_maf, alt_depth,depth


def expected_readlen(bedfile, bedHeader, annoList):
    
    """Calculates the expected length of the reads supporting a variant from the length of amplicon(s)
    supporting it
    @param: the project bed file, 
    @param: header present in project bed? 
    @param: the annovar list for variants
    @return: the expected read length for a variant based on supporting amplicon(s)
    """

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

    
    exp_readlength = []
    of = open('expected_readlen.txt','w')
    for v in annoList:
        chrom = v.split('\t')[0]
        start = v.split('\t')[1]
        nVarAmp = 0
        amplength = []
        #no. of amplicons spanning the chromosome
        namp = len(ampStarts.get(chrom))
        for ia in range(0,namp):
            #no. of amplicons that span the variant
            if (int(start) >= int(ampStarts.get(chrom)[ia])) and (int(start) <= int(ampEnds.get(chrom)[ia])):
                nVarAmp += 1
                amplength.append(int(ampEnds.get(chrom)[ia]) - int(ampStarts.get(chrom)[ia]) + 1)
            # if amplicon start is greater than variant pos, go to next variant
            if int(start) < int(ampStarts.get(chrom)[ia]):
                break
        exp_readlength.append(np.mean(amplength))

        of.write(str(np.mean(amplength)) + '\n')
    
    of.close()
   
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


