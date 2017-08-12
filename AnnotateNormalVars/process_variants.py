"""
Analysis of all variant positions in the N samples
calculate bam-readcount metrics
write a file with for PON variants: these are non-cosmic (n<15) & non-germline variants

"""


from __future__ import division

__Author__ =    "Rahul K. Das"
__Date__ = "May 14, 2016"
__Version__ = "2.0"
__LaseModified__ = "Jan 23, 2017"


from optparse import OptionParser
import os, sys, subprocess
import datetime
import scipy.stats as stats
import statsmodels.stats.weightstats as smws
import statsmodels.stats.multitest as smm
import numpy as np
from processor1 import *
from processor2 import*
from getSequence import *



def process_variants(normalbam, varListAll, annoListAll, mnp_idx, mnp_AAChange, cosmic, bamreadcount, ref_fasta, ncpu):


    """
    @Get the list of only variants in N samples that are not cosmic (n<15) and not known germline
    @Note, in the somatic pipeline the germline variants are already filtered out and therefore, it
    makes sense to filter out the germline variants here also for the comparison later

    """

    # create output files
    var_anno_file = open('N_Vars_Dups.tsv', 'w') 
    
    # header for output files
    header = ['Chr','Pos','Ref','Alt','VAF','Var_Depth','Total_Depth','Variant_Function','Variant_Gene','Variant_Exonic_Function','AA_Change','Cosmic_Counts','p_Strand_Bias',\
            'HP_Lenth', 'Diff_Q', 'Diff_MMQS', 'Diff_Read_Length','Position', 'Distance_to_3p']
    var_anno_file.write('\t'.join(header)+'\n')

    
    # Remove known germline
    varList, annoList, idxList, ngermline, germ_maf = removeGerm(varListAll, annoListAll)

    
    #remove the cosmic entries with counts > 15 (keep this count a bit high to
    # filter out FPs that were present in paired N samples but were kept for being cosmic entry;
    # if these cosmic variants appear in multiple N samples, these are sequencing error and not
    # that the N sample was contaminated with the T

    newAnnoList = []; newidxList = []
    for idx, el in enumerate(annoList):
        if el.split('\t')[13] == '.':
            #non-cosmic
            newAnnoList.append(el)
            newidxList.append(idxList[idx])
        else:
            #cosmic with < n counts
            cosmicHits = el.split('\t')[13].split(';')[1].split('=')[1].split(',')
            counts = 0
            for item in cosmicHits:
                start = item.find( '(' )
                counts += int(item[:start])
            if counts < 15:
                newAnnoList.append(el)
                newidxList.append(idxList[idx])

    annoList = newAnnoList
    idxList = newidxList
    
    
    # get strand bias metric from TVC VCF for stage-2 filtering
    # this will be removed later
    if os.path.exists('AllVars.vcf'):
        tempVCF = 'AllVars.vcf'
    else:
        print "The VCF file with all unfiltered variants not found"
        sys.exit(1)
    
    
    vcfDict = {}; stbp = []
    with open(tempVCF,'r') as f:
        for il, line in enumerate(f):
            if il in idxList:
                for el in line.split('\t')[7].split(';'):
                    if len(el.split('=')) == 2:
                        vcfDict[el.split('=')[0]] = el.split('=')[1]
                if 'STBP' in vcfDict:
                    stbp.append(vcfDict['STBP'])
                else:
                    stbp.append('.')


    
    # get metrics by running bam-readcount on Normal bam files
    diff_q, diff_mmqs, diff_readlen, pos, dist3p, alt_maf, alt_depth, depth = processor2(annoList, bamreadcount, ref_fasta, normalbam, ncpu)

    
    
    #get HP length
    with open('temp.var','w') as f:
        for el in annoList:
            f.write('\t'.join(el.split('\t')[:2]) + '\n')

    motifList2 = getSequence('temp.var',ref_fasta,8) #for HPLen 
    os.remove('temp.var')
    
    
    
    #-----------------------------------------------------------------------
    # write all the variants with all the annotations  
    #-----------------------------------------------------------------------

    for idx, el in enumerate(annoList):
        ref = el.split('\t')[3]
        alt = el.split('\t')[4]

        #HP length
        hrun = str(max(sum(1 for i in g) for k,g in groupby(motifList2[idx])))


        #counts of cosmic entries
        if el.split('\t')[13] != '.' and el.split('\t')[13] != '*':
            cosmicHits = el.split('\t')[13].split(';')[1].split('=')[1].split(',')   
            counts = 0
            for item in cosmicHits:
                start = item.find( '(' )
                counts += int(item[:start])
        elif el.split('\t')[13] == '*':
            counts = '*'
        else:
            counts = 0
        
        #-----------------------now write all the relevant columns for the variants----------------
        if 'exonic' in el.split('\t')[5] and el.split('\t')[9] != 'UNKNOWN'\
            and idxList[idx] not in mnp_idx:
            #variants that have aa change info in annovar 1st run
            if len(el.split('\t')[9].split(',')[0].split(':')[-1].split('.')) == 2:
                var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                    el.split('\t')[3:5]+
                    map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                    el.split('\t')[5:7]+\
                    ['_'.join(el.split('\t')[8].split(' '))]+\
                    [el.split('\t')[9].split(',')[0].split(':')[-1].split('.')[1]]+\
                    [str(counts)]+[stbp[idx],hrun]+\
                    map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx]]))+'\n')
        
            #variants that don't have aa change info in both annovar's 1st & 2nd run
            else:
                var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                    el.split('\t')[3:5]+
                    map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                    el.split('\t')[5:7]+\
                    ['_'.join(el.split('\t')[8].split(' '))]+\
                    ['.']+\
                    [str(counts)]+[stbp[idx],hrun]+\
                    map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx]]))+'\n')
        
        elif 'exonic' in el.split('\t')[5] and el.split('\t')[9] == 'UNKNOWN'\
            and idxList[idx] not in mnp_idx:
            var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                el.split('\t')[3:5]+
                map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                el.split('\t')[5:7]+\
                ['_'.join(el.split('\t')[8].split(' '))]+\
                [el.split('\t')[9]]+\
                [str(counts)]+[stbp[idx],hrun]+\
                map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx]]))+'\n')
                    
        # mnps and nfs substitutions for which annovar was run 2nd time
        elif 'exonic' in el.split('\t')[5] and idxList[idx] in mnp_idx:
            # mnps
            if len(el.split('\t')[3]) == len(el.split('\t')[4]):
                var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                el.split('\t')[3:5]+
                map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                el.split('\t')[5:7]+\
                ['nonsynonymous_MNV']+\
                mnp_AAChange[mnp_idx.index(idxList[idx])]+\
                [str(counts)]+[stbp[idx],hrun]+\
                map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx]]))+'\n')
            else:
                #nfs indels
                var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                el.split('\t')[3:5]+
                map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                el.split('\t')[5:7]+\
                ['nonframeshift_substitution']+\
                mnp_AAChange[mnp_idx.index(idxList[idx])]+\
                [str(counts)]+[stbp[idx],hrun]+\
                map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx]]))+'\n')
        else:
             var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                el.split('\t')[3:5]+
                map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                el.split('\t')[5:7]+\
                ['.','.']+\
                [str(counts)]+[stbp[idx],hrun]+\
                map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx]]))+'\n')
        
        
    var_anno_file.close()

    
    #-------------Remove Duplicate Variants that may appear due to semi-duplicate entries in VCF file"
    os.system("awk '!x[$0]++' %s > %s" %('N_Vars_Dups.tsv','N_Vars_noDups.tsv'))

    
    #-------------Change Diff_Q, Diff_MMQS, Diff_Read_Length columns to '.' for HOMs with AF = 1 in T
    oldtempfile = open('N_Vars_noDups.tsv','r')
    newtempfile = open('N_Vars.tsv','w')

    for il,line in enumerate(oldtempfile):
        # write header
        if il == 0:
            newtempfile.write(line)
        elif il > 0:
        # fix for HOM with MF=1
            if float(line.split('\t')[4]) == 1:
                newline = '\t'.join(line.split('\t')[:14]+['.','.','.']+line.split('\t')[-2:])
                newtempfile.write(newline)
            else:
                newtempfile.write(line)

    oldtempfile.close(); newtempfile.close()
   
    os.remove('N_Vars_Dups.tsv');os.remove('N_Vars_noDups.tsv')
