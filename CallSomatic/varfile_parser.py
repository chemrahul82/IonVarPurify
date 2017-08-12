"""
@ intersects the N and T VCF files to construct the matched variant list
@ parses the matched_variants file to extract the potential somatic variants based on our hard AF cutoffs;
@ Makes the vcf file of this subset of variants

"""

__Author__ =    "Rahul K. Das"
__Date__ = "May 4, 2016"
__Version__ = "2.0"
__LastModified__ = "Nov 10, 2016"


from optparse import OptionParser
import os, sys, time
import scipy.stats as stats
import statsmodels.stats.weightstats as smws
import statsmodels.stats.multitest as smm
import numpy as np
from vcf2tsv import *
from collections import OrderedDict


def get_matched_variants(normalVCF,tumorVCF,run_mode, varcall_mode):
    
    """
    Construct the matched normal-tumor variant file
    @param: normalVCF: VCF file for Normal
    @param: tumorVCF: VCF file for Tumor
    @param: run_mode: paired or unpaired
    @param: varcall_mode: hotspot or not-hotspot
    @return: the matched variant file
    """
   
    #------------paired mode-------------------
    if run_mode == "paired":
        vcf2tsv(normalVCF, 'Normal_Variants.tsv')
        vcf2tsv(tumorVCF, 'Tumor_Variants.tsv')

        with open('Normal_Variants.tsv', 'r') as f:
            next(f)
            nVars = [line for line in f]

        with open('Tumor_Variants.tsv', 'r') as f:
            next(f)
            tVars = [line for line in f]

        nDict = OrderedDict(); tDict = OrderedDict()
        for el in nVars:
            nDict['_'.join(el.strip().split('\t')[:2])] = el.strip().split('\t')

        for el in tVars:
            tDict['_'.join(el.strip().split('\t')[:2])] = el.strip().split('\t')

        #For hotspot mode: write only those positions which are in both T & N VCF:
        #For non-hotspot mode: if a position in T VCF is not in N VCF, annotate that as WT in N;
        #   in non-hotspot mode a missing position in N VCF may also mean a 'NoCall' at that position: this may increase FPs

        outfile = open('matched_variants.tsv', 'w')
        header = '\t'.join(['Chr','Pos','Ref','Alt','GT_N','VAF_N','Alt_Depth_N','Ref_Depth_N','GT_T','VAF_T','Alt_Depth_T','Ref_Depth_T']) + '\n'
        outfile.write(header)
        for key in tDict:
            if key in nDict:
                outfile.write('\t'.join(nDict[key] + tDict[key][4:]) + '\n')
            else:
                if varcall_mode == 'not_hotspot':
                    outfile.write('\t'.join(tDict[key][:4] + ['WT','0.0','0','0'] + tDict[key][4:]) + '\n')

        outfile.close()


    #-------------unpaired mode-----------------
    elif run_mode == "unpaired":
        vcf2tsv(tumorVCF, 'Tumor_Variants.tsv')

        outfile = open('matched_variants.tsv', 'w')
        header = '\t'.join(['Chr','Pos','Ref','Alt','GT_N','VAF_N','Alt_Depth_N','Ref_Depth_N','GT_T','VAF_T','Alt_Depth_T','Ref_Depth_T']) + '\n'
        outfile.write(header)

        with open('Tumor_Variants.tsv', 'r') as f:
            next(f) 
            for line in f:
                outfile.write('\t'.join(line.strip().split('\t')[:4] +  ['.','.','.','.'] + line.strip().split('\t')[4:]) + '\n')

        outfile.close()


def extract_raw_somatic(normalvcf, tumorvcf, run_mode, varcall_mode):

    """ Extract the variants that have: 0<=N_AF<0.2 & T_AF>0

    NOTE: The matched_variants file was generated by comparing hotspot VCFs of N & T runs;

    @param: normal VCF
    @param: tumor VCF
    @param: run mode: paired/unpaired
    @param: varcall_mode: hotspot/not_hotspot
    @return: list with the info from matched_variants file for the above variants

    """

    get_matched_variants(normalvcf,tumorvcf,run_mode,varcall_mode)
    matchedVarFile = 'matched_variants.tsv'

    #output file to write all potential somatic variants
    outfile = 'Somatic_Unfiltered.txt'
    outvcf = open('Somatic_Unfiltered.vcf', 'w')
    
    #extract the variants that have 0<=N_AF<0.2 & T_AF>=0
    if run_mode == 'paired':
        command = "awk '{if ($6 >= 0 && $6 < 0.2 && $10>0) print}' %s > %s" %(matchedVarFile, outfile)
    else:
        command = "awk '{if (NR>1 && $10>0) print}' %s > %s" %(matchedVarFile, outfile)
    os.system(command)
    
    #store these potential somatic variants in a list
    with open(outfile, 'r') as f:
        somatic_all = [line.strip('\n') for line in f]
    
    # position of these raw somatic
    somaticpos = ['_'.join(el.split('\t')[:2]) for el in somatic_all]
    
    # store lines of all variants in tumor vcf in a list
    with open(tumorvcf, 'r') as f:
        allvars = [line for line in f]
    
    # position of all variants in tumor vcf
    allpos = ['_'.join(el.split('\t')[:2]) for el in allvars]
    
    # make the vcf for these potential somatic vars
    index_all = dict((pos,idx) for idx, pos in enumerate(allpos))
    [outvcf.write(allvars[index_all[pos]]) for pos in somaticpos]
    
    outvcf.close()
    
    return somatic_all


#def get_vcfInfo(tumor_vcf, varList):
#    
#    """ Get info from VCF file for a list of variants
#    
#    Input: tumor VCF
#    Returns: list with relevant variant info extracted from T-VCF
#
#    """
#
#    with open(tumor_vcf, 'r') as tvcf:
#        vcf_info_all = [line for line in tvcf if '#' not in line]
#   
#    stbp = []
#    hrun = []
#    for var1 in varList:
#        for idx, var2 in enumerate(vcf_info_all):
#            if var1.split('\t')[0] == var2.split('\t')[0]\
#                and var1.split('\t')[1] == var2.split('\t')[1]:
#                fields = var2.split('\t')[7].split(';')
#                stbp.append([i.split('=')[1] for i in fields if 'STBP' in i])
#                hrun.append([i.split('=')[1] for i in fields if 'HRUN' in i])
#                break
#    
#    return stbp, hrun
