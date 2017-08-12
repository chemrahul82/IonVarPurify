"""
Parses the vcf file to extract the positions with alternative AF > 0
& Makes the vcf file of this subset of variants

"""

from __future__ import division

__Author__ =    "Rahul K. Das"
__Date__ = "May 4, 2016"
__Version__ = "1.0"
__LastModified__ = "Sept 19, 2016"


from optparse import OptionParser
import os, sys, time
import scipy.stats as stats
import statsmodels.stats.weightstats as smws
import statsmodels.stats.multitest as smm
import numpy as np


def extract_variants(vcf):

    """ Extract the variants with AF>0 & make a vcf

    Input: VCF with all positions

    Returns: list with all variants AF>0

    """

    varListAll = []
    outvcf = open('AllVars.vcf', 'w')
    
    with open(vcf) as f:
        for line in f:
            if len(line.split('\t')) > 8:
                vcffields = line.split('\t')[7].split(';')
                
                if len(vcffields[0].split('=')) > 1:
                    
                    if vcffields[0].split('=')[0] == 'AF':
                        maf = vcffields[0].split('=')[1]
                    
                    else:
                        #if maf not present calculate maf=AO/DP
                        maf = float(vcffields[0].split('=')[1])/float(vcffields[1].split('=')[1])
                    
                    if len(vcffields[0].split('=')[1].split(',')) == 1 and float(maf)>0: 
                        outvcf.write(line)
                        varListAll.append(line)
                    #elif len(vcffields[0].split('=')[1].split(',')) > 1:
                        #multiple variants in the same locus: NEED TO IMPLEMENT THIS
                    #    outvcf.write(line)
                    #    varListAll.append(line)
    
    outvcf.close()

    return varListAll

