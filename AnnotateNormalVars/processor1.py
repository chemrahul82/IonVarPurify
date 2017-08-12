"""
Performs statistical tests and variants's sequencing metrics 
related filtering on the potential somatic variants based on hard AF cutoffs
1. Fisher's exact test (on pair where no alternative reads are present in N)
2. Binomial test (It is helpful if we are using a nonzero WT/HET AF cutoff for T samples 
    and want to ask the significance level of the observed AF in T;
    As we are including all non-zero AF T-variants as potential somatic, we can skip this test;
    The bam-file centric qc metrics complements this filter)
3. Remove known germline variants
4. Filter out non-germline variants with >0 MAF in N (Not calling the function from here)
5. Strand bias filter (Not calling this function from here)
6. Homopolymer length filter (Not calling this function from here)
7. Rescue non-germline cosmic entries

"""
from __future__ import division

__Author__ =    "Rahul K. Das"
__Date__ = "May 4, 2016"
__Version_ = "1.1"
__LastModified_ = "Jan 23, 2017"


from optparse import OptionParser
import os, sys
import scipy.stats as stats
import statsmodels.stats.weightstats as smws
import statsmodels.stats.multitest as smm
import numpy as np
#from varfile_parser import *


def annovar_annotate(software_dir, annovar_dir, var_vcf, ncpu):

    """ Run annovar on potential variants to get annotations
    
    Input: 1. Path of the software
           2. Path of annovar directory
           3. list with variants
           4. no. of CPU to run annovar
    
    Output: Generates the annotation file in the running directory
    returns: four lists: varlist with annotations, indices of mnps &
        nfs substitutions, AA changes for these variants,
        indices of variants that are cosmic entries but not 
        dbsnp, 1000g entries 
    
    """

    oldAnnoFile = '_tmp_/annovar_variants_annotation_table.txt'
    oldAnnoMnpLog = '_tmp_/annovar_variants_table_mnp.log'

    
    if not os.path.exists(oldAnnoFile):
        command = '%s/run_annovar.sh %s %s %s >/dev/null 2>&1' %(software_dir, annovar_dir, var_vcf, ncpu)
        os.system(command)
        AnnoFile = 'annovar_variants_annotation_table.txt'
        AnnoMnpLog = 'annovar_variants_table_mnp.log'

    else:
        print " annovar has already been run, delete the output file to rerun"
        AnnoFile = oldAnnoFile
        AnnoMnpLog = oldAnnoMnpLog



    with open(AnnoFile, 'r') as annofile:
        next(annofile)
        annolist = [line.strip('\n') for line in annofile]
    
   #get the indices of variants with mnp and their annotation
    mnp_idx = []
    mnp_AAChange = []

    with open(AnnoMnpLog, 'r') as mnpfile:
        for line in mnpfile:
            if 'position' in line:
                attributes = line.split(' ')
                lineno = int(attributes[0].replace('>line',''))
                
                refb = annolist[lineno-1].split('\t')[3]
                altb = annolist[lineno-1].split('\t')[4]
                
                #mnp that results in single AA change
                if len(refb) == len(altb) and refb != '-' and altb != '-' \
                    and 'exonic' in annolist[lineno-1].split('\t')[5]:
                    mnp_idx.append(lineno-1)
                    mnp_AAChange.append([i+j+k for i,j,k in zip(attributes[9],attributes[6].split('-')[0:],attributes[11])])
                
                # nfs substitution
                elif len(refb) != len(altb) and \
                    'nonframeshift substitution' in annolist[lineno-1].split('\t')[8]: 
                    mnp_idx.append(lineno-1)
                    mnp_AAChange.append([i+j+k for i,j,k in zip(attributes[9],attributes[6].split('-')[0:],['nfs'])])

                # fs substitution
                elif len(refb) != len(altb) and \
                    annolist[lineno-1].split('\t')[8] == 'frameshift substitution': 
                    mnp_idx.append(lineno-1)
                    mnp_AAChange.append([i+j+k for i,j,k in zip(attributes[9],attributes[6].split('-')[0:],['fs'])])


                
    # get cosmic hits that are not known germline
    cosmic = [idx for idx, var in enumerate(annolist) if\
            var.split('\t')[10].strip(' ') == '.' and \
            var.split('\t')[11].strip(' ') == '.' and \
            var.split('\t')[13].strip(' ') != '.']
    
    return annolist, mnp_idx, mnp_AAChange, cosmic


def damagingStats(annolist,idxlist,mnp_idx):
    
    """ get Damaging/Tolerant predictions from LJB database
    
    Input: list of variants with annovar annotations, index list of these variants,
        index list of mnps/nfs
    
    Returns: list with 'D'/'T'/'NA'/'.' annotations for input variants

    """
    
    DList = []
    
    for idx, el in enumerate(annolist):
        if el.split('\t')[5] == 'exonic' and 'nonsynonymous' in el.split('\t')[8]:
            ljblist = el.split('\t')[20:37:2]
            totalD = ljblist.count('D')+ljblist.count('A')+ljblist.count('H')+ljblist.count('M')
            Dscore = totalD/9
            if totalD >= 3:
                DList.append('D'+':'+'%3.2f'%(Dscore))
            else:
                DList.append('T'+':'+'%3.2f'%(Dscore))
        elif len(el.split('\t')[3]) == len(el.split('\t')[4]) and \
                len(el.split('\t')[3]) == 2 and idxlist[idx] in mnp_idx:
            DList.append('NA')
        else:
            DList.append('.')

    return DList


def fdr_test(pval_raw, alpha):

    """ Benjamini/Hochberg multiple testing correction (FDR)
    
    Input: list with raw p-values, alpha
    Returns: Boolean list:True(accept),False(reject)

    """
   
    rej_acc, pval_corr = smm.multipletests(pval_raw, alpha = alpha, method = 'fdr_bh', returnsorted = False)[:2]
    
    return rej_acc, pval_corr


def removeGerm(varList, annoList):
    
    """ Filter 1: Remove known germline variants
    
    Input: list of variants with information from matched variant file,
           list of variants with annovar annotations
    
    Returns: The above lists with germline removed,
             A list with indices of filtered variants that map with original list
             Total known germline that were removed
             MAFs of the germline variants
        
    """
   
    ngermline = 0
    germ_maf = []
    PassVarList = []
    PassAnnoList = []
    idxPassList = []

    #save 1000g_hit, dbsnp_hit, cosmic_hit info in a list
    _1000g = [var.split('\t')[10].strip(' ') for var in annoList]
    dbsnp = [var.split('\t')[11].strip(' ') for var in annoList]
    cosmic = [var.split('\t')[13] for var in annoList]
    
    for idx in range(len(annoList)):
        #known germline
        if (_1000g[idx] != '.' or dbsnp[idx] != '.'):
            ngermline += 1
            germ_maf.append(float(annoList[idx].split('\t')[54].split(';')[0].split('=')[1]))
        else:
            idxPassList.append(idx)
            PassVarList.append(varList[idx])
            PassAnnoList.append(annoList[idx])

    return PassVarList, PassAnnoList, idxPassList, ngermline, germ_maf
 

def FishTest(varList, annoList, idxList, alpha):
   
    """ Filter-2: Fisher's two-sided exact test 
    
    Input: list of variants with information from matched variant file,
           list of variants with annovar annotations
           list with indices of variants that map them with the original list
           significance level (alpha)
    
    Returns: The above lists with variants that fail the test removed,
             Total variants that were removed
     
    """       
    
    p_fish_raw = []
    PassVarList = []
    PassAnnoList = []
    idxPassList = []
    nFishTestFail = 0

    Nnzeroidx = []
    for idx in range(len(varList)):
        alt_N = int(varList[idx].split('\t')[6])
        ref_N = int(varList[idx].split('\t')[7])
        alt_T = int(varList[idx].split('\t')[10])
        ref_T = int(varList[idx].split('\t')[11])

        #test only when alt. reads are present in N because:
        #applying test when no alt reads in N and very low MAF in T can 
        #lead to false negative for low quality T samples
        if alt_N > 0:
            Nnzeroidx.append(idx)
            p_val = stats.fisher_exact([[ref_N, ref_T], [alt_N, alt_T]], alternative = 'greater')[1]
            p_fish_raw.append(p_val)


    #Apply Benjamini-Hochberg multiple testing correction
    if len(p_fish_raw) > 0:
        rej_acc_fish, p_fish_corr = fdr_test(p_fish_raw, alpha)
    
    #if all variants pass Fisher's test, set rej_acc_status to True for all
    else:
        rej_acc_fish = [True]*len(varList)

    for idx in range(len(varList)):
        if idx in Nnzeroidx:
            if rej_acc_fish[Nnzeroidx.index(idx)] == True:
                idxPassList.append(idxList[idx])
                PassVarList.append(varList[idx])
                PassAnnoList.append(annoList[idx])
                #print "Pass"
            else:
                nFishTestFail += 1
                #print "Fail"
        else:
            idxPassList.append(idxList[idx])
            PassVarList.append(varList[idx])
            PassAnnoList.append(annoList[idx])
            #print "Not applied"



    return PassVarList, PassAnnoList, idxPassList, nFishTestFail


def BinomTest(varList, annoList, idxList, alpha):
    
    """ Filter 3: BINOMIAL TEST
    
    Input: list of variants with information from matched variant file,
           list of variants with annovar annotations
           list with indices of variants that map them with the original list
           significance level (alpha)
    
    Returns: The above lists with variants that fail the test removed,
             Total variants that were removed
    
    """

    p_bin1_raw = []
    idx1 = []
    nBinomFail = 0
    ncount = 0
    PassVarList = []
    PassAnnoList = []
    idxPassList = []

    for idx in range(len(varList)):
        alt_N = int(varList[idx].split('\t')[6])
        ref_N = int(varList[idx].split('\t')[7])
        alt_T = int(varList[idx].split('\t')[10])
        ref_T = int(varList[idx].split('\t')[11])
        maf_T = float(varList[idx].split('\t')[9])
    
    # FILTER OUT variants for which MAF in N is 0,
    # but for T is close to the WT/HET cutoff for somatic
    # These are most likely WT-WT
    # Peform bionomial/ z-test on T-reads
    # We are dealing with categorical binary data, where one-sample t-test can't be used
       
        if alt_N == 0:
            idx1.append(idx)
            p_bin1_raw.append(stats.binom_test(alt_T, alt_T+ref_T, 0.04))
                
    #Apply Benjamini-Hochberg multiple testing correction to p-values from Binomial test-1
        if len(p_bin1_raw) > 0:
            rej_acc_bin1, p_bin1_corr = fdr_test(p_bin1_raw, alpha)
        else:
            rej_acc_bin1 = []

    for idx in range(len(varList)):
        # keep variants which bionomial test was not applied on
        if idx not in idx1:
            idxPassList.append(idxList[idx])
            PassVarList.append(varList[idx])
            PassAnnoList.append(annoList[idx])
        else:
            # keep variants that pass binomial test
            ncount += 1
            if rej_acc_bin1[ncount-1] == True:
                idxPassList.append(idxList[idx])
                PassVarList.append(varList[idx])
                PassAnnoList.append(annoList[idx])
            else:
                nBinomFail += 1
    
    return PassVarList, PassAnnoList, idxPassList, nBinomFail  


def removeNoise(varList, annoList, idxList):
    
    """ Filter 4: Remove variants with N_MAF > 0: NOVEL GERMLINE OR SEQUENCING ERRORS
    
    Input: list of variants with information from matched variant file,
           list of variants with annovar annotations
           list with indices of variants that map them with the original list
    
    Returns: The above lists with variants that fail the test removed,
             Total variants that were removed
    
    """

    nnzero_N_MAF = 0
    PassVarList = []
    PassAnnoList = []
    idxPassList = []

    for idx in range(len(varList)):
        alt_N = int(varList[idx].split('\t')[6])
        # filter out
        if alt_N > 0:
            nnzero_N_MAF += 1
        else:
            idxPassList.append(idxList[idx])
            PassVarList.append(varList[idx])
            PassAnnoList.append(annoList[idx])

    return PassVarList, PassAnnoList, idxPassList, nnzero_N_MAF


def strbias_filter(varList, annoList, idxList, tumorvcf, p_stbias):
    
    """ Function for filtering based on strand bias
    
    Input: list of variants with information from matched variant file,
           list of variants with annovar annotations
           list with indices of variants that map them with the original list
           tumor vcf file
           p-value threshold for filtering
    
    Returns: The above lists with variants that fail the test removed,
             Total variants that were removed
    
    """

    PassVarList = []
    PassAnnoList = []
    idxPassList = []
    nhighstbias = 0

    stbp = get_vcfInfo(tumorvcf, varList)[0]
    for idx in range(len(varList)):
        # apply filter
        # for some variants, stbp is not reported, take care of those
        if len(stbp[idx]) != 0 and float(stbp[idx][0]) < p_stbias:
            nhighstbias += 1
        else:
            idxPassList.append(idxList[idx])
            PassVarList.append(varList[idx])
            PassAnnoList.append(annoList[idx])

    return PassVarList, PassAnnoList, idxPassList, nhighstbias


def hplen_filter(varList, annoList, idxList, tumorvcf, hplen):
    
    """Function for filtering based on homopolymer length

    Input: list of variants with information from matched variant file,
           list of variants with annovar annotations
           list with indices of variants that map them with the original list
           tumor vcf file
           p-value threshold for filtering
    
    Returns: The above lists with variants that fail the test removed,
             Total variants that were removed
    
    """

    PassVarList = []
    PassAnnoList = []
    idxPassList = []
    nlonghp = 0

    hrun = get_vcfInfo(tumorvcf, varList)[1]
    for idx in range(len(varList)):
        # apply filter
        # for some variants, hrun is not reported, take care of those
        if len(hrun[idx]) != 0 and int(hrun[idx][0]) > hplen:
            nlonghp += 1
        else:
            idxPassList.append(idxList[idx])
            PassVarList.append(varList[idx])
            PassAnnoList.append(annoList[idx])
        
    return PassVarList, PassAnnoList, idxPassList, nlonghp


def rescue_cosmic(varListAll, annoListAll, varList, annoList, idxList, cosmic):
    
    """Function for rescuing non-germline cosmic entries that were filtered out

    Input: list of all potential variants with information from matched variant file
           list of all potential variants with annovar annotations
           list of variants with information from matched variant file,
           list of variants with annovar annotations
           list with indices of variants that map them with the original list
           list with cosmic annotations
    
    Returns: The three running lists with known cosmic entries (if any & was removed) added
             Total variants that were whitlisted
    
    """

    fc = open('/rawdata/software/annovar_Feb2016/humandb/hg19_cosmic70.txt', 'r')
    
    ann_idx = {}; #dictionary to map chrNo, bin, and start and end of bins (in bits)
    chrBinDict = {} #dictionary to map chr no. and bin indices
    for i in range(22):
        chrBinDict[str(i+1)] = []

    chrBinDict['X'] = []; chrBinDict['Y'] = []

    #cosmic index file
    with open('/rawdata/software/annovar_Feb2016/humandb/hg19_cosmic70.txt.idx', 'r') as fi:
        next(fi)
        for line in fi:
            chrom = line.split('\t')[0]
            bin_idx = line.split('\t')[1]
            #chrNo & bin mapping
            chrBinDict[chrom].append(int(bin_idx))

            chr_bin = chrom + '_' + bin_idx
            file_pos = line.strip('\n').split('\t')[2:4]
            #chr+bin and bin range mapping
            ann_idx[chr_bin] = file_pos

    #with open('/rawdata/software/annovar_Feb2016/humandb/hg19_cosmic70.txt', 'r') as f:
    #    cosmicdat = ['chr'+line.split('\t')[0]+'_'+line.split('\t')[1]+'_'+line.split('\t')[2] for line in f]
   
    nvar0 = len(varList)
    
    for cosidx in cosmic:
        if cosidx not in idxList:
            varList.append(varListAll[cosidx])
            annoList.append(annoListAll[cosidx])
            idxList.append(cosidx)
    
    print "Adding cosmic annotations to exonic variants that are potentially missed by Annovar...."
    for idx, el in enumerate(annoListAll):

        # only look fot exonic variants that are in the passed list, don't already have cosmic annotation
        # from annovar runs and are not known germline vars and have N_AF = 0
        if idx in idxList and idx not in cosmic and el.split('\t')[5] == 'exonic' \
            and el.split('\t')[10] == '.' and el.split('\t')[11] == '.'\
            and float(varListAll[idx].split('\t')[5]) == 0:
            
            chrNo = el.split('\t')[0].replace('chr','')
            varPosStart = int(el.split('\t')[1])
            varPosEnd = int(el.split('\t')[2])

            #map variant position to the bin in the cosmic index file
            bin_idx = int(varPosStart/1000)*1000
            #print chrNo,el.split('\t')[1],bin_idx

            search_key = chrNo + '_' + str(bin_idx)
            binRange = [value for key,value in ann_idx.items() if key == search_key]
            if len(binRange) > 0:
                bitPosStart = binRange[0][0]
                bitPosEnd = binRange[0][1]

                bitPos = int(bitPosStart)
            
                while bitPos < int(bitPosEnd):
                    fc.seek(int(bitPos))
                    cosStart = int(fc.readline().split('\t')[1])
                    cosEnd = int(fc.readline().split('\t')[2])
                    #print fc.readline()
                
                    if varPosStart < cosStart:
                        break
                
                    if varPosStart >= cosStart and varPosEnd <= cosEnd:
                        # update the cosmic field with '*'
                        annoNew = '\t'.join(annoListAll[idx].split('\t')[:13]+['*']+annoListAll[idx].split('\t')[14:])
                        annoList[idxList.index(idx)] = annoNew
                        idxList.append(idx)
                        break

                    bitPos = fc.tell()


    nvar1 = len(varList)
    #no. of variant whitelisted & rescued
    nwhite = nvar1 - nvar0

    return varList, annoList, idxList, nwhite


