"""
1. Apply all the filters one-by-one;
2. Write out variants after applying all filters
"""

from __future__ import division 

__Author__ =    "Rahul K. Das"
__Date__ = "May 4, 2016"
__Version__ = "2.0"
__LastModified__ = "Feb 25, 2017"


import os, sys, subprocess
import warnings
import scipy.stats as stats
from itertools import groupby
from processor1 import *
from varfile_parser import *
from processor2 import*
from qc_pcr_artifacts import *
from getSequence import *
from cal_hplen import *
import time




def process_somatic(tumorvcf, tumorbam, bedfile, bedHeader, blackfile, softdir, bamrcPath, ref_fasta, \
        varListAll, annoListAll, mnp_idx, mnp_AAChange, cosmic, \
        readq, baseq, alpha1, p_stbias, hplen_indel, hplen_snp, hplen_mnp, diffq, diffmmq, diffread1, diffread2, readpos, delpos, dist3, vardepth, varfreq, ncoshit, contextLen, run_mode, ncpu):

    """
    Apply all filters one-by-one, write different tier variant info with all QC metrics and annotations,
    
    """

    # create output files
    filt_var_anno_file = open('Somatic_Filtered_Dups.tsv', 'w') # output file to write stage-1 filter-passed variants
    
    # header for output files
    header = ['Chr','Pos','Context','Ref','Alt','Normal_AF','Normal_Alt_Depth','Normal_Tot_Depth','Tumor_AF','Tumor_Alt_Depth',\
            'Tumor_Tot_Depth','Variant_Function','Variant_Gene','Variant_Exonic_Function','AA_Change','Damaging?','Cosmic_Counts','p_Strand_Bias',\
            'HP_Length', 'Q:Ref-Alt', 'MMQS:Alt-Ref', 'ReadLength:Ref-Alt','ReadLength:Exp-Ref','Position','Pos:Ref-Alt', 'Distance_to_3p']
    filt_var_anno_file.write('\t'.join(header)+'\n')

    
    #-------------------------------------------------------------------------------- 
    # -------------------apply filters one-after-another-----------------------------
    #-------------------------------------------------------------------------------- 

    # feed the output of the previous filter to the next filter
    
    #remove_germline = 'no'
    # Remove known germline
    #if remove_germline == 'yes':
    varList, annoList, idxList, ngermline, germ_maf = removeGerm(varListAll, annoListAll)
    
    #else:
    #    ngermline = 0
    #    varList = varListAll
    #    annoList = annoListAll
    #    idxList = range(len(varListAll))
    

    # Fisher's exact test
    if run_mode == "paired":
        varList, annoList, idxList, nFishFail = FishTest(varList, annoList, idxList, alpha1)
    else:
        nFishFail = 0

    # Binomial Test: For now set the p-value to 1 to disable this
    #varList, annoList, idxList, nBinomFail = BinomTest(varList, annoList, idxList, alpha2)
    
    # rescue_cosmic
    #varList, annoList, idxList, tmp = rescue_cosmic(varListAll, annoListAll, varList, annoList, idxList, cosmic)
    varList, annoList, idxList, tmp = rescue_cosmic(varListAll, annoListAll, idxList, cosmic)

    
    #varList, annoList, idxList, nhighstbias = strbias_filter(varList, annoList, idxList, tumorvcf, p_stbias)
    #varList, annoList, idxList, nlonghp = hplen_filter(varList, annoList, idxList, tumorvcf, hplen)
   
    
    # get strand bias metric from TVC VCF for stage-2 filtering
    # this will be removed later
    if os.path.exists('Somatic_Unfiltered.vcf'):
        tempVCF = 'Somatic_Unfiltered.vcf'
    elif os.path.exists('_tmp_/Somatic_Unfiltered.vcf'):
        tempVCF = '_tmp_/Somatic_Unfiltered.vcf'
    else:
        print "The VCF file with all unfiltered variants not found"
        sys.exit()


    vcfDict = {}; stbp = []; hrun_tvc = []
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

                #if 'HRUN' in vcfDict:
                #    hrun_tvc.append(vcfDict['HRUN'])
                #else:
                #    hrun_tvc.append('.')


    #stbp_final = get_vcfInfo(tumorvcf, varList)[0] # strand bias for final vars
    #print stbp_final 
    
     
    print " Calculating variant qc metrics from bam..."
    diff_q, diff_mmqs, diff_readlen, altpos, diffpos, dist3p, alt_readlen, diff_expec_readlen, alt_maf, alt_depth, depth = \
            processor2(readq, baseq, bamrcPath, ref_fasta, tumorbam, bedfile, bedHeader, varList, annoList, ncpu) # other varscan filters

    
    # get the Damaging/Tolerant stats for these variants
    DList = damagingStats(annoList,idxList, mnp_idx)
   
    
    # get the motif sequence and homopolymer length around the variants
    with open('temp.var','w') as f:
        for el in annoList:
            f.write('\t'.join(el.split('\t')[:2]) + '\n')

    
    motifList1 = getSequence('temp.var',ref_fasta,contextLen) #for priting motifs
    
    scanLength = 20
    motifList2 = getSequence('temp.var',ref_fasta,scanLength) #for HPLen 
    os.remove('temp.var')
    

    #------------------------------------------------------------------------------
    # write the stage-1 filtered variants (non-germline) list with all the info
    #------------------------------------------------------------------------------

    for idx, el in enumerate(annoList):
        tvc_depth = varList[idx]
        if run_mode == 'paired':
            tvc_N_depth = int(tvc_depth.split('\t')[6]) + int(tvc_depth.split('\t')[7])
        else:
            tvc_N_depth = '.'

        ref = el.split('\t')[3]
        alt = el.split('\t')[4]
        temppos = el.split('\t')[1]
        
        
        ##Sequence Context##
        motifseq = motifList1[idx][:contextLen] + motifList1[idx][contextLen].lower() + motifList1[idx][-contextLen:]
        
        ##adjusted HPLEN##
        #if len(ref) > 20:
        #    scanLength = len(ref) + 20

        if len(ref) < 15:
            hrun_adj = cal_hplen(motifList2[idx],scanLength,ref,alt)
        else:
            #need to think how to handle weird long indels
            hrun_adj = '-1'

        #if len(ref) == 1 and len(alt) == 1 and ref != '-' and alt != '-': 
        #    #convert the snp base to small letter
        #    motifseq = motifList1[idx][:contextLen] + motifList1[idx][contextLen].lower() + motifList1[idx][-contextLen:]
        #else:
        #    motifseq = '.'
        
        
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
       


        #-------------------------------------------------------------------------------- 
        #-------------------write all the relevant columns for the variants--------------
        #--------------------------------------------------------------------------------

        if el.split('\t')[5] == 'exonic' and el.split('\t')[9] != 'UNKNOWN'\
            and idxList[idx] not in mnp_idx:
            #variants that have aa change info in annovar 1st run    
            if len(el.split('\t')[9].split(',')[0].split(':')[-1].split('.')) == 2:
                filt_var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                    [motifseq]+\
                    el.split('\t')[3:5]+\
                    tvc_depth.split('\t')[5:7]+ [str(tvc_N_depth)] + \
                    map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                    el.split('\t')[5:7]+\
                    ['_'.join(el.split('\t')[8].split(' '))]+\
                    [el.split('\t')[9].split(',')[0].split(':')[-1].split('.')[1]]+\
                    [DList[idx]]+[str(counts)]+[stbp[idx],hrun_adj]+\
                    map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], diff_expec_readlen[idx], altpos[idx], diffpos[idx], dist3p[idx]]))+'\n')
            #variants that don't have aa change info in both annovar's 1st & 2nd run
            else:
                filt_var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                    [motifseq]+\
                    el.split('\t')[3:5]+\
                    tvc_depth.split('\t')[5:7]+ [str(tvc_N_depth)] + \
                    map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                    el.split('\t')[5:7]+\
                    ['_'.join(el.split('\t')[8].split(' '))]+\
                    ['.']+\
                    [DList[idx]]+[str(counts)]+[stbp[idx],hrun_adj]+\
                    map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], diff_expec_readlen[idx], altpos[idx], diffpos[idx], dist3p[idx]]))+'\n')
        
        elif el.split('\t')[5] == 'exonic' and el.split('\t')[9] == 'UNKNOWN'\
            and idxList[idx] not in mnp_idx:
            filt_var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                [motifseq]+\
                el.split('\t')[3:5]+\
                tvc_depth.split('\t')[5:7]+ [str(tvc_N_depth)]+ \
                map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                el.split('\t')[5:7]+\
                ['_'.join(el.split('\t')[8].split(' '))]+\
                [el.split('\t')[9]]+\
                [DList[idx]]+[str(counts)]+[stbp[idx],hrun_adj]+\
                map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], diff_expec_readlen[idx], altpos[idx], diffpos[idx], dist3p[idx]]))+'\n')
                    
        # mnps and nfs substitutions for which annovar was run 2nd time
        elif el.split('\t')[5] == 'exonic' and idxList[idx] in mnp_idx:
            # mnps
            if len(el.split('\t')[3]) == len(el.split('\t')[4]):
                filt_var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                [motifseq]+\
                el.split('\t')[3:5]+\
                tvc_depth.split('\t')[5:7]+ [str(tvc_N_depth)] + \
                map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                el.split('\t')[5:7]+\
                ['nonsynonymous_MNV']+\
                mnp_AAChange[mnp_idx.index(idxList[idx])]+\
                [DList[idx]]+[str(counts)]+[stbp[idx],hrun_adj]+\
                map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], diff_expec_readlen[idx], altpos[idx], diffpos[idx], dist3p[idx]]))+'\n')
            else:
                #nfs indels
                filt_var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                [motifseq]+\
                el.split('\t')[3:5]+\
                tvc_depth.split('\t')[5:7]+ [str(tvc_N_depth)] + \
                map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                el.split('\t')[5:7]+\
                ['nonframeshift_substitution']+\
                mnp_AAChange[mnp_idx.index(idxList[idx])]+\
                [DList[idx]]+[str(counts)]+[stbp[idx],hrun_adj]+\
                map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], diff_expec_readlen[idx], altpos[idx], diffpos[idx], dist3p[idx]]))+'\n')
        else:
             filt_var_anno_file.write('\t'.join(el.split('\t')[:2]+\
                [motifseq]+\
                el.split('\t')[3:5]+\
                tvc_depth.split('\t')[5:7]+ [str(tvc_N_depth)] + \
                map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
                el.split('\t')[5:7]+\
                ['.','.']+\
                [DList[idx]]+[str(counts)]+[stbp[idx],hrun_adj]+\
                map(str,[diff_q[idx], diff_mmqs[idx], diff_readlen[idx], diff_expec_readlen[idx], altpos[idx], diffpos[idx], dist3p[idx]]))+'\n')
        
        
    filt_var_anno_file.close()

    #-------------Remove Duplicate Variants that may appear due to semi-duplicate entries in VCF file"
    os.system("awk '!x[$0]++' %s > %s" %('Somatic_Filtered_Dups.tsv','Somatic_Filtered_noDups.tsv'))


    #-------------Change Diff_Q, Diff_MMQS, Diff_readlen1,Diff_pos column to '.' for HOMs with AF = 1 in T
    oldtempfile = open('Somatic_Filtered_noDups.tsv','r')
    newtempfile = open('Somatic_Filtered_temp.tsv','w')
    #newtempfile = open('Somatic_Filtered_Unsorted.tsv','w')

    all_maf = [] # mafs of non-germline variants
    for il,line in enumerate(oldtempfile):
        # write header
        if il == 0:
            newtempfile.write(line)
        else:
            #Most likely HET
            if float(line.split('\t')[8]) < 0.8:
                newtempfile.write(line)
                all_maf.append(float(line.split('\t')[8]))
            else:
            # fix for HOM 
               #ref allele depth >=10; difference makes sense, write those
                if (float(line.split('\t')[10]) - float(line.split('\t')[9])) >= 10:
                    newtempfile.write(line)
                    all_maf.append(float(line.split('\t')[8]))
                else:
                    #ref allele depth too low to calculate meaningful stats; fill with '.'
                    newline = '\t'.join(line.split('\t')[:19]+['.','.','.']+[line.split('\t')[-4],line.split('\t')[-3],'.',line.split('\t')[-1]])
                    newtempfile.write(newline)
                    all_maf.append(float(line.split('\t')[8]))

    oldtempfile.close(); newtempfile.close()


    #-----------Append percentile score of variants on the basis of their MAF value
    with open('Somatic_Filtered_temp.tsv','r') as f:
        next(f)
        all_maf = [float(line.split('\t')[8]) for line in f]
    
    newtempfile = open('Somatic_Filtered.tsv','w')
    newheader = header + ['Percentile_Score']
    newtempfile.write('\t'.join(newheader) + '\n')

    oldtempfile = open('Somatic_Filtered_temp.tsv','r')
    next(oldtempfile)
    for line in oldtempfile:
        maf = float(line.split('\t')[8])
        pctile_score = '%4.1f' %stats.percentileofscore(sorted(all_maf),maf,kind='mean')
        newline = '\t'.join(line.strip('\n').split('\t') + [str(pctile_score)]) + '\n'
        newtempfile.write(newline)

    oldtempfile.close(); newtempfile.close()
                

    # sort the variants by chromosome no. because the rescued cosmic variants are appened
    #os.system('sort -t, -k1 %s > %s' %('Somatic_Filtered_Unsorted.tsv','Somatic_Filtered.tsv'))

    os.remove('Somatic_Filtered_Dups.tsv'); os.remove('Somatic_Filtered_noDups.tsv'); os.remove('Somatic_Filtered_temp.tsv')
    #os.remove('Somatic_Filtered_Unsorted.tsv')

    
    
    #-----------------PCR errors: Mispriming or amplicon bias---------------------------------------
    #NOTE! a lot of pcr error FPs are already filtered out by the read length diff filter
    #   so, if the PCR filter does not filter out any at this stage, that doesn't mean that there were no PCR artifacts

    #get the actual no. of PCR artifacts from the Somatic_Filtered.tsv file
    os.system('''awk '{if (NR>1) print $1"\t"$2"\t"$4"\t"$5}' Somatic_Filtered.tsv > temp1.region''')
    pass_fail_status_1 = qc_pcr_artifacts('temp1.region', tumorbam, bedfile, bedHeader)
    n_misprime = len([el for el in pass_fail_status_1 if el == "mispriming"])
    n_ampbias = len([el for el in pass_fail_status_1 if el == "amplicon_bias"])
    os.remove('temp1.region')
    
    #append the mispriming/amplicon bias tag
    tempList = []
    with open('Somatic_Filtered.tsv','r') as f:
        next(f)
        for line in f:
            tempList.append(line.strip().split('\t'))

    of = open('Somatic_Filtered_Detailed.tsv','w')
    newheader1 = newheader + ['PCR_error']
    of.write('\t'.join(newheader1) + '\n')
    for i in range(len(tempList)):
        of.write('\t'.join(tempList[i] + [pass_fail_status_1[i]]) + '\n')

    of.close()
    os.remove('Somatic_Filtered.tsv')
    os.system('mv Somatic_Filtered_Detailed.tsv Somatic_Filtered.tsv')



    
    #---------------------Apply varscan FP filters--------------------------------------------------
    print 'Applying bamfile-based FP filters....'
    
    # these are the subsets:
    #1. cosmic entries: synonymous SNVs that are post-annotated (*) are excluded; alt. reads >= 5
    #2. HOMs with MAF=1: hplen,readlen diff w.r.t. the panel, & position filters apply (when these values are available) 
    #3. variants for which str. bias info not present, other qc metrics are used; alt. reads >= 5
    #4. variants for which bam-readcount metrics are not present, hplen & readlen diff w.r.t the panel apply; alt. reads >= 5
    #5: variants for which all fields are present
    ##the total recomputed read depth in the T has to be >0 for all the filters
   
    #if ncoshit > 0:
    #    command = '''awk '{if (NR==1 || \
    #        ($11>0 && $17>=%s) || ($11>0 && $10>=5 && $17=="*" && $14!="synonymous_SNV") || \
    #        ($11>0 && $9==1 && $10>=%s && $19!="." && $19<=%s && $23!="." && $23<=%s && $24!="." && $24>=%s) || \
    #        ($11>0 && $10>=%s && $18=="." && $19<=%s && $20<=%s && $21<=%s && $22>=-%s && $22<=%s && $23<=%s && $24>=%s && $25>=-%s && $25<=%s && $26>=%s) || \
    #        ($11>0 && $10>=%s && $20=="." && $19<=%s && $23<=%s) || \
    #        ($11>0 && $10>=%s && $18<=%s && $19<=%s && $20<=%s && $21<=%s && $22>=-%s && $22<=%s && $23<=%s && $24>=%s && $25>=-%s && $25<=%s && $26>=%s)) print }'\
    #        Somatic_Filtered.tsv > Somatic_HQ.tsv''' %(
    #                ncoshit,
    #                vardepth,hplen,diffread2,readpos,
    #                vardepth,hplen,diffq,diffmmq,diffread1,diffread1,diffread2,readpos,delpos,delpos,dist3,
    #                vardepth,hplen,diffread2,
    #                vardepth,p_stbias,hplen,diffq,diffmmq,diffread1,diffread1,diffread2,readpos,delpos,delpos,dist3)
    #
    ##if user does not want to use cosmic based whitelisting
    #elif ncoshit == 0:
    #    command = '''awk '{if (NR==1 || \
    #        ($11>0 && $9==1 && $10>=%s && $19!="." && $19<=%s && $23!="." && $23<=%s && $24!="." && $24>=%s) || \
    #        ($11>0 && $10>=%s && $18=="." && $19<=%s && $20<=%s && $21<=%s && $22>=-%s && $22<=%s && $23<=%s && $24>=%s && $25>=-%s && $25<=%s && $26>=%s) || \
    #        ($11>0 && $10>=%s && $20=="." && $19<=%s && $23<=%s) || \
    #        ($11>0 && $10>=%s && $18<=%s && $19<=%s && $20<=%s && $21<=%s && $22>=-%s && $22<=%s && $23<=%s && $24>=%s && $25>=-%s && $25<=%s && $26>=%s)) print }'\
    #        Somatic_Filtered.tsv > Somatic_HQ.tsv''' %(
    #                vardepth,hplen,diffread2,readpos,
    #                vardepth,hplen,diffq,diffmmq,diffread1,diffread1,diffread2,readpos,delpos,delpos,dist3,
    #                vardepth,hplen,diffread2,
    #                vardepth,p_stbias,hplen,diffq,diffmmq,diffread1,diffread1,diffread2,readpos,delpos,delpos,dist3)
    #      
    #os.system(command)
    

   
    nvarDepth = 0
    nPCR_error = 0
    nwhite = 0
    n_hp_indel = 0
    n_hp_snp = 0
    n_hp_mnp = 0


    nf = open('Somatic_HQ_temp.tsv','w')
    qf = open('Somatic_QC_Failed.tsv','w')
    nf.write('\t'.join(newheader1) + '\n')
    qf.write('\t'.join(newheader1) + '\n')
    with open('Somatic_Filtered.tsv','r') as of:
        next(of)
        for line in of:
            fields = line.strip().split('\t')
            #only investigate positions with > 0 total depth
            if int(fields[10]) > 0:
                #whitelist cosmic hits
                if (ncoshit > 0 and int(fields[16]) > ncoshit):
                    nf.write('\t'.join(fields) + '\n')
                    nwhite += 1
                else:
                    #only include positions that don't have mispriming or amplicon bias
                    if fields[-1] != 'amplicon_bias' and fields[-1] != 'mispriming':
                        #variant allele depth filter
                        if int(fields[9]) >= vardepth:
                            #treat homozygous
                            if float(fields[8]) == 1:
                                if (fields[22] != '.' and float(fields[22]) <= diffread2 and \
                                    fields[23] != '.' and float(fields[23]) >= readpos):
                                        nf.write('\t'.join(fields) + '\n')
                                else:
                                    qf.write('\t'.join(fields) + '\n')
                            else:
                                #no bamrc info
                                if fields[19] == '.':
                                    if float(fields[22]) <= diffread2:
                                        nf.write('\t'.join(fields) + '\n')
                                    else:
                                        qf.write('\t'.join(fields) + '\n')
                                else:
                                    #no strand bias info
                                    if fields[17] == '.':
                                        if (abs(float(fields[19])) <= diffq and float(fields[20]) <= diffmmq and \
                                            abs(float(fields[21])) <= diffread1 and float(fields[22]) <= diffread2 and \
                                            float(fields[23]) >= readpos and abs(float(fields[24])) <= delpos and float(fields[25]) >= dist3):
                                                nf.write('\t'.join(fields) + '\n')
                                        else:
                                            qf.write('\t'.join(fields) + '\n')
                                    else:
                                        #all fields present and not HOM
                                        if (float(fields[17]) < p_stbias and abs(float(fields[19])) <= diffq and \
                                            float(fields[20]) <= diffmmq and abs(float(fields[21])) <= diffread1 and float(fields[22]) <= diffread2 and \
                                            float(fields[23]) >= readpos and abs(float(fields[24])) <= delpos and float(fields[25]) >= dist3):
                                                nf.write('\t'.join(fields) + '\n')
                                        else:
                                            qf.write('\t'.join(fields) + '\n')
                        else:
                            nvarDepth += 1
                            qf.write('\t'.join(fields) + '\n')
                    else:
                        nPCR_error += 1
                        qf.write('\t'.join(fields) + '\n')
    
    nf.close()
    
    
    #now deal with homopolymeric FPs and indel FPs with hplen filter along with the vardepth and VAF filters
    #logic: if homopolymer above a threshold or indel: then both vardepth and vaf filters need to be passed
    nf = open('Somatic_HQ.tsv','w')
    nf.write('\t'.join(newheader1) + '\n')
    with open('Somatic_HQ_temp.tsv','r') as of:
        next(of)
        for line in of:
            #whitelist cosmic
            fields = line.strip().split('\t')
            if (ncoshit > 0 and int(fields[16]) > ncoshit):
                nf.write('\t'.join(fields) + '\n')
            else:
                #snps
                if len(fields[3]) == len(fields[4]) and len(fields[3]) == 1 and fields[3] != '-' and fields[4] != '-':
                    #AAAAAaTTTTTT --> AAAAAtTTTTTT
                    if (fields[2][contextLen - 1] == fields[4] or fields[2][contextLen + 1] == fields[4]):
                        if int(fields[18]) <= hplen_snp and (int(fields[9]) >= 15 and float(fields[8]) >= 0.1):
                            nf.write('\t'.join(fields) + '\n')
                        else:
                            n_hp_snp += 1
                            qf.write('\t'.join(fields) + '\n')
                    else:
                        nf.write('\t'.join(fields) + '\n')
                #mnps
                elif len(fields[3]) == len(fields[4]) and len(fields[3]) == 2:
                    #TCTCAgAAAAA --> TCTCAaGAAAA
                    if int(fields[18]) <= hplen_mnp and (int(fields[9]) >= 15 and float(fields[8]) >= 0.1):
                        nf.write('\t'.join(fields) + '\n')
                    else:
                        n_hp_mnp += 1
                        qf.write('\t'.join(fields) + '\n')
                else:
                    #indels
                    if int(fields[18]) <= hplen_indel and (int(fields[9]) >= 15 and float(fields[8]) >= 0.1):
                        nf.write('\t'.join(fields) + '\n')
                    else:
                        n_hp_indel += 1
                        qf.write('\t'.join(fields) + '\n')

                                
    nf.close()
    qf.close()


    
    #--------------Remove variants that appear in the paired N sample: Novel germline or sequencing errors--------------
    # the filter doesn't apply to the cosmic hits with counts <=5 or those annotated by '*'
    if run_mode == 'paired':
        command = '''awk '{if (NR==1 || $7 == 0 || $17>=5) print}' Somatic_HQ.tsv > temp.tsv'''
        os.system(command)
    else:
        command = 'cp Somatic_HQ.tsv temp.tsv'
        #command = '''awk '{if (NR>1) print}' Somatic_HQ.tsv > temp.tsv'''
        os.system(command)


    #--------------Sequencing error: Remove variants that appear in the non-paired N sample-----------
    
    # the blacklist consists of non-cosmic & non-germline variants that appreared in N sample(s) [PON]
    if blackfile != '':
        with open(blackfile,'r') as bf:
            next(bf)
            blackList = ['_'.join(line.split('\t')[1:5]) for line in bf]

        newfile1 = open('Somatic_UHQ.tsv','w')
        newfile2 = open('Somatic_removed_blackList.tsv','w')
        newfile1.write('\t'.join(newheader1)+'\n')

        # filter out variants that overlap with blacklisted position
        # match contig, pos, ref, alt
        with open('temp.tsv','r') as vf:
            next(vf)
            for line in vf:
                if '_'.join(line.split('\t')[:2] + line.split('\t')[3:5]) not in blackList:
                    newfile1.write(line)
                else:
                    newfile2.write(line)

        newfile1.close(); newfile2.close()
        try:
            nfc1 = int(subprocess.check_output("grep -c '' temp.tsv", shell=True))
        except subprocess.CalledProcessError:
            nfc1 = 0
        try:
            nfc2 = int(subprocess.check_output("grep -c '' Somatic_UHQ.tsv", shell=True))
        except subprocess.CalledProcessError:
            nfc2 = 0

        nblack = nfc1 - nfc2
        os.remove('temp.tsv')
    else:
        os.rename('temp.tsv','Somatic_UHQ.tsv')
        nblack = 0

    os.system("grep -E 'Context|exonic|splicing' Somatic_UHQ.tsv > Somatic_exonic_splicing_UHQ.tsv") 
    
    #-----------------PCR errors: Mispriming or amplicon bias---------------------------------------

    #get the stats for no. of PCR FPs left after read QC and filtered out at the final stage
    #os.system('''awk '{print $1"\t"$2"\t"$4"\t"$5}' Somatic_UHQ_temp.tsv > temp.region''')
    #pass_fail_status = qc_pcr_artifacts('temp.region', tumorbam, bedfile, bedHeader)
    #os.remove('temp.region')

    ##print pass_fail_status
    #newfile = open('Somatic_UHQ.tsv','w')
    #newfile.write('\t'.join(newheader)+'\n')
    #with open('Somatic_UHQ_temp.tsv','r') as f:
    #    for ix, line in enumerate(f):
    #        if pass_fail_status[ix] != "mispriming" or pass_fail_status[ix] != "amplicon_bias":
    #            newfile.write(line)
    #newfile.close()
    #
    #try:
    #    n1 = int(subprocess.check_output("grep -c '' Somatic_UHQ_temp.tsv", shell=True))
    #except subprocess.CalledProcessError:
    #    n1 = 0
    #try:
    #    n2 = int(subprocess.check_output("grep -c '' Somatic_UHQ.tsv", shell=True))
    #except subprocess.CalledProcessError:
    #    n2 = 0
    #nPCR_error_final = n1 + 1- n2
    #
    #os.remove('Somatic_UHQ_temp.tsv')


    return ngermline, nFishFail, nwhite, nblack, n_misprime, n_ampbias, nPCR_error, nvarDepth, n_hp_snp, n_hp_indel, n_hp_mnp, all_maf
