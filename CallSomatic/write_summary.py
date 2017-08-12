"""
Compute the summary statistics & write both in text & json formats

"""

from __future__ import division 

__Author__ =    "Rahul K. Das"
__Date__ = "May 4, 2016"
__Version__ = "2.0"
__LaseModified__ = "Feb 25, 2017"


from optparse import OptionParser
import os, sys, subprocess
import warnings
import datetime
import numpy as np
from collections import OrderedDict
import json

def write_summary(varListAll, ngermline, nFishFail, nwhite1, nwhite, nblack, nDup, n_misprime, n_ampbias, nPCR_error_final, nvarDepth, nzeroVarDepth, n_hp_snp, n_hp_indel, n_hp_mnp, all_maf, normalvcf, tumorvcf, tumorbam, reffasta, bedfile, blackfile, softdir,\
        readq, baseq, alpha1, pstbias, hplen_indel, hplen_snp, hplen_mnp, diffq, diffmmq, diffread1, diffread2, readpos, delpos, dist3, vardepth, varfreq, ncoshit, contextLen, run_mode, sample_type):

    
    #--------------------------------------------------------
    #-----------------Compute Statistics---------------------
    #--------------------------------------------------------

    # no. of variants at different stages of filtration
    nvar1 = subprocess.check_output("grep -c '' Somatic_Germline_Filtered.tsv", shell=True)
    nvar2 = subprocess.check_output("grep -c '' Somatic_HQ.tsv", shell=True)
    nvar3 = subprocess.check_output("grep -c '' Somatic_UHQ.tsv", shell=True)
   
    #get summary stat from final Ultra-High Quality filtered variants list
    with open('Somatic_UHQ.tsv', 'r') as f:
        next(f)
        finalAnnoList = [line.strip('\n') for line in f]

    nTotal = len(varListAll) # no. of all potential somatic vars
    nFinal = len(finalAnnoList) # no. of UHQ vars
    
    #UHQ maf stats
    uhq_maf = [float(el.split('\t')[8]) for el in finalAnnoList]
    if len(uhq_maf) > 0:
        uhq_maf_mean = '%3.2f' %(np.mean(uhq_maf))
        uhq_maf_median = '%3.2f' %(np.median(uhq_maf))
        uhq_maf_std = '%3.2f' %(np.std(uhq_maf))
    else:
        uhq_maf_mean = 0; uhq_maf_median = 0; uhq_maf_std = 0

    #all somatic maf stats
    #if len(all_maf) > 0:
    #    all_maf_mean = np.mean(all_maf)
    #    all_maf_median = np.median(all_maf)
    #    all_maf_std = np.std(all_maf)
    #else:
    #    all_maf_mean = 0; all_maf_median = 0; all_maf_std = 0


    ncos = 0; ncos25 = 0; ncos100 = 0
    ndamage = 0; nnadamage = 0; ntol = 0;

    # counts for six types of substitutions only for the exonic
    n_CA = 0; n_CG = 0; n_CT = 0; n_TA = 0; n_TC = 0; n_TG = 0
    nTv = 0; nTs = 0 

    for idx, el in enumerate(finalAnnoList):
        if 'D' in el.split('\t')[15]:
            ndamage += 1
        if 'T' in el.split('\t')[15]:
            ntol += 1
        if 'NA' in el.split('\t')[15]:
            nnadamage += 1
        if el.split('\t')[16] != '*' and int(el.split('\t')[16]) > 0:
            ncos += 1
        if el.split('\t')[16] != '*' and int(el.split('\t')[16]) >25:
            ncos25 += 1
        if el.split('\t')[16] != '*' and int(el.split('\t')[16]) >100:
            ncos100 += 1
        
        
    # counts for exonic, intronic, nonsyn etc..
    try:
        nexonic = subprocess.check_output("grep -c 'exonic' Somatic_UHQ.tsv", shell=True)
    except subprocess.CalledProcessError:
        nexonic = 0
    try:
        nintronic = subprocess.check_output("grep -c 'intronic' Somatic_UHQ.tsv", shell=True)
    except subprocess.CalledProcessError:
        nintronic = 0
    try:
        nsplice = subprocess.check_output("grep -c 'splicing' Somatic_UHQ.tsv", shell=True)
    except subprocess.CalledProcessError:
        nsplice = 0
    try:
        nsyn = subprocess.check_output('''awk '{if ($14=="synonymous_SNV") print}' Somatic_UHQ.tsv | wc -l''', shell=True)
    except subprocess.CalledProcessError:
        nsyn = 0
    try:
        nmiss = subprocess.check_output("grep -c 'nonsynonymous' Somatic_UHQ.tsv", shell=True)
    except subprocess.CalledProcessError:
        nmiss = 0
    try:
        nfsi = subprocess.check_output("grep -c '\<frameshift_insertion\>' Somatic_UHQ.tsv", shell=True)
    except subprocess.CalledProcessError:
        nfsi = 0
    try:
        nfsd = subprocess.check_output("grep -c '\<frameshift_deletion\>' Somatic_UHQ.tsv", shell=True)
    except subprocess.CalledProcessError:
        nfsd = 0
    try:
        nfss = subprocess.check_output("grep -c '\<frameshift_substitution\>' Somatic_UHQ.tsv", shell=True)
    except subprocess.CalledProcessError:
        nfss = 0
    try:    
        nnfsi = subprocess.check_output("grep -c '\<nonframeshift_insertion\>' Somatic_UHQ.tsv", shell=True) 
    except subprocess.CalledProcessError:
        nnfsi = 0
    try:
        nnfsd = subprocess.check_output("grep -c '\<nonframeshift_deletion\>' Somatic_UHQ.tsv", shell=True) 
    except subprocess.CalledProcessError:
        nnfsd = 0
    try:
        nnfss = subprocess.check_output("grep -c '\<nonframeshift_substitution\>' Somatic_UHQ.tsv", shell=True)
    except subprocess.CalledProcessError:
        nnfss = 0
    try:
        nstop = subprocess.check_output("grep -E -c 'stopgain|stoploss' Somatic_UHQ.tsv", shell=True) 
    except subprocess.CalledProcessError:
        nstop = 0
    
    #six types of SNPs
    try:
        nCA = subprocess.check_output('''awk '{if (($12=="exonic") && (($4=="C" && $5=="A") || ($4=="G" && $5=="T"))) print}' Somatic_UHQ.tsv | wc -l''', shell=True)
    except subprocess.CalledProcessError:
        nCA = 0
    try:
        nCG = subprocess.check_output('''awk '{if (($12=="exonic") && (($4=="C" && $5=="G") || ($4=="G" && $5=="C"))) print}' Somatic_UHQ.tsv | wc -l''', shell=True)
    except subprocess.CalledProcessError:
        nCG = 0
    try:
        nCT = subprocess.check_output('''awk '{if (($12=="exonic") && (($4=="C" && $5=="T") || ($4=="G" && $5=="A"))) print}' Somatic_UHQ.tsv | wc -l''', shell=True)
    except subprocess.CalledProcessError:
        nCT = 0
    try:
        nTA = subprocess.check_output('''awk '{if (($12=="exonic") && (($4=="T" && $5=="A") || ($4=="A" && $5=="T"))) print}' Somatic_UHQ.tsv | wc -l''', shell=True)
    except subprocess.CalledProcessError:
        nTA = 0
    try:
        nTC = subprocess.check_output('''awk '{if (($12=="exonic") && (($4=="T" && $5=="C") || ($4=="A" && $5=="G"))) print}' Somatic_UHQ.tsv | wc -l''', shell=True)
    except subprocess.CalledProcessError:
        nTC = 0
    try:
        nTG = subprocess.check_output('''awk '{if (($12=="exonic") && (($4=="T" && $5=="G") || ($4=="A" && $5=="C"))) print}' Somatic_UHQ.tsv | wc -l''', shell=True)
    except subprocess.CalledProcessError:
        nTG = 0

    # transversions and transitions
    nTv = (int(nCA) + int(nCG) + int(nTA) + int(nTG))
    nTs = (int(nCT) + int(nTC))
    if nTv>0:
        Ts_Tv = '%3.2f' %(nTs/nTv)
    else:
        Ts_Tv = 0


    #------------------------------------------------------------------------
    #-----------Store stats in dictionary & write to json--------------------
    #------------------------------------------------------------------------

    #initialize
    summary_stats = OrderedDict() 
    summary_stats['analysis_date'] = str(datetime.date.today())
    summary_stats['run_mode']= run_mode
    summary_stats['sample_type'] = sample_type
    summary_stats['software_dir'] = softdir
    summary_stats['input_files'] = OrderedDict()
    summary_stats['filtering_params'] = OrderedDict()
    summary_stats['results'] = OrderedDict()
    summary_stats['results']['stats_filtering_tree'] = OrderedDict()
    summary_stats['results']['final_somatic_stats'] = OrderedDict()

    summary_stats['input_files']['tumor_bam'] = tumorbam
    summary_stats['input_files']['tumor_vcf'] = tumorvcf
    if run_mode == 'paired':
        summary_stats['input_files']['normal_vcf'] = normalvcf
    summary_stats['input_files']['reference_fasta'] = reffasta
    summary_stats['input_files']['bed_file'] = bedfile
    summary_stats['input_files']['blacklist_file'] = blackfile

    summary_stats['filtering_params']['variant allele depth'] = vardepth
    summary_stats['filtering_params']['p_fisher_exact_test'] = alpha1
    summary_stats['filtering_params']['p_strand_bias'] = pstbias
    summary_stats['filtering_params']['homopolymer_length_snp'] = hplen_snp
    summary_stats['filtering_params']['homopolymer_length_mnp'] = hplen_mnp
    summary_stats['filtering_params']['homopolymer_length_indel'] = hplen_indel 
    summary_stats['filtering_params']['ref_alt_read_Q_difference'] = diffq
    summary_stats['filtering_params']['alt_ref_MMQS_difference'] = diffmmq
    summary_stats['filtering_params']['ref_alt_read_length_difference'] = diffread1
    summary_stats['filtering_params']['expected_alt_read_length_difference'] = diffread2
    summary_stats['filtering_params']['normalized_variant_position'] = readpos
    summary_stats['filtering_params']['ref_alt_position_difference'] = delpos
    summary_stats['filtering_params']['distance_to_3p_end'] = dist3
    summary_stats['filtering_params']['counts_cosmic_for_whitelisting'] = ncoshit

    summary_stats['results']['stats_filtering_tree']['total_candidate_variants'] = nTotal
    summary_stats['results']['stats_filtering_tree']['annotated_germline'] = ngermline
    summary_stats['results']['stats_filtering_tree']['fail_fisher'] = nFishFail
    summary_stats['results']['stats_filtering_tree']['rescue_cosmic'] = nwhite1
    summary_stats['results']['stats_filtering_tree']['duplicate_vars'] = nDup
    summary_stats['results']['stats_filtering_tree']['before_read_QC'] = OrderedDict()
    summary_stats['results']['stats_filtering_tree']['before_read_QC']['total'] = int(nvar1) - 1
    summary_stats['results']['stats_filtering_tree']['before_read_QC']['whitelist_cosmic'] = nwhite
    summary_stats['results']['stats_filtering_tree']['before_read_QC']['mispriming'] = int(n_misprime)
    summary_stats['results']['stats_filtering_tree']['before_read_QC']['amplicon_bias'] = int(n_ampbias)
    summary_stats['results']['stats_filtering_tree']['fail_read_QC'] = OrderedDict()
    summary_stats['results']['stats_filtering_tree']['fail_read_QC']['total'] = int(nvar1)-int(nvar2) -nblack
    summary_stats['results']['stats_filtering_tree']['fail_read_QC']['pcr_filter'] = int(nPCR_error_final)
    summary_stats['results']['stats_filtering_tree']['fail_read_QC']['zero_vardepth'] = int(nzeroVarDepth)
    summary_stats['results']['stats_filtering_tree']['fail_read_QC']['hp'] = OrderedDict()
    summary_stats['results']['stats_filtering_tree']['fail_read_QC']['hp']['snp'] = int(n_hp_snp)
    summary_stats['results']['stats_filtering_tree']['fail_read_QC']['hp']['mnp'] = int(n_hp_mnp)
    summary_stats['results']['stats_filtering_tree']['fail_read_QC']['hp']['indel'] = int(n_hp_indel)
    summary_stats['results']['stats_filtering_tree']['fail_read_QC']['vardepth'] = int(nvarDepth)
    summary_stats['results']['stats_filtering_tree']['in_PON'] = nblack
    summary_stats['results']['stats_filtering_tree']['in_paired_normal'] = int(nvar2)-int(nvar3)

    summary_stats['results']['final_somatic_stats']['total'] = nFinal
    
    summary_stats['results']['final_somatic_stats']['exonic'] = OrderedDict()
    summary_stats['results']['final_somatic_stats']['exonic']['total']= int(nexonic)
    summary_stats['results']['final_somatic_stats']['exonic']['synonymous'] = int(nsyn)

    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous'] = OrderedDict()
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['total'] = int(nmiss)+int(nfsi)+int(nfsd)+int(nfss)++int(nstop)+int(nnfsi)+int(nnfsd)+int(nnfss) 
    
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['missense'] = OrderedDict()
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['missense']['total'] = int(nmiss)
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['missense']['damaging'] = ndamage
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['missense']['tolerant'] = ntol
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['missense']['mnp'] = nnadamage
    
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['frameshift'] = OrderedDict()
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['frameshift']['total'] = int(nfsi)+int(nfsd)+int(nfss)
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['frameshift']['insertion'] = int(nfsi)
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['frameshift']['deletion'] = int(nfsd)
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['frameshift']['substitution'] = int(nfss)

    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['nonframeshift'] = OrderedDict()
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['nonframeshift']['total'] = int(nnfsi)+int(nnfsd)+int(nnfss)
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['nonframeshift']['insertion'] = int(nnfsi)
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['nonframeshift']['deletion'] = int(nnfsd)
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['nonframeshift']['substitution'] = int(nnfss)
    summary_stats['results']['final_somatic_stats']['exonic']['nonsynonymous']['stopgain_stoploss'] = int(nstop) 
    
    summary_stats['results']['final_somatic_stats']['exonic']['snp_type'] = OrderedDict()
    summary_stats['results']['final_somatic_stats']['exonic']['snp_type']['C>A'] = int(nCA)
    summary_stats['results']['final_somatic_stats']['exonic']['snp_type']['C>G'] = int(nCG)
    summary_stats['results']['final_somatic_stats']['exonic']['snp_type']['C>T'] = int(nCT) 
    summary_stats['results']['final_somatic_stats']['exonic']['snp_type']['T>A'] = int(nTA)
    summary_stats['results']['final_somatic_stats']['exonic']['snp_type']['T>C'] = int(nTC)
    summary_stats['results']['final_somatic_stats']['exonic']['snp_type']['T>G'] = int(nTG)

    summary_stats['results']['final_somatic_stats']['exonic']['Ts_Tv'] = float(Ts_Tv)
    summary_stats['results']['final_somatic_stats']['splicing'] = int(nsplice)
    summary_stats['results']['final_somatic_stats']['intronic'] = int(nintronic)
    summary_stats['results']['final_somatic_stats']['cosmic'] = OrderedDict()
    summary_stats['results']['final_somatic_stats']['cosmic']['>0'] = ncos
    summary_stats['results']['final_somatic_stats']['cosmic']['>25'] = ncos25
    summary_stats['results']['final_somatic_stats']['cosmic']['>100'] = ncos100
    
    summary_stats['results']['final_somatic_stats']['maf_stats'] = OrderedDict()
    summary_stats['results']['final_somatic_stats']['maf_stats']['mean'] = float(uhq_maf_mean)
    summary_stats['results']['final_somatic_stats']['maf_stats']['stdev'] = float(uhq_maf_std)
    summary_stats['results']['final_somatic_stats']['maf_stats']['median'] = float(uhq_maf_median)
    #summary_stats['results']['final_somatic_stats']['maf_stats']['all_somatic']['mean'] = all_maf_mean
    #summary_stats['results']['final_somatic_stats']['maf_stats']['all_somatic']['stdev'] = all_maf_std
    #summary_stats['results']['final_somatic_stats']['maf_stats']['all_somatic']['median'] = all_maf_median


    out_json = 'Summary_Stats.json'
    with open(out_json, 'w') as f:
        json.dump(summary_stats, f, sort_keys=False, indent = 4)



    #-------------------------------------------------------------------
    #-----------------Write summary statistics--------------------------
    #-------------------------------------------------------------------

    #output file
    #summfile = open('Summary_Stats.txt' ,'w')
    #summfile.write("************************************************************************************\n")
    #summfile.write("*               Somatic Variant Calling from Ion Torrent Data                      *\n")
    #summfile.write("*               Version: 2.0                                                       *\n")
    #summfile.write("*               Author: Rahul K. Das                                               *\n")
    #summfile.write("*               Analysis Date: %s                                          *\n" %datetime.date.today())                              
    #summfile.write("************************************************************************************\n\n")
    #summfile.write("Run-mode: %s\n" %run_mode)
    #summfile.write("Sample-type: %s\n\n" %sample_type)
    #summfile.write("Software directory: %s\n\n" %softdir)
    #summfile.write("INPUT FILES....\n")
    #summfile.write("T-bam file: %s\n" %tumorbam)
    #summfile.write("T-VCF file: %s\n" %tumorvcf)
    #if run_mode == 'paired':
    #    summfile.write("N-VCF file: %s\n" %normalvcf)
    #summfile.write("Reference fasta file: %s\n" %reffasta)
    #summfile.write("BED file: %s\n" %bedfile)
    #summfile.write("Blacklist file: %s\n\n" %blackfile)
    #summfile.write("FILTERING PARAMETERS....\n")
    #summfile.write("p-value threshold for Fisher's exact test: %s\n" %alpha1)
    #summfile.write("p-value threshold for strand bias: %s\n" %pstbias)
    #summfile.write("Homopolymer length threshold: %s\n" %hplen)
    #summfile.write("Ref/Alt read quality difference threshold: %s\n" %diffq)
    #summfile.write("Alt/Ref mismatch quality difference threshold: %s\n" %diffmmq)
    #summfile.write("Ref/Alt read length difference threshold: |%s|\n" %diffread1)
    #summfile.write("Expected/Alt read length difference threshold: %s\n" %diffread2)
    #summfile.write("Average variant\'s position on the read threshold: %s\n" %readpos)
    #summfile.write("Ref/Alt position difference threshold: |%s|\n" %delpos)
    #summfile.write("Threshold for closest distance to 3\'-end of read: %s\n" %(dist3))
    #summfile.write("Threshold for no. of cosmit hits to whitelist: %s\n" %ncoshit)
    #summfile.write("************************************************************************************\n\n")
    #summfile.write("NODE STATS OF FILTERING NETWORK....\n")
    #summfile.write('Total potential somatic variants based on hard AF cutoff: %s\n' %str(nTotal))
    #summfile.write("Total known germline variants that were filtered out (1000g, dbsnp): %s\n" %ngermline)
    #summfile.write("Total variants that failed Fisher's exact test: %s\n" %str(nFishFail))
    #summfile.write("Total variants that did not pass QC metrics (Str. bias, HP lengths etc..): %s\n" %(int(nvar1)-int(nvar2)))
    #summfile.write("Total variants that were filtered out for MAF>0 in paired/unpaired N: %s\n" %(int(nvar2)-int(nvar3))) 
    #summfile.write("Total variants that were whitelisted due to cosmic entries: %s\n" %str(nwhite))
    #summfile.write("************************************************************************************\n\n")
    #summfile.write("SUMMARY OF ULTRA_HIGH_QUALITY SOMATIC VARIANTS....\n")
    #summfile.write('Total=%s\n' %str(nFinal))
    #summfile.write('Exonic=%s\n' %int(nexonic))
    #summfile.write('    Missense_SNP_MNP=%s;Damaging=%s,Tolerant=%s,MNP=%s\n' %(int(nmiss),ndamage,ntol,nnadamage))
    #summfile.write('    Frameshift=%s;Insertion=%s,Deletion=%s,Substitution=%s\n' %(int(nfsi)+int(nfsd)+int(nfss),int(nfsi),int(nfsd),int(nfss)))
    #summfile.write('    NonFrameshift=%s;Insertion=%s,Deletion=%s,Substitution=%s\n' %(int(nnfsi)+int(nnfsd)+int(nnfss),int(nnfsi),int(nnfsd),int(nnfss)))
    #summfile.write('    Stopgain_Stoploss=%s\n' %int(nstop))
    #summfile.write('    Total_Nonsynonymous=%s\n' %(int(nmiss)+int(nfsi)+int(nfsd)+int(nfss)++int(nstop)+int(nnfsi)+int(nnfsd)+int(nnfss)))
    #summfile.write('C>A=%s;C>G=%s;C>T=%s\n' %(n_CA, n_CG, n_CT))
    #summfile.write('T>A=%s;T>C=%s;T>G=%s\n' %(n_TA, n_TC, n_TG))
    #summfile.write('Ts=%s;Tv=%s;Ts/Tv=%s\n' %(nTs, nTv, Ts_Tv))
    #summfile.write('Splicing=%s\n' %int(nsplice))
    #summfile.write('Intronic=%s\n' %int(nintronic))
    #summfile.write("Cosmic_entries=%s;atleast 50=%s\n\n" %(ncos,ncos25))
    #summfile.write('VARIANTS MAF STATS IN THE T SAMPLES....\n')
    ##if len(germ_maf) > 0:
    ##    summfile.write('Germline:Mean=%3.2f;Std=%3.2f;Median=%3.2f\n' %(np.mean(germ_maf),np.std(germ_maf),np.median(germ_maf)))
    ##else:
    ##    summfile.write('Germline:Mean=0;Std=0;Median=0\n')
    ##if len(all_maf) > 0:
    ##    summfile.write('All_Somatic:Mean=%3.2f;Std=%3.2f;Median=%3.2f\n' %(np.mean(all_maf),np.std(all_maf),np.median(all_maf)))
    ##else:
    ##    summfile.write('All_Somatic:Mean=0;Std=0;Median=0\n')
    #if len(uhq_maf) > 0:
    #    summfile.write('UHQ_Somatic:Mean=%3.2f;Std=%3.2f;Median=%3.2f\n' %(np.mean(uhq_maf),np.std(uhq_maf),np.median(uhq_maf)))
    #else:
    #    summfile.write('UHQ_Somatic:Mean=0;Std=0;Median=0\n')
    #
    #
    #summfile.close()
    
