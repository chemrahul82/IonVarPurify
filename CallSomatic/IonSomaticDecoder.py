#!/usr/bin/env python

"""
Parse the input arguments from command line and
run IonSomaticDecoder

"""

__Author__ =    "Rahul K. Das"
__Date__ = "May 4, 2016"
__Version__ = "2.0"
__LastUpdated__ = "Feb 25, 2017"

import argparse
import os, sys
import time
import shutil

from varfile_parser import *
from processor1 import *
#from process_somatic import *
from process_somatic_new import *
from write_summary import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 2.0')
    parser.add_argument('-q', '--readq', dest='readq', type=int, default=30, help='Read quality threshold (Default=30)', metavar='')
    parser.add_argument('-Q', '--baseq', dest='baseq', type=int, default=20, help='Base quality threshold (Default=20)', metavar='')
    parser.add_argument('-a', '--alpha1', dest='alpha1', type=float, default=0.005, help='p-value threshold for Fisher\'s Exact Test (Default=0.005)', metavar='')
    parser.add_argument('-p', '--pstbias', dest='pstbias', type=float, default=1.0, help='p-value threshold for Strand Bias (Default=1.0)', metavar='')
    parser.add_argument('-l', '--hplen_indel', dest='hplen_indel', type=int, default=2, help='homopolymer length threshold for indels (Default=2)', metavar='')
    parser.add_argument('-L', '--hplen_snp', dest='hplen_snp', type=int, default=3, help='homopolymer length threshold for snps (Default=3)', metavar='')
    parser.add_argument('-G', '--hplen_mnp', dest='hplen_mnp', type=int, default=3, help='homopolymer length threshold for mnps (Default=3)', metavar='')
    parser.add_argument('-D', '--diffq', dest='diffq', type=int, default=20, help='Map Quality Difference threshold (default=20)', metavar='')
    parser.add_argument('-M', '--diffmmq', dest='diffmmq', type=int, default=40, help='MMQS Difference threshold (Default=40)', metavar='')
    parser.add_argument('-r', '--diffread1', dest='diffread1', type=int, default=25, help='|Ref/var Read Length Difference| threshold (Default=25)', metavar='')
    parser.add_argument('-R', '--diffread2', dest='diffread2', type=int, default=35, help='Expected/var Read Length Difference threshold (Default=35)', metavar='')
    parser.add_argument('-x', '--readpos', dest='readpos', type=float, default=0.1, help='Variant Read Position threshold (Default=0.1)', metavar='')
    parser.add_argument('-X', '--delpos', dest='delpos', type=float, default=0.25, help='|Ref/Var Read Position Difference| threshold (Default=0.25)', metavar='')
    parser.add_argument('-d', '--dist3', dest='dist3', type=float, default=0.1, help='Distance from 3\'-end threshold (Default=0.1)', metavar='')
    parser.add_argument('-H', '--vardepth', dest='vardepth', type=int, default=2, help='Variant read depth threshold (Default=2)', metavar='')
    parser.add_argument('-F', '--varfreq', dest='varfreq', type=float, default=0.0, help='Variant AF threshold (Default=0.0)', metavar='')
    parser.add_argument('-W', '--ncoshit', dest='ncoshit', type=int, default=25, help='Threshold for no. of cosmic hits to whitlist a variant (Default=25)', metavar='')
    parser.add_argument('-T', '--tumorvcf', dest='tumorvcf', help='Tumor VCF from TVC', metavar='')
    parser.add_argument('-N', '--normalvcf', dest='normalvcf', help='Paired Normal VCF from TVC', metavar='')
    parser.add_argument('-b', '--tumorbam', dest='tumorbam', help='Final Tumor bam', metavar='')
    parser.add_argument('-o', '--outdir', dest='outdir', help='Output directory', metavar='') 
    parser.add_argument('-s', '--softdir', dest='softdir', help='Software source directory', metavar='')
    parser.add_argument('-c', '--bamrcpath', dest='bamrcpath', help='bam-readcount executable path', metavar='')
    parser.add_argument('-A', '--annovar', dest='annovar', help='Annovar directory path', metavar='')
    parser.add_argument('-f', '--reffasta', dest='reffasta', help='Reference fasta path', metavar='')
    parser.add_argument('-B', '--bedfile', dest='bedfile', help='the BED file', metavar='')
    parser.add_argument('-Y', '--bedHeader', dest='bedHeader', choices=['yes','no'], help='header present in BED file', metavar='')
    parser.add_argument('-k', '--blackfile', dest='blackfile', default='', help='the Panel of Normal/ blacklisted variants', metavar='')    
    parser.add_argument('-C', '--contextLen', dest='contextLen', type=int, default=5, help='the desired length of sequence context', metavar='')
    parser.add_argument('-U', '--sample_type', dest='sample_type', choices=['fresh','ffpe','unknown'], default='unknown', help='the sample source: fresh/ffpe', metavar='')
    parser.add_argument('-Z', '--run_mode', dest='run_mode', choices=['paired','unpaired'], default='unpaired', help='the run mode: paired/unpaired', metavar='')
    parser.add_argument('-E', '--varcall_mode', dest='varcall_mode', choices=['hotspot','not_hotspot'], default='not_hotspot', help='the variant calling mode: hotspot/not_hotspot', metavar='')
    parser.add_argument('-P', '--ncpu', dest='ncpu', type=int, default=1, help='number of CPUs requested')
    
    options = parser.parse_args()
    

    readq = int(options.readq)
    baseq = int(options.baseq)
    alpha1 = float(options.alpha1)
    pstbias = float(options.pstbias)
    hplen_indel = int(options.hplen_indel)
    hplen_snp = int(options.hplen_snp)
    hplen_mnp = int(options.hplen_mnp)
    diffq = int(options.diffq)
    diffmmq = int(options.diffmmq)
    diffread1 = int(options.diffread1)
    diffread2 = int(options.diffread2)
    readpos = float(options.readpos)
    delpos = float(options.delpos)
    dist3 = float(options.dist3)
    vardepth = int(options.vardepth)
    varfreq = float(options.varfreq)
    ncoshit = int(options.ncoshit)
    normalvcf = options.normalvcf
    tumorvcf = options.tumorvcf
    tumorbam = options.tumorbam
    outdir = options.outdir
    softdir = options.softdir
    bedfile = options.bedfile 
    bedHeader = options.bedHeader
    blackfile = options.blackfile
    reffasta = options.reffasta
    bamrcpath = options.bamrcpath
    annovar = options.annovar
    contextLen = options.contextLen
    sample_type = options.sample_type
    run_mode = options.run_mode
    varcall_mode = options.varcall_mode
    ncpu = options.ncpu
    
    #sanity checks
    if run_mode == 'paired' and not options.normalvcf:
        print "ERROR:paired run mode requested but the normal VCF is not provided!"
        sys.exit()
    if run_mode == 'unpaired' and options.normalvcf is not None:
        print "WARNING: unpaired run mode requested, the normal VCF will be ignored!"
    if not options.bedHeader:
        print "ERROR: Whether the project bed file has header or not??"
        sys.exit()
    
    print time.asctime( time.localtime(time.time()) ) + ": Workflow started..."
    #go to output directory, make directory and and output everything there
    workdir = os.path.join(outdir,'Somatic_Variants_UHQ')
    #if os.path.exists(workdir):
    #    shutil.rmtree(workdir)
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    else:
        print "WARNING! Working directory exists, not overwriting; If the previous run was in different mode, delete the directory"
    os.chdir(workdir)

    # create the temp directory, where all files except the summary file and UHQ file will move
    if not os.path.exists('_tmp_'): 
        os.mkdir('_tmp_')


    #get the list of all unfiltered somatic vars
    print "Extracting raw unfiltered somatic variants...."
    varListAll = extract_raw_somatic(normalvcf, tumorvcf, run_mode, varcall_mode)
    print "%s positions extracted for evaluation..." %len(varListAll)

    print " Getting annotations of unfiltered potential somatic variants..."
    annoListAll, mnp_idx, mnp_AAChange, cosmic = annovar_annotate(softdir, annovar, 'Somatic_Unfiltered.vcf', ncpu)
   
    print "Started false positive filtration process..."
    ngermline, nFishFail, nwhite1, nwhite, nblack, nDup, n_misprime, n_ampbias, nPCR_error_final, nvarDepth, nzeroVarDepth, n_hp_snp, n_hp_indel, n_hp_mnp, all_maf = process_somatic(tumorvcf, tumorbam, bedfile, bedHeader, blackfile, softdir, bamrcpath, reffasta, varListAll, annoListAll, \
            mnp_idx, mnp_AAChange, cosmic, readq, baseq, alpha1, pstbias, hplen_indel, hplen_snp, hplen_mnp, diffq, diffmmq, diffread1, diffread2, readpos, delpos, dist3, vardepth, varfreq, ncoshit, contextLen, run_mode, ncpu)
  
    write_summary(varListAll,ngermline, nFishFail, nwhite1, nwhite, nblack, nDup, n_misprime, n_ampbias, nPCR_error_final, nvarDepth, nzeroVarDepth, n_hp_snp, n_hp_indel, n_hp_mnp, all_maf, normalvcf, tumorvcf, tumorbam, reffasta, bedfile, blackfile, softdir, \
            readq, baseq, alpha1, pstbias, hplen_indel, hplen_snp, hplen_mnp, diffq, diffmmq, diffread1, diffread2, readpos, delpos, dist3, vardepth, varfreq, ncoshit, contextLen, run_mode, sample_type)

    # move all files except summary and UHQ variant files to tmp directory
    command = "ls | egrep -v '(%s|%s|%s|%s|%s)' | xargs -i mv {} %s/" %('Somatic_UHQ.tsv','Somatic_exonic_splicing_UHQ.tsv','Summary_Stats.json','Summary_Stats.txt','_tmp_','_tmp_')
    os.system(command)
    print time.asctime( time.localtime(time.time()) )+": Finished Successfully..."

    
