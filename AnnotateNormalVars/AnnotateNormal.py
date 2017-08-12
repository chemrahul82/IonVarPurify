#!/usr/bin/env python

__Author__ =    "Rahul K. Das"
__Date__ = "May 4, 2016"
__Version__ = "2.0"
__LastModified__ = "Jan 23, 2017"



import argparse
import os, sys, datetime
from processor1 import *
from processor2 import *
from process_variants import *
from varfile_parser import *



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version', version='%(prog)s 2.0')
    parser.add_argument('-N', '--normalvcf', dest='normalvcf', help='Normal VCF from TVC', metavar='')
    parser.add_argument('-b', '--normalbam', dest='normalbam', help='Final Normal bam', metavar='')
    parser.add_argument('-o', '--outdir', dest='outdir', help='Output directory', metavar='')
    parser.add_argument('-A', '--annovar', dest='annovar', help='Annovar directory path', metavar='')
    parser.add_argument('-s', '--softdir', dest='softdir', help='Software source directory', metavar='')
    parser.add_argument('-c', '--bamreadcount', dest='bamreadcount', help='bam-readcount executable path', metavar='')
    parser.add_argument('-f', '--reffasta', dest='reffasta', help='Reference fasta path', metavar='')
    parser.add_argument('-P', '--ncpu', dest='ncpu', type=int, default=1, help='number of CPUs requested')


    results = parser.parse_args()
    
    normalvcf = results.normalvcf
    normalbam = results.normalbam
    outdir = results.outdir
    softdir = results.softdir
    annovar = results.annovar
    bamreadcount = results.bamreadcount
    reffasta = results.reffasta
    ncpu = results.ncpu

    print time.asctime( time.localtime(time.time()) ) + ": Workflow started...."
    #go to output directory, and output everything there
    os.chdir(outdir)

    # create the temp directory
    if not os.path.exists('_tmp_'):
       os.mkdir('_tmp_')
    
    ##get the list of all variants and make the vcf with only variants
    print " Extracting all variants...."
    # keep only the MAF>0 variants
    varListAll = extract_variants(normalvcf)
    
    #annotate by Annovar
    print " Getting annovar annotations of all variants...."
    varvcf = 'AllVars.vcf'
    annoListAll, mnp_idx, mnp_AAChange, cosmic = annovar_annotate(softdir,annovar,varvcf,ncpu)
    
    #process
    print " Processing variants and generating list of potential error positions...."
    process_variants(normalbam, varListAll, annoListAll, mnp_idx, mnp_AAChange, cosmic, bamreadcount, reffasta, ncpu)

    # move all files except summary and variant files to tmp directory
    command = "ls | egrep -v '(%s|%s|%s|%s)' | xargs -i mv {} %s/" %('AllVars.vcf','Summary.txt','N_Vars.tsv','_tmp_','_tmp_')
    os.system(command)

    print time.asctime( time.localtime(time.time()) )+": Finished Successfully...."
    
    
    #-----------------Write summary statistics------------------------------
    summfile = open('Summary.txt','w')    
    summfile.write("************************************************************************************\n")
    summfile.write("*        Annotation of Non-cosmic & Non-germline variants from Ion Torrent Data    *\n")
    summfile.write("*               Version: 2.0                                                       *\n")
    summfile.write("*               Author: Rahul K. Das                                               *\n")
    summfile.write("*               Analysis Date: %s                                          *\n" %datetime.date.today())
    summfile.write("************************************************************************************\n\n")
    summfile.write("INPUT FILES....\n")
    summfile.write("N-bam file: %s\n" %normalbam)
    summfile.write("N-VCF file: %s\n" %normalvcf)
    summfile.write("Output Dir: %s\n" %outdir)

    summfile.close()
