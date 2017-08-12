__Author__ =	"Rahul K. Das"
__Date__ = "May 4, 2016"
__Version__ = "2.0"
__LastModified__ = "Feb 23, 2017"



import argparse
import os, sys, datetime
from processor1 import *
from processor2 import *
from process_variants import *
from varfile_parser import *



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--version', action='version', version='%(prog)s 2.0')
	parser.add_argument('-v', '--vcf', dest='vcf', help='VCF from TVC', metavar='')
	parser.add_argument('-b', '--bam', dest='bam', help='Final PTRIM bam', metavar='')
	parser.add_argument('-o', '--outdir', dest='outdir', help='Output directory', metavar='')
	parser.add_argument('-s', '--softdir', dest='softdir', help='Software directory', metavar='')
	parser.add_argument('-N', '--annovar', dest='annovar', help='Annovar directory path', metavar='')
	parser.add_argument('-c', '--bamrc', dest='bamrc', help='bam-readcount executable path', metavar='')
	parser.add_argument('-f', '--reffasta', dest='reffasta', help='Reference fasta path', metavar='')
	parser.add_argument('-A', '--alpha', dest='alpha', type=float, default=0.1, help='p-value threshold for Bionomial Tests (Default=0.1)', metavar='')
	parser.add_argument('-p', '--pstbias', dest='pstbias', type=float, default=1.0, help='p-value threshold for Strand Bias (Default=1.0)', metavar='')
	parser.add_argument('-l', '--hplen_indel', dest='hplen_indel', type=int, default=2, help='homopolymer length threshold for indels(Default=2); longer indels will be subjected to MAF check based on 1000g or exome6500 databases ', metavar='')
	parser.add_argument('-L', '--hplen_snp', dest='hplen_snp', type=int, default=3, help='homopolymer length threshold for snps(Default=3); longer snps will be subjected to MAF check based on 1000g or exome6500 databases ', metavar='')
        parser.add_argument('-G', '--hplen_mnp', dest='hplen_mnp', type=int, default=2, help='homopolymer length threshold for mnps (Default=2)', metavar='')
        parser.add_argument('-D', '--diffq', dest='diffq', type=int, default=30, help='Map Quality Difference threshold (default=30)', metavar='')
	parser.add_argument('-M', '--diffmmq', dest='diffmmq', type=int, default=45, help='MMQS Difference threshold (Default=45)', metavar='')
	parser.add_argument('-r', '--diffread', dest='diffread', type=int, default=30, help='|Ref/var Read Length Difference| threshold (Default=30)', metavar='')
	parser.add_argument('-x', '--readpos', dest='readpos', type=float, default=0.1, help='Variant Read Position threshold (Default=0.1)', metavar='')
	parser.add_argument('-d', '--dist3', dest='dist3', type=float, default=0.1, help='Distance from 3\'-end threshold (Default=0.1)', metavar='')
	parser.add_argument('-C', '--contextLen', dest='contextLen', type=int, default=5, help='the desired length of sequence context', metavar='')
        parser.add_argument('-P', '--ncpu', dest='ncpu', type=int, default=1, help='number of CPUs requested')
        parser.add_argument('-B', '--bedfile', dest='bedfile', help='the BED file', metavar='')
        parser.add_argument('-Y', '--bedHeader', dest='bedHeader', choices=['yes','no'], help='header present in BED file', metavar='')
        #parser.add_argument('-T', '--mindepth', dest='mindepth', type=int,  default=30, help='Total depth threshold')	

	results = parser.parse_args()
	
	vcf = results.vcf
	bam = results.bam
	outdir = results.outdir
	softdir = results.softdir
	annovar = results.annovar
	bamrc = results.bamrc
	reffasta = results.reffasta
	alpha = float(results.alpha)
	pstbias = float(results.pstbias)
	hplen_indel = int(results.hplen_indel)
        hplen_snp = int(results.hplen_snp)
        hplen_mnp = int(results.hplen_mnp)
	diffq = int(results.diffq)
	diffmmq = int(results.diffmmq)
	diffread = int(results.diffread)
	readpos = float(results.readpos)
	dist3 = float(results.dist3)
	contextLen = int(results.contextLen)
        ncpu = int(results.ncpu)
        bedfile = results.bedfile
        bedHeader = results.bedHeader
        #mindepth = int(results.mindepth)


	print time.asctime( time.localtime(time.time()) ) + ": Workflow started...."
	#go to output directory, and output everything there
	workdir = os.path.join(outdir,'Germline_Variants_QC')
	#if os.path.exists(workdir):
	#	shutil.rmtree(workdir)
	if not os.path.exists(workdir):
		os.mkdir(workdir)
	os.chdir(workdir)
	

	##get the list of all variants and make the vcf with only variants
	print " Extracting all variants...."
	varListAll = extract_variants(vcf)

	#annotate by Annovar
	print " Getting annovar annotations of all variants in %s locus...." %len(varListAll)
	varvcf = 'Master.vcf'
	annoListAll, mnp_idx, mnp_AAChange, cosmic = annovar_annotate(softdir,annovar,varvcf,ncpu)
	
	#process
	print " Generating QC metrics & Processing variants...."
	nstart,nfilter = process_variants(bam, varListAll, annoListAll, mnp_idx, mnp_AAChange, cosmic, bamrc, reffasta,alpha, pstbias, hplen_indel, hplen_snp, hplen_mnp, diffq, diffmmq, diffread, readpos, dist3, contextLen, ncpu, bedfile, bedHeader)

	nint = subprocess.check_output('''grep -w "intronic" Germline_HQ.tsv | wc -l''', shell=True)
	nexo = subprocess.check_output('''grep -w "exonic" Germline_HQ.tsv | wc -l''', shell=True)	
        nspl = subprocess.check_output('''grep -w "splicing" Germline_HQ.tsv | wc -l''', shell=True)
	nsyn = subprocess.check_output('''grep -w "synonymous_SNV" Germline_HQ.tsv | wc -l''', shell=True)
	nnsyn = subprocess.check_output('''grep -w "nonsynonymous_SNV" Germline_HQ.tsv | wc -l''', shell=True)
	nhet = subprocess.check_output('''grep -w "0/1" Germline_HQ.tsv | wc -l''', shell=True)
	nhom = subprocess.check_output('''grep -w "1/1" Germline_HQ.tsv | wc -l''', shell=True)


	nzero = subprocess.check_output('''awk '{if ($6==0) print }' Germline_Removed.tsv | wc -l''', shell=True)
	
	
	#-----------------Write summary statistics------------------------------
	summfile = open('Summary.txt','w')	  
	summfile.write("************************************************************************************\n")
	summfile.write("*				HQ Germline Variants from Ion Torrent Data                         *\n")
	summfile.write("*				Version: 2.0													   *\n")
	summfile.write("*				Author: Rahul K. Das											   *\n")
	summfile.write("*				Analysis Date: %s										   *\n" %datetime.date.today())
	summfile.write("************************************************************************************\n\n")
	summfile.write("Software directory: %s\n" %softdir)
	summfile.write("Output Dir: %s\n\n" %outdir)
	summfile.write("INPUT FILES....\n")
	summfile.write("bam file: %s\n" %bam)
	summfile.write("VCF file: %s\n" %vcf)
	summfile.write("Reference fasta file: %s\n\n" %reffasta)
	summfile.write("FILTERING PARAMETERS....\n")
	summfile.write("p-value threshold for Binomial test: %s\n" %alpha)
	summfile.write("p-value threshold for strand bias: %s\n" %pstbias)
	summfile.write("Homopolymer length threshold (indels): %s\n" %hplen_indel)
	summfile.write("Homopolymer length threshold (snps): %s\n" %hplen_snp)
        summfile.write("Homopolymer length threshold (mnps): %s\n" %hplen_mnp)
        summfile.write("Ref/Alt read quality difference threshold: %s\n" %diffq)
	summfile.write("Alt/Ref mismatch quality sum difference threshold: %s\n" %diffmmq)
	summfile.write("Ref/Alt read length difference threshold: |%s|\n" %diffread)
	summfile.write("Average variant\'s position on the read threshold: %s\n" %readpos)
	summfile.write("Threshold for closest distance to 3\'-end of read: %s\n\n" %(dist3))
	summfile.write("************************************************************************************\n\n")
	summfile.write("SUMMARY OF GERMLINE VARIANTS....\n")	
	summfile.write("Total variants from VCF: %s\n" %nstart)
	summfile.write("Total variants filtered out: %s\n" %nfilter)
	summfile.write("	zero variant alleles: %s\n" %nzero)
	summfile.write("Total variants that passed QC: %s\n" %(int(nstart)-int(nfilter)))
	summfile.write("Exonic: %s\n" %nexo.strip())
	summfile.write("	Synonymous: %s\n" %nsyn.strip())
	summfile.write("    Nonsynonymous_SNP: %s\n" %nnsyn.strip())
        summfile.write("Splicing: %s\n" %nspl.strip())
        summfile.write("Intronic: %s\n" %nint.strip())
	summfile.write("HET: %s\n" %nhet.strip())
	summfile.write("HOM: %s\n" %nhom.strip())

	summfile.close()

	print time.asctime( time.localtime(time.time()) )+": Finished Successfully...."
