"""
Analysis of all variant positions in the N samples
1.calculate bam-readcount metrics
2.perform binomial test
write a file with high-confident germline variants with annotations

"""


from __future__ import division

__Author__ =	"Rahul K. Das"
__Date__ = "May 14, 2016"
__Version__ = "2.0"
__LaseModified__ = "Feb 28, 2017"


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
from cal_hplen import *
from qc_pcr_artifacts import *
from cal_exp_readlen import *


def process_variants(normalbam, varListAll, annoListAll, mnp_idx, mnp_AAChange, cosmic, bamreadc, ref_fasta,\
	alpha, p_stbias, hplen_indel, hplen_snp, hplen_mnp, diffq, diffmmq, diffread, readpos, dist3, contextLen, ncpu, bedfile, bedHeader):


	"""
	Get the bam-readcount stats for only non-cosmic & non-germline variant positions in N samples and make
	multiple files
	"""

	# create output files
	var_anno_file = open('Germline_Dups.tsv', 'w') 
	
	# header for output files
	header = ['Chr','Pos','Context','Ref','Alt','VAF','Alt_Depth','Total_Depth','GT','Variant_Function','Variant_Gene','Variant_Exonic_Function','AA_Change',\
			'1000gAug15','dbsnp138','esp6500siv2','Cosmic_Counts',\
			'p_Strand_Bias', 'HP_Lenth', 'Alt_Q', 'Diff_Q', 'Diff_MMQS', 'Diff_Read_Length','Position', 'Distance_to_3p','Functional_Effect','Clin_Sig','Clin_Dbn']
	var_anno_file.write('\t'.join(header)+'\n')

	
	# Remove known germline
	# varList, annoList, idxList, ngermline, germ_maf = removeGerm(varListAll, annoListAll)

	
	#remove the cosmic entries with counts > 25 (keep this count a bit high to
	# filter out FPs that were present in paired N samples but were kept for being cosmic entry;
	# if these cosmic variants appear in multiple N samples, these are sequencing error and not
	# that the N sample was contaminated with the T

	#newAnnoList = []; newidxList = []
	#for idx, el in enumerate(annoList):
	#	 if el.split('\t')[13] == '.':
	#		 #non-cosmic
	#		 newAnnoList.append(el)
	#		 newidxList.append(idxList[idx])
	#	 else:
	#		 #cosmic with < n counts
	#		 cosmicHits = el.split('\t')[13].split(';')[1].split('=')[1].split(',')
	#		 counts = 0
	#		 for item in cosmicHits:
	#			 start = item.find( '(' )
	#			 counts += int(item[:start])
	#		 if counts < 5:
	#			 newAnnoList.append(el)
	#			 newidxList.append(idxList[idx])

	varList = varListAll
	annoList = annoListAll
	idxList = range(len(annoListAll))
	
	
        # get the Damaging/Tolerant stats for these variants
        DStats = damagingStats(annoList,idxList, mnp_idx)
        
        # get metrics by running bam-readcount on Normal bam files
	alt_maf, alt_depth, depth, alt_q, diff_q, diff_mmqs, diff_readlen, pos, dist3p = processor2(varList, annoList, normalbam, bamreadc, ref_fasta, ncpu)


	# get the homopolymer length around the variants
	with open('temp.var','w') as f:
		for el in annoList:
			f.write('\t'.join(el.split('\t')[:2]) + '\n')
        
        
        motifList1 = getSequence('temp.var',ref_fasta,5)
	scanLength = 20
        motifList2 = getSequence('temp.var',ref_fasta,scanLength) #for HPLen
	os.remove('temp.var')


	# write all the variants with all the annotations	 
	for idx, el in enumerate(annoList):

		vcffields = el.split('\t')[54].split(';')
		vcfDict = {}
		#strand bias and hplen
		for ee in vcffields:
			if len(ee.split('=')) == 2:
				vcfDict[ee.split('=')[0]] = ee.split('=')[1]
		if 'STBP' in vcfDict:
			stbp = vcfDict['STBP']
		else:
			stbp = '.'
			#if 'HRUN' in ee:
			#	hrun = ee.split('=')[1]
		
		#genotype
		#store the GT
		if 'GT' in el.split('\t')[55].split(':'):
			gt = el.split('\t')[56].split(':')[0]
		else:
			gt = '.'
		
                
                ref = el.split('\t')[3]
                alt = el.split('\t')[4]
                
                ##Sequence Context##
                motifseq = motifList1[idx][:contextLen] + motifList1[idx][contextLen].lower() + motifList1[idx][-contextLen:]

		#compute the homopolymer length
                if len(ref) < 15:
                    hrun = cal_hplen(motifList2[idx],scanLength,ref,alt)
                else:
                    hrun = '-1'

                #hrun = str(max(sum(1 for i in g) for k,g in groupby(motifList2[idx])))


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
		if el.split('\t')[5] == 'exonic' and el.split('\t')[9] != 'UNKNOWN'\
			and idxList[idx] not in mnp_idx:
			#variants that have aa change info in annovar 1st run
			if len(el.split('\t')[9].split(',')[0].split(':')[-1].split('.')) == 2:
				var_anno_file.write('\t'.join(el.split('\t')[:2]+\
					[motifseq]+\
                                        el.split('\t')[3:5]+\
					map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
					[gt]+\
					el.split('\t')[5:7]+\
					['_'.join(el.split('\t')[8].split(' '))]+\
					[el.split('\t')[9].split(',')[0].split(':')[-1].split('.')[1]]+\
					el.split('\t')[10:13]+\
					[str(counts)]+[stbp,hrun]+\
					map(str, [alt_q[idx], diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx], DStats[idx], el.split('\t')[14], el.split('\t')[15]]))+'\n')
		
			#variants that don't have aa change info in both annovar's 1st & 2nd run
			else:
				var_anno_file.write('\t'.join(el.split('\t')[:2]+\
					[motifseq]+\
                                        el.split('\t')[3:5]+
					map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
					[gt]+\
					el.split('\t')[5:7]+\
					['_'.join(el.split('\t')[8].split(' '))]+\
					['.']+\
					el.split('\t')[10:13]+\
					[str(counts)]+[stbp,hrun]+\
					map(str,[alt_q[idx], diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx], DStats[idx], el.split('\t')[14], el.split('\t')[15]]))+'\n')
		
		elif el.split('\t')[5] == 'exonic' and el.split('\t')[9] == 'UNKNOWN'\
			and idxList[idx] not in mnp_idx:
			var_anno_file.write('\t'.join(el.split('\t')[:2]+\
				[motifseq]+\
                                el.split('\t')[3:5]+\
				map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
				[gt]+\
				el.split('\t')[5:7]+\
				['_'.join(el.split('\t')[8].split(' '))]+\
				[el.split('\t')[9]]+\
				el.split('\t')[10:13]+\
				[str(counts)]+[stbp,hrun]+\
				map(str,[alt_q[idx], diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx], DStats[idx], el.split('\t')[14], el.split('\t')[15]]))+'\n')
					
		# mnps and nfs substitutions for which annovar was run 2nd time
		elif el.split('\t')[5] == 'exonic' and idxList[idx] in mnp_idx:
			# mnps
			if len(el.split('\t')[3]) == len(el.split('\t')[4]):
				var_anno_file.write('\t'.join(el.split('\t')[:2]+\
				[motifseq]+\
                                el.split('\t')[3:5]+\
				map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
				[gt]+\
				el.split('\t')[5:7]+\
				['nonsynonymous_MNV']+\
				mnp_AAChange[mnp_idx.index(idxList[idx])]+\
				el.split('\t')[10:13]+\
				[str(counts)]+[stbp,hrun]+\
				map(str,[alt_q[idx], diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx], DStats[idx], el.split('\t')[14], el.split('\t')[15]]))+'\n')
			else:
				#nfs indels
				var_anno_file.write('\t'.join(el.split('\t')[:2]+\
				[motifseq]+\
                                el.split('\t')[3:5]+\
				map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
				[gt]+\
				el.split('\t')[5:7]+\
				['nonframeshift_substitution']+\
				mnp_AAChange[mnp_idx.index(idxList[idx])]+\
				el.split('\t')[10:13]+\
				[str(counts)]+[stbp,hrun]+\
				map(str,[alt_q[idx], diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx], DStats[idx], el.split('\t')[14], el.split('\t')[15]]))+'\n')
		else:
			 var_anno_file.write('\t'.join(el.split('\t')[:2]+\
				[motifseq]+\
                                el.split('\t')[3:5]+\
				map(str,[alt_maf[idx],alt_depth[idx],depth[idx]])+\
				[gt]+\
				el.split('\t')[5:7]+\
				['.','.']+\
				el.split('\t')[10:13]+\
				[str(counts)]+[stbp,hrun]+\
				map(str,[alt_q[idx], diff_q[idx], diff_mmqs[idx], diff_readlen[idx], pos[idx], dist3p[idx], DStats[idx], el.split('\t')[14], el.split('\t')[15]]))+'\n')
		
		
	var_anno_file.close()

	
	#-------------Remove Duplicate Variants that may appear due to semi-duplicate entries in VCF file"
	os.system("awk '!x[$0]++' %s > %s" %('Germline_Dups.tsv','Germline_noDups.tsv'))

	
	#-------------Change Diff_Q, Diff_MMQS, Diff_Read_Length columns to '.' for HOMs with AF = 1 in T
	oldtempfile = open('Germline_noDups.tsv','r')
	newtempfile = open('Germline_All.tsv','w')

	for il,line in enumerate(oldtempfile):
		# write header
		if il == 0:
			newtempfile.write(line)
		elif il > 0:
                        # fix for HOM: 
                        # if altQ = |diff_Q| or
                        #if ref allele depth is very low, doesn't make any sense to calculate statistics based on that
                        if (line.split('\t')[19] != '.' and float(line.split('\t')[19]) == abs(float(line.split('\t')[20]))):
                            newline = '\t'.join(line.split('\t')[:19]+['.','.','.','.']+line.split('\t')[-5:])
			    newtempfile.write(newline)

                        elif ((float(line.split('\t')[7]) - float(line.split('\t')[6])) <= 10):
                            newline = '\t'.join(line.split('\t')[:19]+['.','.','.','.']+line.split('\t')[-5:])
                            newtempfile.write(newline)

			else:
			    newtempfile.write(line)

	oldtempfile.close(); newtempfile.close()


        #-----------------PCR errors: Mispriming or amplicon bias---------------------------------------
        #os.system('''awk '{if (NR>1) print $1"\t"$2"\t"$4"\t"$5}' Germline_All.tsv > temp1.region''')
        #pass_fail_status_1 = qc_pcr_artifacts('temp1.region', normalbam, bedfile, bedHeader)
        #n_misprime = len([el for el in pass_fail_status_1 if el == "mispriming"])
        #n_ampbias = len([el for el in pass_fail_status_1 if el == "amplicon_bias"])
        #os.remove('temp1.region')


	##-------------Filter out variants that fail Binomial test AND QC metrics
	oldtempfile = open('Germline_All.tsv','r')
	
	#temporary file to append the p-values from binomial test to last column 
	newtempfile = open('Germline_temp.tsv','w')

	for il,line in enumerate(oldtempfile):
		if il == 0:
			newtempfile.write('\t'.join(line.strip('\n').split('\t')+['p_binom'])+'\n')
		elif il > 0:
			#position that has total_depth >= alt_depth
			if int(line.split('\t')[7]) >= int(line.split('\t')[6]):
				p_bin=stats.binom_test(int(line.split('\t')[6]),int(line.split('\t')[7]), 0.5, alternative='less')
				newline = '\t'.join(line.strip('\n').split('\t')+[str(p_bin)])+'\n'	
				newtempfile.write(newline)
			else:
				#these are positions where TVC and bamrc has mismatch
				newline = '\t'.join(line.strip('\n').split('\t')+['.'])+'\n'
				newtempfile.write(newline)

	oldtempfile.close(); newtempfile.close()

	
	#1. Keep variants that are not NO CALL, pass binomial test, min. total depth threshold, and have MAF above a certain threshold in either 1000g of exome6500
	#2. HET: keep variants that failed 1 but pass all bam QC filters
	#3. HOM: keep variants that failed 1 but pass all bam QC filters
        #4. Apply homopolymer length filter along with MAFs in 1000g and exome studies
        #NOTE: NO CALL variants are excluded right now	
	
	#command = '''awk '{if (NR==1 ||\
	#	($26>%s && $7>20 && $8>40 && $9!="./.") ||\
	#	($9=="0/1" && ($26<=%s || ($7>1 && $7<=20 && $8>=10 && $8<=40)) && ($18<=%s && $21<=%s && $22<=%s && $23>=-%s && $23<=%s && $24>=%s && $25>=%s)) ||\
	#	($9=="1/1" && ($8>=10 && $8<=40) && ($18<=%s && $24>=%s && $25>=%s))) print }'\
	#	Germline_temp.tsv | awk '{OFS="\t";$NF=""; NF--; print}' > Germline_HQ_temp.tsv''' %(
	#	alpha,
	#	alpha, p_stbias, diffq, diffmmq, diffread, diffread, readpos, dist3,
	#	p_stbias, readpos, dist3)

	#os.system(command)


        nf = open('Germline_HQ_temp.tsv','w')
        nf.write('\t'.join(header) + '\n')
        with open('Germline_temp.tsv','r') as of:
            next(of)
            for line in of:
                fields = line.strip().split('\t')
                if fields[-1] != '.' and float(fields[-1]) > alpha and float(fields[7]) >= 40 and \
                    fields[8] == './.' and \
                    ((fields[13] != '.' and float(fields[13]) > 0.2) or \
                    (fields[15] != '.' and float(fields[15]) > 0.2)):
                        nf.write('\t'.join(fields[:-1]) + '\n')
                else:
                    #HET
                    if fields[8] == '0/1' and float(fields[6]) > 10:
                        #bamrc info present
                        if fields[20] != '.':
                            if float(fields[20]) <= diffq and float(fields[21]) <= diffmmq and \
                                abs(float(fields[22])) <= diffread and float(fields[23]) >= readpos and float(fields[24]) >= dist3:
                                    nf.write('\t'.join(fields[:-1]) + '\n')
                        else:
                            #no bamrc info: will be filtered based on hplen later
                            nf.write('\t'.join(fields[:-1]) + '\n') 
                    
                    #HOM
                    elif fields[8] == '1/1' and float(fields[6]) > 10:
                        #exlen = expected_readlen(bedfile, bedHeader, fields[0], fields[1])
                        #print exlen
                        #bamrc info present
                        if fields[23] != '.':
                            if float(fields[23]) >= readpos and float(fields[24]) >= dist3:
                                nf.write('\t'.join(fields[:-1]) + '\n')
                        else:
                            #no bamrc info: will be filtered based on hplen later
                            nf.write('\t'.join(fields[:-1]) + '\n')

        nf.close()

       
        #now filter out homopolymeric FPs
	nf = open('Germline_HQ.tsv','w')
        nf.write('\t'.join(header) + '\n')
        with open('Germline_HQ_temp.tsv','r') as of:
            next(of)
            for line in of:
                fields = line.strip().split('\t')
                #snps
                if len(fields[3]) == len(fields[4]) and len(fields[3]) == 1 and fields[3] != '-' and fields[4] != '-':
                    #AAAAAaTTTTTT --> AAAAAtTTTTTT
                    if (fields[2][contextLen - 1] == fields[4] or fields[2][contextLen + 1] == fields[4]):
                        #keep those with short hplen
                        if int(fields[18]) <= hplen_snp:
                            nf.write('\t'.join(fields) + '\n')
                        else:
                            #longer hplen but high MAF in 1000g or exome studies; these are probably TPs
                            if (fields[13] != '.' and float(fields[13]) > 0.2) or (fields[15] != '.' and float(fields[15]) > 0.2):
                                nf.write('\t'.join(fields) + '\n')
                            else:
                                #discard longer hplen that don't exist in 1000g or exome studies
                                continue
                                #n_hp_snp += 1
                    else:
                        nf.write('\t'.join(fields) + '\n')
                
                #mnps
                elif len(fields[3]) == len(fields[4]) and len(fields[3]) == 2:
                    #TCTCAgAAAAA --> TCTCAaGAAAA
                    if int(fields[18]) <= hplen_mnp:
                        nf.write('\t'.join(fields) + '\n')
                    else:
                        #longer hplen but high MAF in 1000g or exome studies; these are probably TPs
                        if (fields[13] != '.' and float(fields[13]) > 0.2) or (fields[15] != '.' and float(fields[15]) > 0.2):
                            nf.write('\t'.join(fields) + '\n')
                        else:
                            #discard longer hplen that don't exist in 1000g or exome studies
                            continue
                            
                #indels
                else:
                    #keep those with short hplen
                    if int(fields[18]) <= hplen_indel:
                        nf.write('\t'.join(fields) + '\n')
                    else:
                        #longer hplen but high MAF in 1000g or exome studies; these are probably TPs
                        if (fields[13] != '.' and float(fields[13]) > 0.2) or (fields[15] != '.' and float(fields[15]) > 0.2):
                            nf.write('\t'.join(fields) + '\n')
                        else:
                            #discard longer hplen that don't exist in 1000g or exome studies
                            continue
                            #n_hp_indel += 1

        nf.close()

	nstart = int(subprocess.check_output("grep -c '' Germline_temp.tsv", shell=True)) - 1
	nend = int(subprocess.check_output("grep -c '' Germline_HQ.tsv", shell=True)) - 1

	nfilter = int(nstart) - int(nend)
	
	#Write the QC-removed variants to a file
	with open('Germline_temp.tsv') as f1:
		next(f1)
		line1 = ['\t'.join(line.strip().split('\t')[:-1]) for line in f1]
        with open('Germline_HQ.tsv') as f2:
		next(f2)
		line2 = ['\t'.join(line.strip().split('\t')) for line in f2]
	qc_failed = [x for x in line1 if x not in line2]
	with open('Germline_Removed.tsv' ,'w') as outf:
		for el in qc_failed:
			outf.write(el+'\n')
	
	os.remove('Germline_Dups.tsv');os.remove('Germline_noDups.tsv'); os.remove('Germline_temp.tsv')

	return nstart,nfilter
