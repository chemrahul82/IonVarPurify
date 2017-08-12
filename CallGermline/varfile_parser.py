"""
Parses the vcf file to extract the single-allele positions with alternative AF > 0
& Makes the vcf file of this subset of variants

"""

from __future__ import division

__Author__ =	"Rahul K. Das"
__Date__ = "May 4, 2016"
__Version__ = "2.0"
__LastModified__ = "Feb 24, 2017"


from optparse import OptionParser
import os, sys, time
import scipy.stats as stats
import statsmodels.stats.weightstats as smws
import statsmodels.stats.multitest as smm
import numpy as np


def extract_variants(vcf):

	""" Extract the variants with AF>0 & make a vcf
	@Note: if input VCF is from hotspot run, "NOCALL" positions will be included
	@params vcf: VCF with all positions
	@returns varListAll: list with all single-allele variants with AF>0; when two consecutive lines correspond
						to the same OPOS, the last one is recorded
	@returns Master.vcf: the vcf file with the extracted variants
	"""

	varListAll = []

	opos_prev = '99999999999'
	with open(vcf) as f:
		for line in f:
			if line[0] != '#':
				chrom = line.split('\t')[0]
				pos = line.split('\t')[1]

				infoDict = {}; vaf = ''
				infoFields = line.split('\t')[-2].split(':')
				infoFields1 = line.split('\t')[7].split(';')
				for ix, el in enumerate(infoFields):
					infoDict[el] = line.strip().split('\t')[-1].split(':')[ix]
				
				opos_curr = [it.split('=')[1] for it in infoFields1 if 'OPOS' in it][0]
				#excluding multi-allele sites
				#if AF field present, get the VAF
				if 'AF' in infoDict and len(infoDict['AF'].split(',')) == 1:
					vaf = infoDict['AF']
					
					#get the ref and alt depth
					if 'FRO' in infoDict and 'FAO' in infoDict:
						ref_depth = infoDict['FRO']
						alt_depth = infoDict['FAO']
					else:
						if 'RO' in infoDict and 'AO' in infoDict:
							ref_depth = infoDict['RO']
							alt_depth = infoDict['AO']
				else:
					#if AF field is not present, use FRO and FAO fields for VAF calculation
					if 'FRO' in infoDict and 'FAO' in infoDict and \
						len(infoDict['FRO'].split(',')) == 1 and len(infoDict['FAO'].split(',')) == 1:
						ref_depth = infoDict['FRO']
						alt_depth = infoDict['FAO']
						vaf = int(alt_depth)/(int(ref_depth) + int(alt_depth))
					else:
						#if AF field is not present, use RO and AO fields for VAF calculation
						if 'RO' in infoDict and 'AO' in infoDict and \
							len(infoDict['RO'].split(',')) == 1 and len(infoDict['AO'].split(',')) == 1:
							ref_depth = infoDict['RO']
							alt_depth = infoDict['AO']
							vaf = int(alt_depth)/(int(ref_depth) + int(alt_depth))
			
				if vaf != '' and float(vaf) > 0 and float(vaf) <= 1:
					if opos_prev != opos_curr:
						varListAll.append(line)
					#two variants for the same OPOS in two consecutive lines; record the second one
					else:
						varListAll.pop(-1)
						varListAll.append(line)

				opos_prev = opos_curr

				
				#vcffields = line.split('\t')[7].split(';')
				#
				#if len(vcffields[0].split('=')) > 1:
				#	 if vcffields[0].split('=')[0] == 'AF':
				#		 maf = vcffields[0].split('=')[1]
				#	 else:
				#		 #if maf not present calculate maf=AO/DP
				#		 maf = float(vcffields[0].split('=')[1])/float(vcffields[1].split('=')[1])
				#	
				#	#consider only single-allele variants
				#	 if len(vcffields[0].split('=')[1].split(',')) == 1 and float(maf)>0: 
				#		 outvcf.write(line)
				#		 varListAll.append(line)
	
	
	with open('Master.vcf', 'w') as outvcf:
		for line in varListAll:
			outvcf.write(line)

	
	return varListAll

