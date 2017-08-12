import pysam
from collections import Counter



#C->G; LBR_35_T
#for pileupcolumn in samfile.pileup("chr1", 27107026, 27107027):
#A->C; LBR_35_T
#for pileupcolumn in samfile.pileup("chr6", 117710520, 117710521):
#C->A: LBR_35_T; KRAS
#for pileupcolumn in samfile.pileup("chr12", 25398284, 25398285):
	

def varStartPos(chrom,start,end,ref,var):
	nvar = 0; nref = 0
	varReadStart = []; refReadStart = []

	for pileupcolumn in samfile.pileup(chrom, start, end):
	
		#print pileupcolumn.pos
		#print ("\ncoverage at base %s = %s" %
		#(pileupcolumn.pos, pileupcolumn.n))
		#print pileupcolumn.nsegments
		for pileupread in pileupcolumn.pileups:
			#print pileupread.alignment.cigartuples
			#print pileupread.alignment.cigarstring
			if not pileupread.is_del:
				if pileupcolumn.pos == start:
			
					#print pileupread.alignment.query_sequence[pileupread.query_position]
					#print pileupread.alignment.query_name
					#print pileupread.alignment.query_sequence
					#print ('\tbase in read %s = %s' %
					#	(pileupread.alignment.query_name,
					#	pileupread.alignment.query_sequence[pileupread.query_position]))
					
					if pileupread.alignment.query_sequence[pileupread.query_position] == var:
						#print pileupread.alignment.pos,pileupread.alignment.cigarstring
						#print "var"
						#print pileupread.alignment.pos,pileupread.alignment.qstart,pileupread.alignment.qend,pileupread.alignment.qlen,pileupread.alignment.rlen,pileupread.alignment.tlen
						#print "var"
						#print pileupread.alignment.pos,pileupread.alignment.cigarstring
						nvar += 1
						varReadStart.append(pileupread.alignment.pos)
					
					if pileupread.alignment.query_sequence[pileupread.query_position] == ref:
						#print "ref"
						#print pileupread.alignment.pos,pileupread.alignment.cigarstring
						nref += 1
						refReadStart.append(pileupread.alignment.pos)
						#print pileupread.alignment.pos,pileupread.alignment.qstart,pileupread.alignment.qend,pileupread.alignment.qlen,pileupread.alignment.rlen,pileupread.alignment.tlen

	#print len(Counter(refReadStart).most_common()), len(Counter(varReadStart).most_common())
	#print nref,len(Counter(refReadStart).most_common()),Counter(refReadStart).most_common(), nvar, len(Counter(varReadStart).most_common()), Counter(varReadStart).most_common()
	
	
	
	#return dict(Counter(refReadStart)), dict(Counter(varReadStart))
	return Counter(refReadStart).most_common(), Counter(varReadStart).most_common()

			#print base
		#if not pileupread.is_del and not pileupread.is_refskip:
		#	# query position is None if is_del or is_refskip is set.
		#	print ('\tbase in read %s = %s' %
		#		(pileupread.alignment.query_name,
		#		pileupread.alignment.query_sequence[pileupread.query_position]))


if __name__ == "__main__":
	
	samfile = pysam.AlignmentFile("/mnt/Orcus/projects/LungBio/pair_035/Tumor/Merged/PTRIM.bam", "rb")
	
	for pileupcolumn in samfile.pileup("chr1", 16458630, 16458630):
		print pileupcolumn.pos
	
	##samfile = pysam.AlignmentFile("/mnt/Orcus/projects/LungBio/pair_123/Tumor/T-2/PTRIM.bam", "rb")
	#LBR_049
	#samfile = pysam.AlignmentFile("/mnt/Orcus/projects/Radiogenomics/pair_008/Tumor/Merged/PTRIM.bam", "rb")
	
	#LBR_082
	#samfile = pysam.AlignmentFile("/mnt/Orcus/projects/Radiogenomics/pair_010/Tumor/Merged/PTRIM.bam", "rb")

	#Wales_Case1A02
	#samfile = pysam.AlignmentFile("/data/results/projects/Wales/case1_A02/Run1/tvc_rerun/PTRIM.bam", "rb")

	
	with open('/scratch/support_files/BED/IAD38165_Designed_muc16znf717excluded.bed') as f:
		ampStart = [line.split('\t')[0]+'_'+line.split('\t')[1] for line in f]

	with open('LBR35_region', 'r') as f:
		lines = [line.strip().split('\t') for line in f]

	for line in lines:
		#print line
		chrom = line[0]
		start = int(line[1])
		end = int(line[2])
		ref = line[3]
		var = line[4]

		if start == end and len(ref) == 1 and len(var) == 1 and ref != '-' and var != '-':
			#print chrom, start
			d1,d2 = varStartPos(chrom,start-1,start,ref,var)
			#print d1, d2
			if chrom+'_'+str(d2[0][0]) not in ampStart:
				print "mispriming!"
			if len(d2) == 1:
				if len(d1) > 1:
					#if d2.keys()[0] in d1:
					if d2[0][0] in [el[0] for el in d1]:
						print "KEEP"
					else:
						print "REJECT"


	samfile.close()
