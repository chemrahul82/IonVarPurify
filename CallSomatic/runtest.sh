python IonSomaticDecoder.py -s /rawdata/Rahul/Development/Somatic_Calling_Pipeline/Dev/version2.0/CallSomatic\
	-c /rawdata/software/bam-readcount/bin/bam-readcount\
	-A /rawdata/software/annovar_Feb2016/\
	-B /rawdata/support_files/BED/AmpliSeq-Exome.bed\
	-f /results/referenceLibrary/tmap-f3/hg19/hg19.fasta\
	-T /mnt/Despina/projects/PNET/A_146/QC/allA_146_Normal_Merged_11262014vsA_146_Tumor_Merged_11132014/VCF2_Final.vcf\
	-N /mnt/Despina/projects/PNET/A_146/QC/allA_146_Normal_Merged_11262014vsA_146_Tumor_Merged_11132014/VCF1_Final.vcf\
	-b /mnt/Despina/projects/PNET/A_146/Tumor/Merged/PTRIM.bam\
	-o /rawdata/Rahul/Development/Somatic_Calling_Pipeline/Dev/version2.0/CallSomatic/testing3/\
	-Y no\
	-Z paired\
	-E hotspot\
	-P 16\
	-W 0\
	--alpha1 0.05


#python IonSomaticDecoder.py -s /rawdata/Rahul/Development/Somatic_Calling_Pipeline/Dev/version2.0/CallSomatic\
#	-c /rawdata/software/bam-readcount/bin/bam-readcount\
#	-A /rawdata/software/annovar_Feb2016/\
#	-B /rawdata/support_files/BED/Radiogenomics/IAD38165_Designed_muc16znf717excluded.bed\
#	-f /results/referenceLibrary/tmap-f3/hg19/hg19.fasta\
#	-k /rawdata/Rahul/Analysis/LBR_Analysis/Somatic/PON_BlackList/BlackList_Positions.txt\
#	-T /mnt/Orcus/projects/LungBio/pair_035/QC/allpair_035_Normal_Merged_07172015vspair_035_Tumor_Merged_07202015/VCF2_Final.vcf\
#	-N /mnt/Orcus/projects/LungBio/pair_035/QC/allpair_035_Normal_Merged_07172015vspair_035_Tumor_Merged_07202015/VCF1_Final.vcf\
#	-b /mnt/Orcus/projects/LungBio/pair_035/Tumor/Merged/PTRIM.bam\
#	-o /rawdata/Rahul/Development/Somatic_Calling_Pipeline/Dev/version2.0/CallSomatic/testing1/\
#	-Z paired\
#	-l 4 \
#	--sample_type ffpe \
#	-p 0.8\
#	-Y no\
#	-H 5\
#	-W 25\
#	-E hotspot

#-V /mnt/Orcus/projects/LungBio/pair_035/QC/allpair_035_Normal_Merged_07172015vspair_035_Tumor_Merged_07202015/VCF1_Final.vcf\

