python -u QC_Germline.py \
	-s /rawdata/Rahul/Development/Somatic_Calling_Pipeline/Dev/version2.0/CallGermline/\
	-N /rawdata/software/annovar_Feb2016/\
	-c /rawdata/software/bam-readcount/bin/bam-readcount\
	-f /results/referenceLibrary/tmap-f3/hg19/hg19.fasta\
	-b /mnt/Orcus/projects/LungSquamous/LS_204/Normal/Merged/PTRIM.bam\
	-v /mnt/Orcus/projects/LungSquamous/LS_204/QC/allLS_204_Normal_Merged_01182017vsLS_204_Tumor_Merged_01192017/VCF1_Final.vcf\
	-o testing/ -P 24\
	-B /rawdata/support_files/BED/AmpliSeq-Exome.bed\
	-Y no
