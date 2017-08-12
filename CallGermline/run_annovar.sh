#*****************************************
# Run annovar twice to get the completed
# annotation of all unfiltered variants
# Author: Rahul K. Das
# Date: May 4, 2016
# Last Updated: Feb 11, 2017

#*****************************************

annovar_dir="$1"
somatic_input="$2"
ncpu="$3"

$annovar_dir/table_annovar.pl \
	--thread $ncpu \
	--vcfinput \
	$somatic_input \
	$annovar_dir/humandb/ \
	--outfile annovar_variant_table \
	--buildver hg19 \
	--protocol refGene,1000g2015aug_all,snp138NonFlagged,esp6500siv2_all,cosmic70,clinvar_20160302,ljb26_all \
	--operation g,f,f,f,f,f,f \
	--nastring .  

#run annovar 2nd time to annotate mnp that are not annotated from first run
$annovar_dir/coding_change.pl annovar_variant_table.refGene.exonic_variant_function \
	$annovar_dir/humandb/hg19_refGene.txt \
	$annovar_dir/humandb/hg19_refGeneMrna.fa > annovar_variant_table_mnp.log


mkdir -p annovar_tmp
mv annovar_variant_table.hg19_multianno.txt annovar_variant_annotation_table.txt
mv annovar_variant_table.* annovar_tmp/
rm -r annovar_tmp
