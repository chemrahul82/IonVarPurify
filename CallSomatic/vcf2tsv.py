"""
convert vcf to tsv with relevant fields
"""

from __future__ import division

__Author__ =    "Rahul K. Das"
__Version__ = "2.0"
__LastUpdated__ = "Feb 24, 2017"


import os, sys


def vcf2tsv(vcf_file, out_tsv):

    """
    @param: vcf_file: input VCF
    @return: out_tsv: output tsv file
    """

    outfile = open(out_tsv, 'w')
    header = '\t'.join(['Chr','Pos','Ref','Alt','GT','VAF','Alt_Depth','Ref_Depth']) + '\n'
    outfile.write(header)

    #parse the VCF file
    vf = open(vcf_file, 'r')
    for line in vf:
        if line[0] != "#":
            chrom = line.split('\t')[0]
            pos = line.split('\t')[1]
            ref = line.split('\t')[3]
            alt = line.split('\t')[4]

            infoDict = {}; vaf = ''
            
            infoFields = line.split('\t')[-2].split(':')
            for ix, el in enumerate(infoFields):
                infoDict[el] = line.strip().split('\t')[-1].split(':')[ix]

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
                        
            #genotypes
            if vaf != '' and float(vaf) <= 1:
                if infoDict['GT'] == '0/0':
                    gt = 'WT'
                elif infoDict['GT'] == '0/1':
                    gt = 'HET'
                elif infoDict['GT'] == '1/1':
                    gt = 'HOM'
                else:
                    gt = '.'

                outfile.write('\t'.join(map(str,[chrom,pos,ref,alt,gt,vaf,alt_depth,ref_depth])) + '\n')

    outfile.close()
