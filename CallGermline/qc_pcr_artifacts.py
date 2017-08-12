"""

#Accounting for false positives due to mispriming events or amplicon bias in PCR
#determine the starting mapped positions of reads & group the reference & variant supporting reads

"""

__Author__ =    "Rahul K. Das"
__LastModified_ = "Nov 08, 2016"


import pysam
from collections import Counter
import time

def groupReads(bamFileObj, chrom, start, end, ref, var):
    
    """
    @param: bamFileObj: pysam alignment object
    @param: chr: chromosome no. in the same format as in bam
    @param: start: variant start position
    @param: ref: reference allele
    @param: var: variant allele
    @return: two lists of tuples, one for reference & one for variant reads,
            each tuple contains start position & counts of a read group
    """
  
    varReadStart = []; refReadStart = []
    #group variant and reference reads by start mapped positions:
    for pileupcolumn in bamFileObj.pileup(chrom, start, end):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and pileupcolumn.pos == start:
                if pileupread.alignment.query_sequence[pileupread.query_position] == var:
                    varReadStart.append(pileupread.alignment.pos)
                if pileupread.alignment.query_sequence[pileupread.query_position] == ref:
                    refReadStart.append(pileupread.alignment.pos)

    return Counter(refReadStart).most_common(), Counter(varReadStart).most_common()


def qc_pcr_artifacts(varBed, bam, projectBed, bedHeader):
    
    """
    @param varBed: the list of variants in bed format but +1-based positions
    @param bam: the input bam
    @param projectBed: the project bed files with amplicons specs
    @param bedHeader: header in project bed??
    @return pass_fail_status: list with "pass","amplicon_bias","mispriming" status of variants
    """
    
    #----------------------------------------------
    #--------------amplicon information------------
    #----------------------------------------------

    #chromosomes from project bed
    with open(projectBed, 'r') as f:
        if bedHeader=='yes':
            chroms = [key for key in dict(Counter([line.split('\t')[0] for il, line in enumerate(f) if il>0])).keys()]
        else:
            chroms = [key for key in dict(Counter([line.split('\t')[0] for il, line in enumerate(f)])).keys()]
        #print chroms

    #dictionary: keys are chrom and values are lists with start/end positions
    ampStarts = {}; ampEnds = {}
    for chrom in chroms:
        ampStarts[chrom] = []
        ampEnds[chrom] = []
    
    with open(projectBed, 'r') as f:
        if bedHeader=='yes':
            next(f)    
        for line in f:
            ampStarts[line.split('\t')[0]].append(int(line.split('\t')[1]))

    with open(projectBed, 'r') as f:
        if bedHeader=='yes':
            next(f)
        for line in f:
            ampEnds[line.split('\t')[0]].append(int(line.split('\t')[2]))


    #----------------------------------------------
    #---------evaluate reads for PCR errors--------
    #----------------------------------------------

    #pysam alignment object
    bamFileObj = pysam.AlignmentFile(bam, 'rb')

    #region/bed file for variants
    with open(varBed, 'r') as f:
        lines = [line.strip().split('\t') for line in f]
    

    pass_fail_status = []
    for line in lines:
        chrom = line[0]
        start = int(line[1])
        ref = line[2]
        var = line[3]

        #only for SNPs
        if len(ref) == 1 and len(var) == 1 and ref != '-' and var != '-':    
            refGrp, varGrp = groupReads(bamFileObj,chrom,start-1,start,ref,var) # 0-based indexing
            #print line
            #print refGrp, varGrp
            
            #--------Mispriming-------------
            # check if variant read start positions overlap with amplicon start position
            # and how many groups
            nVarReadGrp = len(varGrp)
            nOverlap = 0
            varReadStartPos = [int(el[0]) for el in varGrp]
            for startpos in varReadStartPos:
                if startpos in ampStarts[chrom]:
                    nOverlap += 1
                
            if nOverlap == 0:
                pass_fail_status.append("mispriming")

            
            #-------Amplicon bias-------------
            
            nVarAmp = 0
            #no. of amplicons spanning the chromosome
            namp = len(ampStarts.get(chrom))
            for ia in range(0,namp):
                #no. of amplicons that span the variant
                if (int(start) >= int(ampStarts.get(chrom)[ia])) and (int(start) <= int(ampEnds.get(chrom)[ia])):
                    nVarAmp += 1
                # if amplicon start is greater than variant pos, go to next variant
                if int(start) < int(ampStarts.get(chrom)[ia]):
                    break
            
            #one match out of >1 amplicons
            if (nOverlap == 1):
                if (nVarAmp > 1):
                    pass_fail_status.append("amplicon_bias")
                elif (nVarAmp == 1):
                    pass_fail_status.append('pass')

            if (nOverlap > 1):
                pass_fail_status.append('pass')

            #print nOverlap, nVarReadGrp, nVarAmp
        else:
            pass_fail_status.append('.')
   
    bamFileObj.close()

    return pass_fail_status
