__Author__ =    "Rahul K. Das"
__Version__ = "2.0"
__LastModified__ = "Feb 11, 2017"


import os, sys, subprocess
from itertools import groupby

def cal_hplen(sequence, scanLength, ref, alt):
    
    """
    @param sequence: input sequence
    @param scanLength: number of bases to scan each side around the reference position
    @param ref: reference allele
    @param alt: alternative alele
    @returns hrun_adj: adjusted homopolymer length

    """

    #----------------------compute the homopolymer length statistics------------------
    
    #get the different homopolymeric runs, their lengths, and component base
    grps = [list(g) for k,g in groupby(sequence)]
    grps_base = [list(k) for k,g in groupby(sequence)]
    grplen = [len(i) for i in grps]
   

    #----------------------------------------------#
    #               snp and 1-bp indels            #
    #----------------------------------------------#

    if len(ref) == 1 and len(alt) == 1:
        #-----------------snps--------------------
        if ref != '-' and alt != '-':
            query_pos = scanLength + 1
            nsum = 0
            for i,l in enumerate(grplen):
                nsum += l
                if nsum >= query_pos:
                    bidx = i
                    break
           
            #default HRUN value
            hrun = str(grplen[bidx])
            
            if int(hrun) == 1:
                #simplest snp
                hrun_adj = hrun
            
            #    if grps_base[bidx + 1][0] != alt and grps_base[bidx - 1][0] != alt:
            #        hrun_adj = hrun
            #    else:
            #        #<----signature-1---->
            #        #TCATCtCCAGC 
            #        #TCATCcCCAGC
            #        if grps_base[bidx + 1][0] == alt:
            #            hrun_right = grplen[bidx + 1]
            #        else:
            #            hrun_right = 0
            #        if grps_base[bidx - 1][0] == alt:
            #            hrun_left = grplen[bidx - 1]
            #        else:
            #            hrun_left = 0
            #
            #        hrun_adj = str(int(hrun) + hrun_left + hrun_right)
            
            if int(hrun) > 1:
                #<----signature-2---->
                #TTTTtAAAAAA or TTTTTaAAAAA
                #TTTTaAAAAAA    TTTTTtAAAAA
                #hprun of ref allele >1 & ref-1 position has alt allele with hprun >1

                if sequence[query_pos] == ref and sequence[query_pos - 2] == alt and grplen[bidx - 1] > 1:
                    hrun_adj = str(int(hrun) + grplen[bidx - 1])
                elif sequence[query_pos - 2] == ref and sequence[query_pos] == alt and  grplen[bidx + 1] > 1:
                    hrun_adj = str(int(hrun) + grplen[bidx + 1])
                else:
                    hrun_adj = hrun

            #else:
            #    hrun_adj = hrun

       
        #-------------deletion---------------
        if alt == '-':
            query_pos = scanLength + 1
            nsum = 0
            for i,l in enumerate(grplen):
                nsum += l
                if nsum >= query_pos:
                    bidx = i
                    break
           
            hrun = str(grplen[bidx])
            hrun_adj = hrun


        #---------------insertion-------------
        if ref == '-':
            query_pos = scanLength + 1
            nsum = 0
            for i,l in enumerate(grplen):
                nsum += l
                if nsum >= query_pos:
                    bidx = i
                    break


            #set the hrun based on whether ref or ref+1 match the inserted base 
            if grps_base[bidx][0] == alt:
                hrun = str(grplen[bidx])
            else:
                if grps_base[bidx + 1][0] == alt:
                    hrun = str(grplen[bidx + 1])
                else:
                    hrun = '0'

            hrun_adj = hrun

        
    
    #----------------------------------------------#
    #                    mnps                      #
    #----------------------------------------------#
    elif len(ref) > 2 and len(ref) == len(alt):
        #query only ref position
        query_pos = scanLength + 1

        nsum = 0
        for i,l in enumerate(grplen):
            nsum += l
            if (nsum >= query_pos):
                bidx = i
                break

        if grps_base[bidx][0] == ref[0]:
            hrun = str(grplen[bidx])
        #else:
        #    hrun = '0'

        hrun_adj = hrun

    #<-----Signature-3: diblock mnps---------->
    #examples: ATATAtAAAAA {tA-->AT}; CTCCTcTTTTT {cT-->TC}

    elif len(ref) == 2 and len(ref) == len(alt):
        ref_left = ref[0]
        ref_right = ref[-1]

        if ref_left == ref_right:
            hrun_adj = '2'
        else:
            query_pos_left = scanLength + 1
            query_pos_right = scanLength + 2
        
            #get left base index
            nsum = 0
            for i,l in enumerate(grplen):
                nsum += l
                if (nsum >= query_pos_left):
                    bidx_left = i
                    break
        
            #get right base index
            nsum = 0
            for i,l in enumerate(grplen):
                nsum += l
                if (nsum >= query_pos_right):
                    bidx_right = i
                    break

            #hprun w.r.t. the 1st ref. position
            hrun1 = grplen[bidx_left]

            #but consider hprun w.r.t the 2nd ref position too
            hrun2 = grplen[bidx_right]

            if hrun1 > 1 and hrun2 > 1:
                hrun_adj = str(hrun1 + hrun2)
            elif hrun1 == 1 and hrun2 == 1:
                hrun_adj = str(hrun1)
            else:
                hrun_adj = str(max([hrun1,hrun2]))
        
        #vlen = len(ref)
        #hrun_list = []
        #query_pos = range(scanLength + 1, scanLength + vlen + 2)
        #for query in query_pos:
        #    nsum = 0
        #    for i,l in enumerate(grplen):
        #        nsum += l
        #        if nsum >= query:
        #            bidx = i
        #            break

        #    hrun_list.append(grplen[bidx])
        #hrun = str(max(hrun_list))
    

    
    #----------------------------------------------# 
    #           deletions of length > 1            #
    #     examples: 'TA' --> '-' or 'TAT' --> 'A'  #
    #----------------------------------------------#
    elif len(ref) > 1 and len(ref) > len(alt):
        vlen = len(ref) 
        ref_left = ref[0]
        ref_right = ref[-1]

        query_pos_left = scanLength + 1
        query_pos_right = scanLength + vlen

        #get left base index
        nsum = 0
        for i,l in enumerate(grplen):
            nsum += l
            if (nsum >= query_pos_left):
                bidx_left = i
                break
        
        #get right base index
        nsum = 0
        for i,l in enumerate(grplen):
            nsum += l
            if (nsum >= query_pos_right):
                bidx_right = i
                break

        #hprun w.r.t. the ref. position; this is what is in VCF
        if grps_base[bidx_left][0] == ref_left:
            hrun1 = grplen[bidx_left]
        else:
            hrun1 = 0

        #but consider hprun w.r.t the right most base of the region where deletion occurs
        if grps_base[bidx_right][0] == ref_right:
            hrun2 = grplen[bidx_right]
        else:
            hrun2 = 0

        #hrun = str(hrun1)
        if bidx_left == bidx_right:#AA --> '-'
            hrun_adj = str(hrun1)
        else:
            hrun_adj = str(hrun1 + hrun2)



    #----------------------------------------------# 
    #              insertion of len >1             #
    #examples: '-' --> 'TA' or 'A' --> 'TAT'       #
    #        or '-' --> AA                         #
    #----------------------------------------------#
    elif len(alt) > len(ref):
        hrun_list = [];bidx_list = []
        vlen = len(alt)

        # if of 2nd type, no. of non-terminal bps
        paddlen = vlen - 2 

        #terminal bps
        alt_left = alt[0]
        alt_right = alt[-1]

        #get hplen of lhs
        query_pos_left = scanLength + 1 - paddlen
        query_pos_right = scanLength + 2
        #query_pos = range(scanLength + 2, scanLength + 4)

        #get left base index
        nsum = 0
        for i,l in enumerate(grplen):
            nsum += l
            if (nsum >= query_pos_left):
                bidx_left = i
                break
        
        #get right base index
        nsum = 0
        for i,l in enumerate(grplen):
            nsum += l
            if (nsum >= query_pos_right):
                bidx_right = i
                break


        if grps_base[bidx_left][0] == alt_left:
            hrun1 = grplen[bidx_left]
        else:
            hrun1 = 0

        if grps_base[bidx_right][0] == alt_right:
            hrun2 = grplen[bidx_right]
        else:
            hrun2 = 0


        if bidx_left == bidx_right: # '-' --> AA
            hrun_adj = str(hrun1)
        else:
            hrun_adj = str(hrun1 + hrun2)   
        
        
        #else:
        #    vlen = len(ref)
        #    stop_pos = scanLength + vlen
        
        #for query in query_pos:
        #    nsum = 0
        #    for i,l in enumerate(grplen):
        #        nsum += l
        #        if nsum >= query:
        #            bidx = i
        #            break
        #    hrun_list.append(grplen[bidx])
        #    bidx_list.append(bidx)

        #hrun = str(max(hrun_list))

        ##if 

    #else:
    #    hrun = '.'

   
    return hrun_adj

    #longest homopolymer length within a sequence context of given size around the reference position
    #hrun = str(max(sum(1 for i in g) for k,g in groupby(motifList2[idx])))
