#!/usr/bin/env python

import sys
import re

"""
purpose: ncbi disallows annotations that are all 'similar to' something.
    Here, I definitively identify the annotations with genevalidator scores >= 90
usage: increase-confidence.py annotations.gff maker-to-ncbi-conversion-table.tsv list-of-confident-annotations.txt > annotations-wConfidence.gff
"""

def main():
    gff=sys.argv[1]
    converto=sys.argv[2]
    confident=sys.argv[3]
    
    transl={}
    best=set()
    
    with open(converto,"r") as convo:
        for line in convo:
            rec=line.strip("\n").split("\t")
            if rec[0] == 'mRNA':
                transl[rec[2]]=rec[1]
    
    with open(confident,"r") as conf:
        for line in conf:
            best.add(line.strip("\n"))
    
    with open(gff,"r") as annots:
        for line in annots:
            rec=line.strip("\n").split("\t")
            if rec[2] == "mRNA":
                reco=rec[8].split(";")[0][3:]
                ido=transl[reco]
                if ido in best:
                    sys.stdout.write(re.sub('similar to ','',line))
                else:
                    sys.stdout.write(line)
            else:
                sys.stdout.write(line)

if __name__ == "__main__":
    main()

### EOF ###
