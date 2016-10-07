#!/usr/bin/env python

import sys
import re


def main():

    """increase-confidence-tbl.py

    purpose: ncbi disallows annotations that are all 'similar to' something.
    Here, I definitively identify the annotations with genevalidator scores >= 90
    usage: increase-confidence-tbl.py annotations.tbl annotations-locus_tag-conversion-key.tsv list-of-confident-annotations.txt > annotations-wConfidence.gff

    """

    tbl = sys.argv[1]
    converto = sys.argv[2]
    confident = sys.argv[3]
    outname = re.sub(".tbl","-wConfidence.tbl",tbl)

    transl = {}
    best = set()
    
    with open(converto,"r") as convo:
        for line in convo:
            rec = line.strip("\n").split("\t")
            if re.findall("mRNA",rec[0]):
                transl[rec[1]] = rec[0]
    
    with open(confident,"r") as conf:
        for line in conf:
            best.add(line.strip("\n"))

    outfile = open(outname,"w")

    # first pass through file to get names
    protNames = {}
    with open(tbl,"r") as annots:
        goodAnnot = False
        for line in annots:
            # only protein coding genes are currently "annotated"
            # so no need to modify ncRNA records
            if re.findall("product",line):
                locus = line.strip("\n").split("\t")[-1]
                if locus in transl:
                    orig = transl[locus]
                    if orig in best:
                        goodAnnot = True
            elif re.findall("prot_desc",line) and goodAnnot:
                # if the annotation is good, then no prot_desc required
                fixedLine = re.sub("similar to ","",line)
                fixedLine = re.sub("prot_desc","product",fixedLine)
                protNames[locus] = fixedLine
                goodAnnot = False
                continue

    with open(tbl,"r") as annots:
        goodAnnot = False
        for line in annots:
            # only protein coding genes are currently "annotated"
            # so no need to modify ncRNA records
            if re.findall("product",line):
                locus = line.strip("\n").split("\t")[-1]
                if locus in protNames:
                    outfile.write(protNames[locus])
                    goodAnnot = True
                else:
                    outfile.write(line)
            elif re.findall("prot_desc",line) and goodAnnot:
                if goodAnnot:
                    continue
                else:
                    outfile.write(line)
            else:
                outfile.write(line)
                goodAnnot = False
    outfile.close()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print main.__doc__
        sys.exit(1)
    main()

### EOF ###
