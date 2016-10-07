#!/usr/bin/env python

"""Select desired genes from a gff based on transcript IDs

If a given gene has even a single transcript in the to-select list,
then all isoforms will be collected

Usage: select-gene-records.py list-of-transcript-ids.txt gene-records.gff
Output: gene-records-selected.gff

"""

import sys
import re

import gfftools

def main():
    idlist = sys.argv[1]
    gff = sys.argv[2]
    outnam = re.sub(".gff","-selected.gff",gff)

    ids = set()
    with open(idlist,"r") as file_object:
        for line in file_object:
            ids.add(line.strip("\n"))
            deriv_geneID = re.sub("-mRNA-[0-9]","",line.strip("\n"))
            ids.add(deriv_geneID)

    outfile = open(outnam,"w")

    with open(gff,"r") as file_object:
        for line in file_object:
            if line[0] == "#":
                outfile.write(line)
                continue
            feature = gfftools.GFF(line)
            if feature.id in ids:
                outfile.write(line)
            else:
                for papa in feature.parent:
                    if papa in ids:
                        outfile.write(line)
                        break
    outfile.close()

if __name__ == "__main__":
    main()

### EOF ###
