#!/usr/bin/env python

import sys
import re

from gfftools import GFF
from fastatools import fasta_iter

def main():
    gff = sys.argv[1]
    fasta = sys.argv[2]

    seqs = set()

    with open(fasta,"r") as fastafile:
        for rec in fasta_iter(fastafile):
            seqs.add(rec[0])

    # clumsy way of getting base name of fasta file
    fastaname = re.sub(".fa","",re.sub(".fasta","",fasta))

    outfile = open(fastaname + ".gff","w")

    with open(gff,"r") as file_object:
        for line in file_object:
            if line[0] == "#":
                continue

            feature = GFF(line)

            if feature.seqid in seqs:
                outfile.write(line)

    outfile.close()

if __name__ == "__main__":
    main()

### EOF ###
