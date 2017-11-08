#!/usr/bin/env python

import sys

from fastatools import fasta_iter


fasta = sys.argv[1]

sys.stdout.write("seqname\tseqlength\tnumN\tpctN\n")

with open(fasta,"r") as infile:
    for rec in fasta_iter(infile):
        nam = rec[0]
        seqn = rec[1]
        seqL = len(seqn)
        numN = 0
        for char in seqn:
            if char is "N":
                numN += 1
        sys.stdout.write(nam + "\t")
        sys.stdout.write(str(seqL) + "\t" + str(numN) + "\t")
        if numN == 0:
            sys.stdout.write("0\n")
        pctN = float(numN)/float(seqL)
        sys.stdout.write(str(pctN) + "\n")

### EOF ###
