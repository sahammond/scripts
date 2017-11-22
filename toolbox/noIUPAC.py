#!/usr/bin/env python

import sys

from fastatools import fasta_iter

infile = sys.argv[1]

bueno = ["A","C","G","T","N"]

with open(infile,"r") as fasta:
    for rec in fasta_iter(fasta):
        seqn = rec[1].upper()
        nam = ">" + rec[0]
        print nam
        for base in seqn:
            if base not in bueno:
                sys.stdout.write("N")
                continue
            sys.stdout.write(base)
        sys.stdout.write("\n")
