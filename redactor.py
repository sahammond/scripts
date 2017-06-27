#!/usr/bin/env python

import sys
import re


def redact(scaffold, min_length):
    """'redact' bases in scaftigs less than 50 bp long by replacing with N"""

    scaftigs = scaffold.strip().upper().split("N")
    scaftigs = [x + "N" for x in scaftigs]
    scaftigs[-1] = scaftigs[-1][:-1] # correct extra terminal N
    outbuff = []

    for tig in scaftigs:
        # allow for the extra N
        if len(tig) < (min_length + 1):
            buff = []
            for nuc in tig:
                buff.append("N")
            outbuff.append("".join(buff))
        else:
            outbuff.append(tig)

    to_return = "".join(outbuff)
    return to_return
    

MINL = 50 # NCBI minimum scaftig length
seqfile = sys.argv[1] # fasta file of scaffolds

with open(seqfile,"r") as file_object:
    for line in file_object:
        if line[0] == ">":
            sys.stdout.write(line)
        else:
            print redact(line, MINL)
