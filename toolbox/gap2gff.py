#!/usr/bin/env python

# Collect runs of N in scaffolds and write out positions as gff

import argparse
from itertools import groupby

from fastatools import fasta_iter as f

parser = argparse.ArgumentParser(description="Collect runs of N in scaffolds and write out positions as gff.")
parser.add_argument('scaffolds', action='store', help='The fasta file of scaffolds')
parser.add_argument('--outfile', '-o', action='store', help='Name to use for output file [gaps.gff]', default='gaps.gff')

args = parser.parse_args()

## gap2gff
#
#   take a sequence and a position
#   output a gff line
def gap2gff(scaffold, startpos, gap_sequence, label=''):
    endpos = int(startpos) - len(gap_sequence)
    if not label:
        label = "gap" + str(endpos)
    gfflist = [scaffold, "gap2gff", "identity", str(endpos), str(startpos), ".", "+", ".", "ID=" + label]
    gffline = "\t".join(gfflist)

    return gffline


## getgaps
#
#    split input string at N
#    return list of gaps
def getgaps(seqname, seq):
    #for base in string add N to a growing string, break if base not N
    gaps = []
    thisgap = ''
    pos = 1 # gff is 1-based
    ingap = 'F'
    for base in seq:
        if base != 'N':
            if ingap == 'T':
                gapline = gap2gff(seqname, pos, thisgap)
                gaps.append(gapline)
#                gaps.append(thisgap)
                thisgap = ''
                ingap = 'F'
        elif ingap == 'T':
            thisgap = thisgap + base
        elif base == 'N':
            ingap = 'T'
            thisgap = thisgap + base
        pos += 1
    else:
        if len(thisgap) > 0:
            gaps.append(gapline)
#            gaps.append(thisgap)
    return gaps

## main
#
scaffolds = []
with open(args.scaffolds,'r') as fasta:
    for scaf in f(fasta):
        scaffolds.append(scaf)

with open(args.outfile,"w") as outfile:
    for scaf in scaffolds:
        scafgaps = getgaps(scaf[0], scaf[1])
        for i in scafgaps:
            outfile.write(i + "\n")

### EOF ###
