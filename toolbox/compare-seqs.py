#!/usr/bin/env python

import sys

import argparse

from fastatools import fasta_iter


parser = argparse.ArgumentParser(
    description="Compare 2 FASTA files. Outputs tsv of OldID and NewID. Collisions while loading the OldIDs are written to collisions.txt")
parser.add_argument('OldFASTA',action='store',help="FASTA of Old or reference sequences")
parser.add_argument('NewFASTA',action='store',help="FASTA of New or query sequences",nargs="?")

args = parser.parse_args()

oldsqn = args.OldFASTA
if args.NewFASTA:
    newsqn = args.NewFASTA

# read in old scaffs with sqn as key and id as value
# then read in new scaffs and compare sqn to keys and write out old and new ids upon match
# or, if only one FASTA given, output unique sequences

seqdict = {}

collisions = open("collisions.txt","w")

with open(oldsqn,"r") as infile:
    for rec in fasta_iter(infile):
        seqid = rec[0]
        sqn = rec[1]
        if sqn in seqdict:
            print >> collisions, "Collision between " + seqid + " AND previously added " + seqdict[sqn]
        else:
            seqdict[sqn] = seqid
collisions.close()

if args.NewFASTA:
    print "OldID\tNewID"

    with open(newsqn,"r") as infile:
        for rec in fasta_iter(infile):
            seqid = rec[0]
            sqn = rec[1]
            if sqn in seqdict:
                print seqdict[sqn] + "\t" + seqid
else:
    # write out unique sequences
    for sqn, sid in seqdict.items():
        print ">" + sid
        print sqn


### EOF ###
