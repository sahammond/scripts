#!/usr/bin/env python

# usage: gff3-to-tbl.py genes.gff genome.fa alignment.blastp first_locus_number

import sys
import re
import copy

from fastatools import fasta_iter as f
import gfftools

############
# parse arguments #TODO replace with argparse

gff = sys.argv[1]
fasta = sys.argv[2]
aln = sys.argv[3]
loc = sys.argv[4]

############
# declarations

# NCBI-assigned locus_tag
LOCUS="AB205"
# start number for locus id iteration
#LOCS=50
LOCS=int(loc)
# width for locus ids (i.e. total with padding)
LOCW=7
# jump size between loci (good idea to allow for genes identified between current ones)
LOCJ=10
# miniumum percent identity to accept alignment for annotation
MIN_IDENT=25
# minimum percent coverage to accept alignment for annotation
MIN_COV=50
# string to split reference protein names with. Script expects swissprot-style protein
#  names in reference.fa, use 'OS='
BRK='OS=' 
# field sep for product line
SEPO=' '
# string to prefix description with
PRE='similar to'
# label to give queries that have no good annotation
UNK='hypothetical protein'
# alphabet to use for isoform labeling
ISOALPHA=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S',
            'T','U','V','W','X','Y','Z']

#gffout=open(re.sub(".gff","-prepared.gff",gff),"w")
tblout=open(re.sub(".gff","-prepared.tbl",gff),"w")
lokey=open(re.sub(".gff","-locus_tag-conversion-key.tsv",gff),"w")

###################
# MAIN

genes = {}
# read in the entries from the gff file
with open(gff,"r") as infile:
    print "Reading the annotations"
    for line in infile:
        if line[0] == "#" or line[0] == "-":
            continue

        feature = gfftools.GFF(line)

        if feature.type == "gene":
            genes[feature.id] = gfftools.Gene(feature)

        elif feature.type == "mRNA" or feature.type == "ncRNA":
            # transcripts can only have one parent, so access it explicitly
            genes[feature.parent[0]].add_transcript(feature)

        elif feature.type == "exon":
            if feature.name:
                genes[feature.name].transcript[feature.parent[0]].add_exon(feature)
                continue
            for parent in feature.parent:
                for gene in genes:
                    if re.findall(genes[gene].id,parent):
                        genes[gene].transcript[parent].add_exon(feature)

        elif feature.type == "CDS":
            if feature.name:
                genes[feature.name].transcript[feature.parent[0]].add_exon(feature)
                continue
            for parent in feature.parent:
                for gene in genes:
                    if re.findall(genes[gene].id,parent):
                        genes[gene].transcript[parent].add_cds(feature)

###################
# adjust small intron boundaries (ncbi min intron length is 10 bp)
# use steps of 3 to preserve frame
#  this will effectively delete up to 4 amino acids per adjustment (max if intron len was 1)
print "Checking intron lengths"
for entry in genes:
    for rec in genes[entry].transcript:
#        print "Checking intron lengths in " + rec
        for j in range(0,len(genes[entry].transcript[rec].exons)):
            dist = 11
            try:
                dist = int(genes[entry].transcript[rec].exons[j+1].start) - int(genes[entry].transcript[rec].exons[j].end)
#                print dist
            except:
                pass
            while dist < 11: # if use 'while dist < 10:' then too-small introns remain...
                print "Correcting intron length between " + rec + " exons " + str(j+1) + " and " + str(j+2)
                genes[entry].transcript[rec].exons[j+1].start = str(int(
                            genes[entry].transcript[rec].exons[j+1].start) + 3)
                dist += 3
        # adjust the corresponding CDS start too, if the CDS start == that exon's start!
        for j in range(0,len(genes[entry].transcript[rec].cds)):
            dist = 11
            try:
                dist = int(genes[entry].transcript[rec].cds[j+1].start) - int(genes[entry].transcript[rec].cds[j].end)
#                print dist
            except:
                pass
            while dist < 11: # if use 'while dist < 10:' then too-small introns remain...
                print "Applying intron length correction between " + rec + " exons " + str(j+1) + " and " + str(j+2) + " to the CDS too"
                genes[entry].transcript[rec].cds[j+1].start = str(int(
                            genes[entry].transcript[rec].cds[j+1].start) + 3)
                dist += 3

###################
# Collect annotations for each transcript
aseen = set()
annots = {}
with open(aln,"r") as blast:
    print "Loading the annotations"
    for rec in blast:
        thisaln=rec.split("\t")
        if thisaln[0] not in aseen:
            cov=thisaln[-1]
            if thisaln[2] >= MIN_IDENT and cov >= MIN_COV:
                thisannot=thisaln[-2].split(" ")
                refnam=''
                for i in range(0,len(thisannot)):
                    if i == 0 and len(re.findall(BRK,thisannot[i])) < 1:
                        refnam+=thisannot[i]
                    elif i > 0 and len(re.findall(BRK,thisannot[i])) < 1:
                        refnam+=SEPO+thisannot[i]
                    else:
                        break
                annots[thisaln[0]] = PRE + SEPO + refnam
#            else:
#                annots[thisaln[0]] = UNK
            aseen.add(thisaln[0])

####################
# Add locus_tag to each gene and annotation to each transcript
scaf_order = {}
print "Ordering the predictions along the genomic scaffolds"
for entry in genes:
    this_gene = genes[entry]
    if this_gene.gene.seqid not in scaf_order:
        if this_gene.gene.strand == "-":
            scaf_order[this_gene.gene.seqid] = [[this_gene.gene.id,this_gene.gene.end]]
        else:
            scaf_order[this_gene.gene.seqid] = [[this_gene.gene.id,this_gene.gene.start]]
    else:
        if this_gene.gene.strand == "-":
            scaf_order[this_gene.gene.seqid].append([this_gene.gene.id,this_gene.gene.end])
        else:
            scaf_order[this_gene.gene.seqid].append([this_gene.gene.id,this_gene.gene.start])

locflag = 0 # flag to indicate first instance of locus ID
lokey_dict = {}
print "Adding locus tags and applying functional annotations"
for scaf in scaf_order:
#    scaf_order[scaf].sort(key = lambda x: x[1])
#    for entry in scaf_order[scaf]:
    for entry in sorted(scaf_order[scaf],key = lambda x: int(x[1])):
        gene_name = entry[0]
        if locflag == 0:
            nam = LOCUS + "_" + str(LOCS).zfill(LOCW)
            locs = LOCS + LOCJ
            locflag = 1
        else:
            nam = LOCUS + "_" + str(locs).zfill(LOCW)
            locs += LOCJ

        genes[gene_name].gene.locus_tag = nam
        lokey.write(gene_name + "\t" + nam + "\n")
        lokey_dict[gene_name] = nam
        # Add locus_tags to each gene's transcripts, with addition of isoform
        #  letter if needed
        if len(genes[gene_name].transcript) == 1:
            tr_id = genes[gene_name].transcript.keys()[0]
#            genes[gene_name].transcript[tr_id].locus_tag = nam
            genes[gene_name].transcript[tr_id].transcript.locus_tag = nam
            lokey.write(tr_id + "\t" + nam + "\n")

            if genes[gene_name].transcript[tr_id].transcript.type == "mRNA":
                if tr_id in annots:
#                    genes[gene_name].transcript[tr_id].product = annots[tr_id]
                    genes[gene_name].transcript[tr_id].transcript.product = annots[tr_id]
                else:
#                    genes[gene_name].transcript[tr_id].product = UNK
                    genes[gene_name].transcript[tr_id].transcript.product = UNK

#            for exon in range(0,len(genes[gene_name].transcript[tr_id].exons)):
#                genes[gene_name].transcript[tr_id].exons[exon].locus_tag = nam
#            if genes[gene_name].transcript[tr_id].cds:
#                for i in range(0,len(genes[gene_name].transcript[tr_id].cds)):
#                    genes[gene_name].transcript[tr_id].cds[i].locus_tag = nam
        else:
            loc_iter = 0
            for prod in sorted(genes[gene_name].transcript):
#                genes[gene_name].transcript[prod].locus_tag = nam + ISOALPHA[loc_iter]
                genes[gene_name].transcript[prod].transcript.locus_tag = nam + ISOALPHA[loc_iter]
                lokey.write(prod + "\t" + nam + ISOALPHA[loc_iter] + "\n")

                if genes[gene_name].transcript[prod].transcript.type == "mRNA":
                    if prod in annots:
#                        genes[gene_name].transcript[prod].product = annots[prod]
                        genes[gene_name].transcript[prod].transcript.product = annots[prod]
                    else:
#                        genes[gene_name].transcript[prod].product = UNK
                        genes[gene_name].transcript[prod].transcript.product = UNK

#                for i in range(0,len(genes[gene_name].transcript[prod].exons)):
#                    genes[gene_name].transcript[prod].exons[i].locus_tag = nam + ISOALPHA[loc_iter]
#                if genes[gene_name].transcript[prod].cds:
#                    for j in range(0,len(genes[gene_name].transcript[prod].cds)):
#                        genes[gene_name].transcript[prod].cds[j].locus_tag = nam + ISOALPHA[loc_iter]
                loc_iter += 1
lokey.close()

###################
#load genomic scaffolds
genome = {}
with open(fasta,"r") as file_object:
    print "Loading genomic scaffolds"
    for line in f(file_object):
        genome[line[0]] = line[1]

###################
# output tbl lines
#scafs = set()
prev_scaf = ''
# order genes by ID; hopefully this corresponds to their position on the scaffold
print "Writing final .tbl file"
for scaf in scaf_order:
#    for entry in scaf_order[scaf]:
    for entry in sorted(scaf_order[scaf],key = lambda x: int(x[1])):
        rec = entry[0]
#        for rec in sorted(genes):
#        if genes[rec].gene.seqid not in scafs:
        if genes[rec].gene.seqid != prev_scaf:
            tblout.write("".join([">Feature ",genes[rec].gene.seqid,"\n"]))
            tblout.write("\t".join(["1",str(len(genome[genes[rec].gene.seqid])),
                                            "REFERENCE\n"]))
            tblout.write("\t".join(["\t\t\tPBARC","12345\n"]))
#            scafs.add(genes[rec].gene.seqid)

        tblout.write(str(genes[rec].print_gene()))
        prev_scaf = genes[rec].gene.seqid

        # sort each gene's transcripts by name
        for prod in sorted(genes[rec].transcript.keys()):
            try:
                # if there is a note about a split feature, fix it
                if genes[rec].transcript[prod].transcript.note:
                    this_note = genes[rec].transcript[prod].transcript.note.split(" ")
                    fixed_note = []
                    for i in this_note:
                        if i in lokey_dict:
                            fixed_note.append(lokey_dict[i])
                        else:
                            fixed_note.append(re.sub("end_","end;",i))
                    genes[rec].transcript[prod].transcript.note = " ".join(fixed_note)
                tblout.write(str(genes[rec].transcript[prod].print_transcript(
                            product_type = genes[rec].transcript[prod].transcript.type,
                            outform = 'tbl',sequence=genome[genes[rec].gene.seqid])))
            except:
                print "Failed to write out "+str(genes[rec].transcript[prod].transcript.id)

tblout.close

###################
# output (repaired) gff lines
#print "Writing .gff equivalent of .tbl file"
#for rec in sorted(genes):
#    gffout.write(str(genes[rec].print_gene(outform = 'gff')))
#
#    for prod in sorted(genes[rec].transcript.keys()):
#        try:
#            gffout.write(str(genes[rec].transcript[prod].print_transcript(
#                        product_type = genes[rec].transcript[prod].transcript.type,
#                        outform = 'gff')))
#        except:
#            print "Failed to write out "+str(genes[rec].transcript[prod].transcript.id)
#
#gffout.close

### EOF ###
