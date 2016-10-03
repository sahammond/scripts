#!/usr/bin/env python

# usage: gff3-to-tbl.py genes.gff genome.fa alignment.blastp

import sys
import re
import copy

from fastatools import fasta_iter as f
import gfftools

############
# parse arguments

gff = sys.argv[1]
fasta = sys.argv[2]
#aln=sys.argv[3] # TODO this part not working yet

############
# declarations

# NCBI-assigned locus_tag
LOCUS="AB205"
# start number for locus id iteration
LOCS=50
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

gffout=open(re.sub(".gff","-prepared.gff",gff),"w")
tblout=open(re.sub(".gff","-prepared.tbl",gff),"w")

###################
# MAIN

genes = {}
# read in the entries from the gff file
with open(gff,"r") as infile:
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
# apparently some maker predictions have CDS entries but no, or too few, exon entries
# This is either my doing, or maker's. I hate maker.
# for those predictions, clone the cds to exons
for entry in genes:
    for rec in genes[entry].transcript:
        if len(genes[entry].transcript[rec].cds) > len(genes[entry].transcript[rec].exons):
            print "Replacing missing exons in "+rec
            genes[entry].transcript[rec].exons = copy.deepcopy(genes[entry].transcript[rec].cds)
            for segment in range(0,len(genes[entry].transcript[rec].exons)):
                genes[entry].transcript[rec].exons[segment].type = "exon"
                genes[entry].transcript[rec].exons[segment].offset = "."
                # fix type and offset in raw, for posterity
                genes[entry].transcript[rec].exons[segment].raw = re.sub("CDS","exon",genes[entry].transcript[rec].exons[segment].raw,flags=re.IGNORECASE)
                genes[entry].transcript[rec].exons[segment].raw = re.sub("\t[0-9]\t","\t.\t",genes[entry].transcript[rec].exons[segment].raw,flags=re.IGNORECASE)
# adjust small intron boundaries (ncbi min intron length is 10 bp)
# use steps of 3 to preserve frame
#  this will effectively delete up to 4 amino acids per adjustment (max if intron len was 1)
        for j in range(0,len(genes[entry].transcript[rec].exons)):
            try:
                dist = (int(genes[entry].transcript[rec].exons[j+1].start)
                        - int(genes[entry].transcript[rec].exons[j].stop))
                while dist < 11: # if use 'while dist < 10:' then too-small introns remain...
                    dist += 3
                    genes[entry].transcript[rec].exons[j+1].start = str(int(
                                genes[entry].transcript[rec].exons[j+1].start) + 3)
                    # adjust the corresponding CDS start too
                    try:
                        for segment in range(0,len(genes[entry].transcript[rec].cds)):
                            if (genes[entry].transcript[rec].exons[j+1].start == 
                                    genes[entry].transcript[rec].cds[segment].start):
                                genes[entry][transcript].cds[segment].start = str(int(genes[entry].transcript[rec].cds[segment].start) + 3)
                    except:
                        continue
            except:
                continue
    # make sure gene coordinates agree with new exons by setting them to the outermost
    # coordinates of all the transcripts
    estarts = set()
    eends = set()
    for rec in genes[entry].transcript:
        for segment in range(0,len(genes[entry].transcript[rec].exons)):
            estarts.add(int(re.sub("[<>]","",genes[entry].transcript[rec].exons[segment].start)))
            eends.add(int(re.sub("[<>]","",genes[entry].transcript[rec].exons[segment].end)))
    genes[entry].gene.start = str(min(estarts))
    genes[entry].gene.end = str(max(eends))

###################

# TODO make sure gene and mrna boundaries agree with each other and
#   constituent exons
###################

# THIS SECTION CURRENTLY NOT WORKING
###################
# Collect annotations for each transcript
#aseen = set()
#annots = {}
#with open(aln,"r") as blast:
#    for rec in blast:
#        thisaln=rec.split("\t")
#        if thisaln[0] not in aseen:
#            cov=thisaln[-1]
#            if thisaln[2] >= MIN_IDENT and cov >= MIN_COV:
#                thisannot=thisaln[-2].split(" ")
#                refnam=''
#                for i in range(0,len(thisannot)):
#                    if i == 0 and len(re.findall(BRK,thisannot[i])) < 1:
#                        refnam+=thisannot[i]
#                    elif i > 0 and len(re.findall(BRK,thisannot[i])) < 1:
#                        refnam+=SEPO+thisannot[i]
#                    else:
#                        break
#                annots[thisaln[0]] = PRE + SEPO + refnam
#            else:
#                annots[thisaln[0]] = UNK
#            aseen.add(thisaln[0])
#
####################
## Add locus_tag to each gene and annotation to each transcript
#scaf_order = {}
#for entry in genes:
#    this_gene = genes[entry]
#    if this_gene.gene.seqid not in scaf_order:
#        scaf_order[this_gene.gene.seqid] = [[this_gene.gene.id,this_gene.gene.start]]
#    else:
#        scaf_order[this_gene.gene.seqid].append([this_gene.gene.id,this_gene.gene.start])
#
#locflag=0 # flag to indicate first instance of locus ID
#for scaf in scaf_order:
#    scaf_order[scaf].sort(key = lambda x: x[1])
#    for entry in scaf_order[scaf]:
#        gene_name = entry[0]
#        if locflag == 0:
#            nam = LOCUS + "_" + str(LOCS).zfill(LOCW)
#            locs = LOCS + LOCJ
#            locflag = 1
#        else:
#            nam = LOCUS + "_" + str(locs).zfill(LOCW)
#            locs += LOCJ
#
#        genes[gene_name].gene.locus_tag = nam
#        # Add locus_tags to each gene's transcripts, with addition of isoform
#        #  letter if needed
#        if len(genes[gene_name].transcript) == 1:
#            tr_id = genes[gene_name].transcript.keys()[0]
#            genes[gene_name].transcript[tr_id].locus_tag = nam
#
#            if tr_id in annots:
#                genes[gene_name].transcript[tr_id].product = annots[tr_id]
#            else:
#                genes[gene_name].transcript[tr_id].product = UNK
#
#            for exon in range(0,len(genes[gene_name].transcript[tr_id].exons)):
#                genes[gene_name].transcript[tr_id].exons[exon].locus_tag = nam
#            if genes[gene_name].transcript[tr_id].cds:
#                for i in range(0,len(genes[gene_name].transcript[tr_id].cds)):
#                    genes[gene_name].transcript[tr_id].cds[i].locus_tag = nam
#        else:
#            loc_iter = 0
#            for prod in sorted(genes[gene_name].transcript):
#                genes[gene_name].transcript[prod].locus_tag = nam + ISOALPHA[loc_iter]
#
#                if prod in annots:
#                    genes[gene_name].transcript[prod].product = annots[prod]
#                else:
#                    genes[gene_name].transcript[prod].product = UNK
#
#                for i in range(0,len(genes[gene_name].transcript[prod].exons)):
#                    genes[gene_name].transcript[prod].exons[i].locus_tag = nam + ISOALPHA[loc_iter]
#                if genes[gene_name].transcript[prod].cds:
#                    for j in range(0,len(genes[gene_name].transcript[prod].cds)):
#                        genes[gene_name].transcript[prod].cds[j].locus_tag = nam + ISOALPHA[loc_iter]

###################
#load genomic scaffolds
genome = {}
with open(fasta,"r") as file_object:
    for line in f(file_object):
        genome[line[0]] = line[1]

###################
# output tbl lines
scafs = set()
# order genes by ID; hopefully this corresponds to their position on the scaffold
for rec in sorted(genes):
    if genes[rec].gene.seqid not in scafs:
        tblout.write("".join([">Feature ",genes[rec].gene.seqid,"\n"]))
        tblout.write("\t".join(["1",str(len(genome[genes[rec].gene.seqid])),
                                        "REFERENCE\n"]))
        tblout.write("\t".join(["\t\t\tPBARC","12345\n"]))
        scafs.add(genes[rec].gene.seqid)

    tblout.write(str(genes[rec].print_gene()))

    # sort each gene's transcripts by name
    for prod in sorted(genes[rec].transcript.keys()):
        try:
            tblout.write(str(genes[rec].transcript[prod].print_transcript(
                        product_type = genes[rec].transcript[prod].transcript.type,
                        outform = 'tbl',sequence=genome[genes[rec].gene.seqid])))
        except:
            print "Failed to write out "+str(genes[rec].transcript[prod].transcript.id)

tblout.close

###################
# output (repaired) gff lines
for rec in sorted(genes):
    gffout.write(str(genes[rec].print_gene(outform = 'gff')))

    for prod in sorted(genes[rec].transcript.keys()):
        try:
            gffout.write(str(genes[rec].transcript[prod].print_transcript(
                        product_type = genes[rec].transcript[prod].transcript.type,
                        outform = 'gff')))
        except:
            print "Failed to write out "+str(genes[rec].transcript[prod].transcript.id)

gffout.close

### EOF ###
