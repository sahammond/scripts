#!/usr/bin/env python

# usage: gff3-to-tbl.py genes.gff genome.fa

import sys
import re
import copy

from fastatools import fasta_iter as f
import gfftools

############
# parse arguments

gff = sys.argv[1]
#aln=s.argv[2]
#ref=s.argv[3]
ref = sys.argv[2]

############
# declarations

LOCUS="AB205" # NCBI-assigned locus_id
LOCS=50 # start number for locus id iteration
LOCW=7 # width for locus ids (i.e. total with padding)
LOCJ=10 # jump size between loci (good idea to allow for genes identified between current ones)
MIN_IDENT=25 # miniumum percent identity to accept alignment for annotation
MIN_COV=50 # minimum percent coverage to accept alignment for annotation
BRK='OS=' # string to split reference protein names with. Script expects swissprot-style protein names in reference.fa, use 'OS='
SEPQ=' ' # query fasta title line field sep
SEPR=' ' # reference fasta title line field sep
SEPO=' ' # field sep for output fasta
PRE='product=similar to' # string to prefix description with
UNK='product=hypothetical protein' # label to give queries that have no good annotation
ISOALPHA=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'] # alphabet to use for isoform labeling
#CTAB=open('maker-to-ncbi-prep-ID-table.txt',"w") # file to write gene and mRNA conversion table to
#SCAFL=len('Rc-01r160223s0000001') # length of your scaffold IDs. must be all the same. must be integer.

############
# declare empty objects

#mrna={} # scaffold -> mRNA -> [exons],[cds]
#seen=set() # set to track observed mRNA IDs
#gene=[] # list of final modified gff lines
#gseen=set() # set to track observed gene lines
#gbound={} # keep the min and max positions across isoforms per gene
#gscaff=set() # set to track observed genomic scaffold IDs
#gseenagain=set() # another set to track observed gene lines
#present=set() # set to track mRNA IDs ### how is it diff from 'seen'?
#refnam={} # dict of reference descriptions
#locflag=0 # flag to indicate first instance of locus ID
#scarfseen=set() # set to track scaffolds as their gene boundaries are corrected
#scarfs={} # dict to hold gene boundaries are they're ordered along the scaffolds
#glocus={} # dict to hold gene info (scaffold,boundaries,locus_tag)
#aseen=set() # set to track observed alignments
#isog={} # dict to track total and 'observed' mRNA isoforms key = gene, value = number of assoc mrnas

gffout=open(re.sub(".gff","-prepared.gff",gff),"w")
tblout=open(re.sub(".gff","-prepared.tbl",gff),"w")

############
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
            for parent in feature.parent:
                for gene in genes:
                    if re.findall(genes[gene].id,parent):
                        genes[gene].transcript[parent].add_exon(feature)
#                        break
#                else:
#                    print "Failed to find gene for "+feature.id

        elif feature.type == "CDS":
            for parent in feature.parent:
                for gene in genes:
                    if re.findall(genes[gene].id,parent):
                        genes[gene].transcript[parent].add_cds(feature)
#                        break
#                else:
#                    print "Failed to find gene for "+feature.id

# TODO make sure gene coordinates agree with new exons 
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
# adjust small intron boundaries
#  ncbi min intron length is 10 bp
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
    for rec in genes[entry].transcript:
        estarts = set()
        eends = set()
        for segment in range(0,len(genes[entry].transcript[rec].exons)):
            estarts.add(int(genes[entry].transcript[rec].exons[segment].start))
            eends.add(int(genes[entry].transcript[rec].exons[segment].end))
        genes[entry].gene.start = str(min(estarts))
        genes[entry].gene.end = str(max(eends))

###################

# TODO append locus id
# TODO add annotations
# TODO make sure gene and mrna boundaries agree with each other and
#   constituent exons
#for i in scarfs:
#
#    scarfs[i].sort(key=lambda q: q[0])
#
#    for j in range(len(scarfs[i])):
#        if locflag == 0:
#            nam=LOCUS+"_"+str(LOCS).zfill(LOCW)
#            locs=LOCS+LOCJ
#            locflag=1
#
#        else:
#            nam=LOCUS+"_"+str(locs).zfill(LOCW)
#            locs+=LOCJ
#
#        glocus[scarfs[i][j][1]]=[i,scarfs[i][j][0],nam]
#
## evaluate alignments, use those that pass defined thresholds to annotate contigs
##  otherwise give them the label in UNK
## get set of mrnas under consideration
#for i in mrna:
#    present.add(i)
#
## add annotations to mrna lines
## GAG will add them to the cds records in the tbl file
#with open(aln,"r") as blast:
#    for rec in blast:
#        # qid = 0, sid = 1, pid = 2, aln length = 3 (includes gaps)
#        thisaln=rec.split("\t")
#        if thisaln[0] not in aseen:
#            if thisaln[0] in present:
##               cov=100*(float(thisaln[3])/float(len(mrna[thisaln[0]])))
#                cov=thisaln[-1]
#                if thisaln[2] >= MIN_IDENT and cov >= MIN_COV:
#                    thisannot=thisaln[-2].split(" ") ##
#                    refnam='' ##
#                    for i in range(0,len(thisannot)): ##
#                        if i == 0 and len(re.findall(BRK,thisannot[i])) < 1: ##
#                            refnam+=thisannot[i] ##
#                        elif i > 0 and len(re.findall(BRK,thisannot[i])) < 1: ##
#                            refnam+=SEPO+thisannot[i] ##
#                        else: ##
#                            break ##
#                    #mrna[thisaln[0]]['annot']=PRE+str(refnam[thisaln[1]])
#                    mrna[thisaln[0]]['annot']=PRE+SEPO+refnam
#                else:
#                    mrna[thisaln[0]]['annot']=UNK
#                aseen.add(thisaln[0])
#
## add UNK to the contigs with no hits at all
#for i in mrna:
#    if i not in aseen:
#        mrna[i]['annot']=UNK
#
## figure out how many isoforms are associated with each gene
##  need this info because mRNAs have to have form of 'isoform A', not 'isoform 1'
#for i in mrna:
#    # reproduce gene line
#    t=''
#    for j in i.split("-"):
#        if j != "mRNA":
#            t+=j+"-"
#        else:
#            gid=t[:-1]
#    
#    if gid not in isog:
#        isog[gid]=[0,1]
#    else:
#        isog[gid][1]+=1

#load genomic scaffolds
genome = {}
with open(ref,"r") as file_object:
    for line in f(file_object):
        genome[line[0]] = line[1]

# output tbl lines
scafs = set()
for rec in sorted(genes):
    if genes[rec].gene.seqid not in scafs:
        tblout.write("".join([">FEATURE ",genes[rec].gene.seqid,"\n"]))
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
                        outform = 'tbl')))
        except:
            print "Failed to write out "+str(genes[rec].transcript[prod].transcript.id)

tblout.close

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

## prepare edited gff lines for output
## GAG take cares of the transcript and protein ID
#for i in mrna:
#    # reproduce gene line
#    t=''
#    for j in i.split("-"):
#        if j != "mRNA":
#            t+=j+"-"
#        else:
#            gid=t[:-1]
#
#    scaf=getScaf(gid,SCAFL)
#
#    gline=[scaf,"maker","gene",glocus[gid][1][0],glocus[gid][1][1],".",mrna[i]['exon'][0][2],".","ID="+glocus[gid][2]]
#
#    # handle extra mrna isoforms
#    if isog[gid][1] > 1:
#        if 'note' in mrna[i]:
##           newnote=re.sub(mrna[i]['note'],
#            mline=[scaf,"maker","mRNA",mrna[i]['exon'][0][0],mrna[i]['exon'][-1][1],".",mrna[i]['exon'][0][2],".","ID="+glocus[gid][2]+ISOALPHA[isog[gid][0]]+";Parent="+glocus[gid][2]+";"+mrna[i]['annot']+";"+mrna[i]['note']]
#        else:
#            mline=[scaf,"maker","mRNA",mrna[i]['exon'][0][0],mrna[i]['exon'][-1][1],".",mrna[i]['exon'][0][2],".","ID="+glocus[gid][2]+ISOALPHA[isog[gid][0]]+";Parent="+glocus[gid][2]+";"+mrna[i]['annot']]
#
#        # write out to ID conversion table
#        CTAB.write("mRNA"+"\t"+i+"\t"+glocus[gid][2]+ISOALPHA[isog[gid][0]]+"\n")
#
#    else:
#        if 'note' in mrna[i]:
#            mline=[scaf,"maker","mRNA",mrna[i]['exon'][0][0],mrna[i]['exon'][-1][1],".",mrna[i]['exon'][0][2],".","ID="+glocus[gid][2]+";Parent="+glocus[gid][2]+";"+mrna[i]['annot']+";"+mrna[i]['note']]
#        else:
#            mline=[scaf,"maker","mRNA",mrna[i]['exon'][0][0],mrna[i]['exon'][-1][1],".",mrna[i]['exon'][0][2],".","ID="+glocus[gid][2]+";Parent="+glocus[gid][2]+";"+mrna[i]['annot']]
#
#        # write out to ID conversion table
#        CTAB.write("mRNA"+"\t"+i+"\t"+glocus[gid][2]+"\n")
#
#    if gid not in gseenagain:
#        # add recapitulated gene line and mrna line to 'gene'
#        gene.append(gline)
#        gene.append(mline)
#        gseenagain.add(gid)
#
#        # write out to ID conversion table
#        CTAB.write("gene"+"\t"+gid+"\t"+glocus[gid][2]+"\n")
#
#    else:
#        gene.append(mline)
#
#    # add exons and cds
#    if isog[gid][1] > 1:
#        for e in range(0,len(mrna[i]['exon'])):
#            # TODO: stop printing "maker", allow for gmap 
#            gene.append([scaf,"maker","exon",mrna[i]['exon'][e][0],mrna[i]['exon'][e][1],".",mrna[i]['exon'][e][2],".","ID="+glocus[gid][2]+ISOALPHA[isog[gid][0]]+":exon;Parent="+glocus[gid][2]+ISOALPHA[isog[gid][0]]])
#        if 'cds' in mrna[i].keys():
#            for c in range(0,len(mrna[i]['cds'])):
#                gene.append([scaf,"maker","CDS",mrna[i]['cds'][c][0],mrna[i]['cds'][c][1],".",mrna[i]['cds'][c][2],mrna[i]['cds'][c][3],"ID="+glocus[gid][2]+ISOALPHA[isog[gid][0]]+":cds;Parent="+glocus[gid][2]+ISOALPHA[isog[gid][0]]])
#
#        # iterate the isoform label
#        #  do it here, and not while reconstructing the mrna line,
#        #    because every mrna will have associated exon(s) and cds(s)
#        isog[gid][0]+=1
#    
#    else:
#        for e in range(0,len(mrna[i]['exon'])):
#            # TODO: stop printing "maker", allow for gmap 
#            gene.append([scaf,"maker","exon",mrna[i]['exon'][e][0],mrna[i]['exon'][e][1],".",mrna[i]['exon'][e][2],".","ID="+glocus[gid][2]+":exon;Parent="+glocus[gid][2]])
#        if 'cds' in mrna[i].keys():
#            for c in range(0,len(mrna[i]['cds'])):
#                gene.append([scaf,"maker","CDS",mrna[i]['cds'][c][0],mrna[i]['cds'][c][1],".",mrna[i]['cds'][c][2],mrna[i]['cds'][c][3],"ID="+glocus[gid][2]+":cds;Parent="+glocus[gid][2]])
#
#CTAB.close()
#    
## write out edited lines as gff
#for i in gene:
#    for j in range(0,len(i)-1):
#        gffout.write(str(i[j])+"\t")
#    else:
#        gffout.write(str(i[-1])+"\n")
#
#gffout.close()

### EOF ###
