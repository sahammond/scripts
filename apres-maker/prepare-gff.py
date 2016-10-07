#!/home/shammond/arch/CentOS_5/bin/python

# usage: sort-gff.py genemodels.gff alignment.blastp swissprot.fa

import sys as s
import re
from fastatools import fasta_iter as f

############
# parse arguments

gff=s.argv[1]
aln=s.argv[2]
ref=s.argv[3]

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
PRE='product=similar to' # string to prefix description with ### modified for definitive naming of genevalidator 90+ HC proteins
UNK='product=hypothetical protein' # label to give queries that have no good annotation
ISOALPHA=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'] # alphabet to use for isoform labeling
CTAB=open('maker-to-ncbi-prep-ID-table.txt',"w") # file to write gene and mRNA conversion table to
SCAFL=len('Rc-01r160223s0000001') # length of your scaffold IDs. must be all the same. must be integer.

############
# declare empty objects

mrna={} # scaffold -> mRNA -> [exons],[cds]
seen=set() # set to track observed mRNA IDs
gene=[] # list of final modified gff lines
gseen=set() # set to track observed gene lines
gbound={} # keep the min and max positions across isoforms per gene
gscaff=set() # set to track observed genomic scaffold IDs
gseenagain=set() # another set to track observed gene lines
present=set() # set to track mRNA IDs ### how is it diff from 'seen'?
refnam={} # dict of reference descriptions
locflag=0 # flag to indicate first instance of locus ID
scarfseen=set() # set to track scaffolds as their gene boundaries are corrected
scarfs={} # dict to hold gene boundaries are they're ordered along the scaffolds
glocus={} # dict to hold gene info (scaffold,boundaries,locus_tag)
aseen=set() # set to track observed alignments
isog={} # dict to track total and 'observed' mRNA isoforms key = gene, value = number of assoc mrnas

############
# define functions

# function to derive scaffold from maker-type mrna or gene IDs
#  will exit if it can't trim off any of the strings in 'prefixes'
prefixes=['maker-','augustus_masked-','genemark-','snap_masked-','gmap-']
def getScaf(seqid,idlen):
    for pref in prefixes:
        preflen=len(pref)
        if seqid[0:preflen] == pref:
            return seqid[preflen:idlen+preflen]
    else:
        print >> s.stderr, 'getScaf failed to derive a scaffold ID from a gene or mRNA ID.'
        print >> s.stderr, 'Please check "prefixes" in '+s.argv[0]
        s.exit(1)

############
# MAIN

# read in the mrna entries from the gff file
with open(gff,"r") as infile:
    for line in infile:
        if line[0] == "#" or line[0] == "-":
            continue
        rec=line.strip("\n").split("\t")
        if rec[2] == "mRNA":
            # get the mRNA's name
            mid=re.sub("ID=",'',rec[-1].split(';')[0])
            if mid not in seen:
                if len(re.findall('note=',line)) > 0:
                    for i in rec[-1].split(";"):
                        if len(re.findall('note=',i)) > 0:
                            mrna[mid]={'note':i.split(";")[0]}
                else:
                    mrna[mid]={}
                seen.add(mid)

# read in the other entries from the gff
with open(gff,"r") as infile:
    for line in infile:
        if line[0] == "#" or line[0] == "-":
            continue
        rec=line.strip("\n").split("\t")
        if rec[2] == "CDS":
            # get the parent mRNA's name
            # 'Parent=' field is second in the comment column.
            comm=rec[-1].split(";")
##          comm=rec[8].split(";")

            # exons and CDS can have more than 1 mRNA as a parent
            mult=comm[1][7:].split(",")

            for papa in range(0,len(mult)):
                mid=mult[papa]

                pos=[int(rec[3]),int(rec[4]),rec[6],rec[7]]

                if mid in seen:
                    if 'cds' not in mrna[mid].keys():
                        mrna[mid]['cds']=[pos]
                    else:
                        t=0
                        for i in mrna[mid]['cds']:
                            if pos[0] < i[0]:
                                mrna[mid]['cds'].insert(t,pos)
                                break
                        else:
                            mrna[mid]['cds'].append(pos)
                else:
                    mrna[mid]={'cds':[pos]}
                    seen.add(mid)

        elif rec[2] == "exon":
            # 'Parent=' field is second in the comment column.
            comm=rec[-1].split(";")
#           comm=rec[8].split(";")

            # exons and CDS can have more than 1 mRNA as a parent
            mult=comm[1][7:].split(",")
#           for field in comm:
#               if len(re.findall('Parent',field)) > 0:
#                   mult=re.sub('Parent=','',field).split(",")
#           else:
#               break
            for papa in range(0,len(mult)):
                mid=mult[papa]
            
                pos=[int(rec[3]),int(rec[4]),rec[6],rec[7]]

                if mid in seen:
                    if 'exon' not in mrna[mid].keys():
                        mrna[mid]['exon']=[pos]
                    else:
                        t=0
                        for i in mrna[mid]['exon']:
                            if pos[0] < i[0]:
                                mrna[mid]['exon'].insert(t,pos)
                                break
                        else:
                            mrna[mid]['exon'].append(pos)
                else:
                    mrna[mid]={'exon':[pos]}
                    seen.add(mid)

### OLD. Now dealt with later on.
# apparently some maker predictions have CDS entries but no exon entries. This is either my doing, or maker's. I hate maker.
# for those predictions, clone the cds to exons
#for i in mrna:
#   if 'exon' not in mrna[i]:
#       mrna[i]['exon']=mrna[i]['cds']

### OLD. Now dealt with later on.
# adjust small intron boundaries
#  ncbi min intron length is 10 bp
# use steps of 3 to preserve frame
#  this will effectively delete 1-4 amino acids per adjustment
#for i in mrna:
#   mrna[i]['exon'].sort(key=lambda q: q[0]) # sort by start position
#
#   if 'cds' in mrna[i].keys():
#       mrna[i]['cds'].sort(key=lambda q: q[0])
#
#   for j in range(0,len(mrna[i]['exon'])):
#       try:
#           dist=int(mrna[i]['exon'][j+1][0] - mrna[i]['exon'][j][1])
#           while dist < 11: # if use 'while dist < 10:' then too-small introns remain...
#               dist+=3
#               mrna[i]['exon'][j+1][0]+=3
#               try:
#                   mrna[i]['cds'][j+1][0]+=3 # adjust cds too!
#               except:
#                   continue
#       except:
#            continue

# get set of scaffold IDs
for i in seen:
    if len(i) < 2:
        continue
    # extract scaffold ID from mRNA ID
    t=''
    for j in i.split("-"):
        if j != "mRNA":
            t+=j+"-"
        else:
            gid=t[:-1]

    gscaff.add(getScaf(gid,SCAFL))

    # get positions of predictions along each scaffold
    if gid not in gseen:
        # add record to gbound (L,R bound)
        gbound[gid]=[mrna[i]['exon'][0][0],mrna[i]['exon'][-1][1]]

        gseen.add(gid)

    else:
        # check gbound record
        if gbound[gid][0] > mrna[i]['exon'][0][0]:
            gbound[gid][0] = mrna[i]['exon'][0][0]
        if gbound[gid][1] < mrna[i]['exon'][-1][1]:
            gbound[gid][1] = mrna[i]['exon'][-1][1]

# order the genes along the scaffolds
for i in gbound:
    # think of a more automatic way to do this
    scafid=getScaf(i,SCAFL) 
    if scafid not in scarfseen:
        scarfs[scafid]=[[gbound[i],i]]
        scarfseen.add(scafid)

    else:
        scarfs[scafid].append([gbound[i],i])

# append locus id
for i in scarfs:

    scarfs[i].sort(key=lambda q: q[0])

    for j in range(len(scarfs[i])):
        if locflag == 0:
            nam=LOCUS+"_"+str(LOCS).zfill(LOCW)
            locs=LOCS+LOCJ
            locflag=1

        else:
            nam=LOCUS+"_"+str(locs).zfill(LOCW)
            locs+=LOCJ

        glocus[scarfs[i][j][1]]=[i,scarfs[i][j][0],nam]

# evaluate alignments, use those that pass defined thresholds to annotate contigs
#  otherwise give them the label in UNK
# get set of mrnas under consideration
for i in mrna:
    present.add(i)

# load reference descriptions
#with open(ref,"r") as fa:
#   for rec in f(fa):
#       title=rec[0].split(SEPR)
#       toap=''
#       i=0
#       for word in title:
#           if word[:len(BRK)] == BRK:
#               break
#           elif word[0] != ">" and i == 1:
#               toap=str(word)
#           elif word[0] != ">" and i > 1:
#               toap=toap+SEPO+str(word)
#           i+=1
#       # store seqnam as key, description as value
#       refnam[title[0]]=toap 

# add annotations to mrna lines
# GAG will add them to the cds records in the tbl file
with open(aln,"r") as blast:
    for rec in blast:
        # qid = 0, sid = 1, pid = 2, aln length = 3 (includes gaps)
        thisaln=rec.split("\t")
        if thisaln[0] not in aseen:
            if thisaln[0] in present:
#               cov=100*(float(thisaln[3])/float(len(mrna[thisaln[0]])))
                cov=thisaln[-1]
                if thisaln[2] >= MIN_IDENT and cov >= MIN_COV:
                    thisannot=thisaln[-2].split(" ") ##
                    refnam='' ##
                    for i in range(0,len(thisannot)): ##
                        if i == 0 and len(re.findall(BRK,thisannot[i])) < 1: ##
                            refnam+=thisannot[i] ##
                        elif i > 0 and len(re.findall(BRK,thisannot[i])) < 1: ##
                            refnam+=SEPO+thisannot[i] ##
                        else: ##
                            break ##
                    #mrna[thisaln[0]]['annot']=PRE+str(refnam[thisaln[1]])
                    mrna[thisaln[0]]['annot']=PRE+SEPO+refnam
                else:
                    mrna[thisaln[0]]['annot']=UNK
                aseen.add(thisaln[0])

# add UNK to the contigs with no hits at all
for i in mrna:
    if i not in aseen:
        mrna[i]['annot']=UNK

# figure out how many isoforms are associated with each gene
#  need this info because mRNAs have to have form of 'isoform A', not 'isoform 1'
for i in mrna:
    # reproduce gene line
    t=''
    for j in i.split("-"):
        if j != "mRNA":
            t+=j+"-"
        else:
            gid=t[:-1]
    
    if gid not in isog:
        isog[gid]=[0,1]
    else:
        isog[gid][1]+=1

# prepare edited gff lines for output
# GAG take cares of the transcript and protein ID
for i in mrna:
    # reproduce gene line
    t=''
    for j in i.split("-"):
        if j != "mRNA":
            t+=j+"-"
        else:
            gid=t[:-1]

    scaf=getScaf(gid,SCAFL)

    gline=[scaf,"maker","gene",glocus[gid][1][0],glocus[gid][1][1],".",mrna[i]['exon'][0][2],".","ID="+glocus[gid][2]]

    # handle extra mrna isoforms
    if isog[gid][1] > 1:
        if 'note' in mrna[i]:
#           newnote=re.sub(mrna[i]['note'],
            mline=[scaf,"maker","mRNA",mrna[i]['exon'][0][0],mrna[i]['exon'][-1][1],".",mrna[i]['exon'][0][2],".","ID="+glocus[gid][2]+ISOALPHA[isog[gid][0]]+";Parent="+glocus[gid][2]+";"+mrna[i]['annot']+";"+mrna[i]['note']]
        else:
            mline=[scaf,"maker","mRNA",mrna[i]['exon'][0][0],mrna[i]['exon'][-1][1],".",mrna[i]['exon'][0][2],".","ID="+glocus[gid][2]+ISOALPHA[isog[gid][0]]+";Parent="+glocus[gid][2]+";"+mrna[i]['annot']]

        # write out to ID conversion table
        CTAB.write("mRNA"+"\t"+i+"\t"+glocus[gid][2]+ISOALPHA[isog[gid][0]]+"\n")

    else:
        if 'note' in mrna[i]:
            mline=[scaf,"maker","mRNA",mrna[i]['exon'][0][0],mrna[i]['exon'][-1][1],".",mrna[i]['exon'][0][2],".","ID="+glocus[gid][2]+";Parent="+glocus[gid][2]+";"+mrna[i]['annot']+";"+mrna[i]['note']]
        else:
            mline=[scaf,"maker","mRNA",mrna[i]['exon'][0][0],mrna[i]['exon'][-1][1],".",mrna[i]['exon'][0][2],".","ID="+glocus[gid][2]+";Parent="+glocus[gid][2]+";"+mrna[i]['annot']]

        # write out to ID conversion table
        CTAB.write("mRNA"+"\t"+i+"\t"+glocus[gid][2]+"\n")

    if gid not in gseenagain:
        # add recapitulated gene line and mrna line to 'gene'
        gene.append(gline)
        gene.append(mline)
        gseenagain.add(gid)

        # write out to ID conversion table
        CTAB.write("gene"+"\t"+gid+"\t"+glocus[gid][2]+"\n")

    else:
        gene.append(mline)

    # add exons and cds
    if isog[gid][1] > 1:
        for e in range(0,len(mrna[i]['exon'])):
            gene.append([scaf,"maker","exon",mrna[i]['exon'][e][0],mrna[i]['exon'][e][1],".",mrna[i]['exon'][e][2],mrna[i]['exon'][e][3],"ID="+glocus[gid][2]+ISOALPHA[isog[gid][0]]+":exon;Parent="+glocus[gid][2]+ISOALPHA[isog[gid][0]]])
        if 'cds' in mrna[i].keys():
            for c in range(0,len(mrna[i]['cds'])):
                gene.append([scaf,"maker","CDS",mrna[i]['cds'][c][0],mrna[i]['cds'][c][1],".",mrna[i]['cds'][c][2],mrna[i]['cds'][c][3],"ID="+glocus[gid][2]+ISOALPHA[isog[gid][0]]+":cds;Parent="+glocus[gid][2]+ISOALPHA[isog[gid][0]]])

        # iterate the isoform label
        #  do it here, and not while reconstructing the mrna line,
        #    because every mrna will have associated exon(s) and cds(s)
        isog[gid][0]+=1
    
    else:
        for e in range(0,len(mrna[i]['exon'])):
            gene.append([scaf,"maker","exon",mrna[i]['exon'][e][0],mrna[i]['exon'][e][1],".",mrna[i]['exon'][e][2],mrna[i]['exon'][e][3],"ID="+glocus[gid][2]+":exon;Parent="+glocus[gid][2]])
        if 'cds' in mrna[i].keys():
            for c in range(0,len(mrna[i]['cds'])):
                gene.append([scaf,"maker","CDS",mrna[i]['cds'][c][0],mrna[i]['cds'][c][1],".",mrna[i]['cds'][c][2],mrna[i]['cds'][c][3],"ID="+glocus[gid][2]+":cds;Parent="+glocus[gid][2]])

CTAB.close()
    
# write out edited lines
for i in gene:
    for j in range(0,len(i)-1):
        s.stdout.write(str(i[j])+"\t")
    else:
        s.stdout.write(str(i[-1])+"\n")
### EOF ###
