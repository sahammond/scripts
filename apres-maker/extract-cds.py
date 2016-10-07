#!/usr/bin/env python

import sys

import gfftools
import fastatools

def main():
    """Extract CDS coordinates for each mRNA in a gff file,
    then print the corresponding sequence.
    The subfeatures MUST be sorted by their start coordinate

    Usage: extract-cds.py genes.gff genome.fa > genes-cds.fa
    """

    if len(sys.argv) is not 3:
        print main.__doc__
        sys.exit(1)

    gff=sys.argv[1]
    fasta=sys.argv[2]

    cds_coords = {}

    with open(gff,"r") as file_object:
        for line in file_object:
            if line[0] == "#":
                sys.stdout.write(line)
                continue

            feature = gfftools.GFF(line)

            if feature.type == "CDS":
                for parent in feature.parent:
                    if parent in cds_coords:
                        if feature.strand == "+":
                            cds_coords[parent]['coords'].append([feature.start,feature.end])
    
                        else:
                            cds_coords[parent]['coords'].insert(0,[feature.start,feature.end])
    
                    else:
                        cds_coords[parent] = {'coords':[[feature.start,feature.end]],'scaf':feature.seqid,'strand':feature.strand}

    sequences={}

    with open(fasta,"r") as file_object:
        for rec in fastatools.fasta_iter(file_object):
            sequences[rec[0]] = rec[1]

    for record in cds_coords:
        whole = ""
        for coords in cds_coords[record]['coords']:
            seqid = cds_coords[record]['scaf']

            if cds_coords[record]['strand'] == "-":
                piece = fastatools.revcomp(sequences[seqid][int(coords[0])-1:int(coords[1])]) ### watch for OBO

            else:
                piece = sequences[seqid][int(coords[0])-1:int(coords[1])] ### watch for OBO

            whole += piece
        else:
            sys.stdout.write(">"+record+"_CDS\n"+whole+"\n")

if __name__ == "__main__":
    main()

### EOF ###
