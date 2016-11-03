#!/usr/bin/env python

import sys
import re

def main():

    """Fix gene IDs in .tbl file notes
    Running gff3-to-tbl.py piecewise on a gff file that includes
    annotations split onto two scaffolds can result in notes
    that contain the original gene IDs instead of the locus_tag.
    This script will replace the gene IDs in these instances
    with their locus_tags.

    Usage: repair-tbl-notes.py annotations.tbl ID-conversion-table.tsv
    Output: annotations-fixedNote.tbl

    Note: ID-conversion-table.tsv must include all the locus tags
        from all pieces of the original gff file.
    """

    tbl = sys.argv[1]
    ids = sys.argv[2]

    convo = {}
    with open(ids,"r") as file_object:
        print "Loading gene ID - locus_tag pairs"
        for line in file_object:
            if re.findall("mRNA",line):
                continue
            else:
                rec = line.strip("\n").split("\t")
                convo[rec[0]] = rec[1]

    outnam = re.sub(".tbl","-fixedNote.tbl",tbl)
    outfile = open(outnam,"w")

    with open(tbl,"r") as file_object:
        print "Reading .tbl file and converting IDs as needed"
        print "Output file will be called " + outnam
        for line in file_object:
            if re.findall("note",line):
                rec = line.strip("\n").split(" ")
                outbuff = []
                for word in rec:
                    if word not in convo:
                        outbuff.append(word)
                    else:
                        outbuff.append(convo[word])
                else:
                    outbuff[-1] += "\n"
                outfile.write(" ".join(outbuff))
            else:
                outfile.write(line)
    outfile.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print main.__doc__
        sys.exit(1)
    main()

### EOF ###
