#!/usr/bin/env python

# usage: add-alt-splice-note.py genome.tbl > genome-altSpliceNote.tbl
import sys
import re

def main():

    tbl = sys.argv[1]

    # read through the table and make the necessary changes
    with open(tbl,"r") as file_object:
        for line in file_object:
            if len(line) < 2:
                continue
            line = re.sub("end_","end;",line)
            # isoforms (via prepare_gff.py) are identified by trailing A-Z
            if re.findall("protein_id",line):
                if re.findall('[A-Z]',line.strip("\n")[-1]):
                    sys.stdout.write("\t\t\tnote\talternatively spliced\n")
                    sys.stdout.write(line)

                else:
                    sys.stdout.write(line)
            else:
                sys.stdout.write(line)
    
if __name__ == "__main__":
    main()

### EOF ###
