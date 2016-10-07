#!/usr/bin/env python

import sys
import re

from gfftools import GFF as GFF

gff = sys.argv[1]
multi = {} # ID: count

with open(gff,"r") as file_object:
    for line in file_object:

        if line[0] == "#":
            continue

        feat = GFF(line)

        if feat.source == "gmap" and feat.type == "mRNA":

            cch = feat.id.split("-")
            deriv = ''
            x = 0

            for i in cch:
                if i == "CCH":
                    deriv = i
                    x = 1
                elif i == "mRNA" or i == "lncRNA":
                    x = 0
                elif x == 1:
                    deriv += "-" + i

            if deriv not in multi:
                multi[deriv] = 1
            else:
                multi[deriv] += 1

mates={}
with open(gff,"r") as file_object:
    for line in file_object:
        if line[0] == "#":
            continue
        feat = GFF(line)
        if feat.source == "maker":
            continue
        elif feat.source == "gmap":
            cch = feat.id.split("-")
            deriv = ''
            x = 0
            for i in cch:
                if i == "CCH":
                    deriv = i
                    x = 1
                elif i == "mRNA" or i == "lncRNA" or i == "gene":
                    x = 0
                elif x == 1:
                    deriv += "-" + i
            if len(re.findall("mRNA-1",feat.id)) > 0:
                mates[deriv]={'1':[re.sub('-mRNA-1','',feat.id.split(":")[0]),feat.seqid]}
            elif len(re.findall("mRNA-2",feat.id)) > 0:
                mates[deriv]['2']=[re.sub('-mRNA-2','',feat.id.split(":")[0]),feat.seqid]

with open(gff,"r") as file_object:
    for line in file_object:
        if line[0] == "#":
            sys.stdout.write(line)
            continue
        feat = GFF(line)
        if feat.source == "maker":
            sys.stdout.write(line)
            continue
        elif feat.source == "gmap":
            cch = feat.id.split("-")
            deriv = ''
            x = 0
            for i in cch:
                if i == "CCH":
                    deriv = i
                    x = 1
                elif i == "mRNA" or i == "lncRNA" or i == "gene":
                    x = 0
                elif x == 1:
                    deriv += "-" + i
            if feat.type != "mRNA" and feat.type != "ncRNA": ##
                sys.stdout.write(line)
            else:
                if '2' not in mates[deriv]:
                    sys.stdout.write(line)
                else:
                    if len(re.findall("mRNA-1",feat.id)) > 0:
                        sys.stdout.write(line.strip("\n")+";note=5' end_ 3' end is gene "+mates[deriv]['2'][0]+" on scaffold "+mates[deriv]['2'][1]+"\n")
                    elif len(re.findall("mRNA-2",feat.id)) > 0:
                        sys.stdout.write(line.strip("\n")+";note=3' end_ 5' end is gene "+mates[deriv]['1'][0]+" on scaffold "+mates[deriv]['1'][1]+"\n")

### EOF ###
