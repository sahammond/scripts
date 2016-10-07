#!/usr/bin/env python

import sys
import re

gff=sys.argv[1]

with open(gff,"r") as file_object:
    for line in file_object:

        if line[0] == "#":
            sys.stdout.write(line)
            continue


        if re.findall('ID=CCH',line):
            rec=line.strip("\n").split("\t")
    
            baseID=rec[8].split(";")[0].split("-")[:3]
            baseID[0]=re.sub("ID=","",baseID[0])
            baseID[2]=re.sub('\.mrna[0-9]','',baseID[2])
            baseID[2]=re.sub('\.exon[0-9]','',baseID[2])
            baseID[2]=re.sub('.cds','',baseID[2])
    
    
            newID=''    
            for i in baseID:
                newID+=i+"-"
    
            newID=newID[:-1]
            newID="gmap-"+rec[0]+"-"+re.sub("ID=","",newID)
        
            oldID=''
            for i in baseID:
                oldID+=i+"-"
            oldID=oldID[:-1]
        
            for field in range(0,len(rec)):
                rec[field]=re.sub(oldID,newID,rec[field])
                rec[field]=re.sub('\.mrna','-mRNA-',rec[field])
                rec[field]=re.sub('\.exon',':exon:',rec[field])
                rec[field]=re.sub('.cds',':cds',rec[field])
        
            for col in range(0,(len(rec)-1)):
                sys.stdout.write(rec[col]+"\t")
            else:
                sys.stdout.write(rec[-1]+"\n")

        else:
            sys.stdout.write(line)
            continue

quit()
