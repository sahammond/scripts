#!/usr/bin/env python

import sys
import re

gff=sys.argv[1]
trans=sys.argv[2]

convo={}
revconvo={}

with open(trans,"r") as file_object:
	for line in file_object:
		rec=line.strip("\n").split("\t")
		if rec[0] == "gene":
			convo[rec[2]]=rec[1]

for key in convo:
#	revconvo[convo[key]]=re.sub('mRNA-[0-9]','',key) # basically gene id from mRNA
	revconvo[convo[key]]=key

with open(gff,"r") as file_object:
	for line in file_object:
		rec=line.strip("\n").split("\t")
		if rec[2] == "mRNA":
			comm=re.sub("ID=",'',rec[8].split(";")[0])
			if len(comm) > 13:
				comm=comm[:13]
			if len(re.findall("gmap",convo[comm])) > 0:
				toout=re.sub(";product=hypothetical protein","",line)
				toout=re.sub("-mRNA-","-lncRNA-",toout)
				toout=re.sub("mRNA","ncRNA",toout)

				if len(re.findall('end is gene',line)) > 0: #catch split records and convert ID
					spacecomm=rec[8].split(" ") ##
					for part in spacecomm: ##
						if len(re.findall('gmap',part)) > 0: ##
							mateid=re.sub('-mRNA-[0-9]','',part.split(":")[0]) ##
							break ##
					toout=re.sub(mateid,revconvo[mateid],toout) ##

				sys.stdout.write(toout)
			else:
				sys.stdout.write(line)
		else:
			sys.stdout.write(line)
quit()
