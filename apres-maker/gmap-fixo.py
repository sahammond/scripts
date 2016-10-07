#!/usr/bin/env python

import sys
import re

for line in sys.stdin:
	rec=line.strip("\n").split("\t")
	if rec[2] == "gene":
		sys.stdout.write(line)
	else:
		for i in rec[0:8]:
			sys.stdout.write(i+"\t")
		else:
			comm=rec[8].split(";")
			x=0
			for j in comm:
				if x == 1:
					x=0
					break
				elif len(re.findall("ID",j)) > 0:
					sys.stdout.write(j+";")
				elif len(re.findall("Parent",j)) > 0:
					sys.stdout.write(j+"\n")
					x=1
quit()
