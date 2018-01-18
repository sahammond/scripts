#!/usr/bin/env python

import argparse
import shutil
import os
import glob


def cleano(script):
    # delete jupyter notebook marks
    sin = '# In['
    for line in script:
        if line[:5] == sin:
            continue
        else:
            yield line


def nodev(script):
    # delete dev cells
    sin1 = '### FOR DEV ###\n'
    sin2 = '### END DEV ###\n'
    sinner = False
    for line in script:
        if sinner:
            if sin2 in line:
#            if line == sin2:
                sinner = False
                continue
            else:
                continue
        elif sin1 in line:
#        elif line == sin1:
            sinner = True
            continue
        else:
            yield line

def respacer(script):
    for line in script:
        if line[:4] == "def ":
            yield "".join(["\n\n", line])
        elif line[:2] == "# ":
            yield "".join(["\n", line])
        elif line[:2] == "#!":
            yield "".join([line, "\n"])
        elif len(line) < 2:
            continue
        else:
            yield line


parser = argparse.ArgumentParser(description='Clean up messy marks in python scripts from jupyter and my own notebook dev sections.')
parser.add_argument('script', help='python script to clean')

args = parser.parse_args()

infile = open(args.script, "r")
temp_out = open(".jupyter-cleaner-cleano.temp", "w")

for ln in cleano(infile):
    if len(ln) > 0:
        temp_out.write(ln)

infile.close()
# backup the original script
shutil.move(args.script, args.script + ".bk")
temp_out.close()

temp_in = open(".jupyter-cleaner-cleano.temp", "r")
temp_out = open(".jupyter-cleaner.temp", "w")
for ln in respacer(temp_in):
    if len(ln) > 0:
        temp_out.write(ln)
temp_out.close()
temp_in.close()


outfile = open(".jupyter-cleaner-pen.temp", "w")
temp_in = open(".jupyter-cleaner.temp", "r")
for ln in nodev(temp_in):
    if len(ln) > 0:
        outfile.write(ln)
outfile.close()
temp_in.close()

# discard text above shebang line
outfile = open(args.script, "w")
with open(".jupyter-cleaner-pen.temp", "r") as temp_in:
    flag = False
    for line in temp_in:
        if flag == True:
            outfile.write(line)
        elif line[:2] == "#!":
            flag = True
            outfile.write(line)

# delete intermediates
for rec in glob.glob(".jupyter-cleaner*.temp"):
    os.remove(rec)

### EOF ###
