#!/usr/bin/python

# functions to work with fasta files
# written by Austin Hammond, BCGSC 2015 (except for fasta_iter)

from itertools import groupby

## fasta_iter
#
#	modified from code written by brentp and retrieved from https://www.biostars.org/p/710/ on May 5, 2015
#	given a fasta file. yield tuples of header, sequence
def fasta_iter(fasta_name):
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fasta_name, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

## dna2aa
#
#	convert dna to amino acid using standard codon table
#	translates all seqs as internal ORFs
#	i.e. doesn't care about start codons, just does straight translations
def dna2aa(seq):
	tbl = {"ATT":"I", "ATC":"I", "ATA":"I", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
        "TTA":"L", "TTG":"L", "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V", "TTT":"F",
        "TTC":"F", "ATG":"M", "TGT":"C", "TGC":"C", "GCT":"A", "GCC":"A", "GCA":"A",
        "GCG":"A", "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G", "CCT":"P", "CCC":"P",
        "CCA":"P", "CCG":"P", "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T", "TCT":"S",
        "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S", "TAT":"Y", "TAC":"Y",
        "TGG":"W", "CAA":"Q", "CAG":"Q", "AAT":"N", "AAC":"N", "CAT":"H", "CAC":"H",
        "GAA":"E", "GAG":"E", "GAT":"D", "GAC":"D", "AAA":"K", "AAG":"K", "CGT":"R",
        "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R", "TAA":"-", "TAG":"-",
        "TGA":"-"}
	j = 0 # for iteration
	aa = ''
	seql = len(seq)
	if seql % 3 != 0:
		seq = seq[0:(seql-seql%3)] #trim to multiple of 3, but should be already
	seql = len(seq) #recalculate
	for j in range(0,seql,3):
		j+=1
		if j == (seql-1):
			codon = seq[j:]
			try:
				aa = aa + tbl[codon]
			except:
				aa = aa + "X"
		else:
			codon = seq[j-1:j+2]
			try:
				aa = aa + tbl[codon]
			except:
				aa = aa+"X"
#	print codon
	return aa

## revcomp
#
#	reverse-complement a dna sequence
def revcomp(seq):
	tbl = {"T":"A","C":"G","G":"C","A":"T","U":"A","R":"Y","Y":"R","S":"S","W":"W","K":"M",
        "M":"K","B":"V","V":"B","D":"H","H":"D","N":"N","-":"-","*":"*"}
	fwd = ''
	for i in seq:
		fwd += tbl[i]
	return fwd[::-1]

## sixframe
#
#	perform naive 6-frame translation of an input DNA sequence
#	i.e. translate through stop codons, no alternate starts
def sixframe(seq):
	one = dna2aa(seq)
	two = dna2aa(seq[1:])
	three = dna2aa(seq[2:])
	none = dna2aa(revcomp(seq))
	ntwo = dna2aa(revcomp(seq)[1:])
	nthree = dna2aa(revcomp(seq)[2:])

	return [one, two, three, none, ntwo, nthree]	

## trimpep
#
#	trim a peptide sequence to the first N-terminal M
def trimpep(seq):
	pos = 0 # for iteration
	tseq = ''
	for res in seq:
		pos += 1
		if res == 'M':
			tseq = seq[pos-1:]
			break

	if len(tseq) > 1:
		return tseq
	else:
		return seq

### EOF ###
