import pandas as pd
import intervaltree as it
import sys, os, re
from os import listdir
from os.path import isfile, join
from utils import *
import pysam
from Bio import pairwise2
from Bio.Seq import Seq
import numpy as np

def hamming (x1, x2):
    j = 0
    for i in range(len(x1)):
        if x1[i] != x2[i]:
            j += 1
    return (j)

def main(inputtable, referencefile, window, outputdir, outputtable):
	outputdir += "/"

	if not os.path.exists(outputdir):
		os.makedirs(outputdir)

	megatable = pd.read_table(inputtable)
	reference = pysam.Fastafile(referencefile)

	table = open(outputdir + outputtable, 'w')
	table.write('\t'.join(list(megatable.columns.values)) + '\t' + 'A_MASK' + '\t' + 'A_SHIFT' + '\t' + 'A_HAMMING' + '\n')

	k = 0
	for index, row in log_progress(megatable.iterrows(), name = inputtable, every = 1, size = len(megatable)):
		alu = str(row['ALU_BEST'])
		k += 1
		pos = int(row['POS'])
		if row['STRAND'] == '+':
			start = pos
			end = pos + window
			seq = Seq(reference.fetch(row['CHR'], start, end))
			seq = seq.reverse_complement()
			seq = str(seq).upper()
		else:
			start = pos - window - 1
			end = pos - 1
			seq = reference.fetch(row['CHR'], start, end)
			seq = seq.upper()
		score = []
		for i in range(len(seq) - len(alu) + 1):
			seq_new = seq[i : len(alu) + i]
			score.append(hamming(seq_new, alu))
		ham = min(score)
		shift = score.index(min(score))
		mask = []
		mask.append("'")
		seq_new = seq[shift : len(alu) + shift]
		for i in range(len(alu)):
			if seq_new[i] == alu[i]:
				mask.append('1')
			else:
				mask.append('0')
		table.write('\t'.join(str(x) for x in list(row)) + '\t' + ''.join(x for x in mask) + '\t' + 
			str(shift) + '\t' + str(ham) + '\n')
	table.close()