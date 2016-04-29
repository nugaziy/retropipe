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

def main(inputtable, referencefile, primer, window, outputdir, outputtable):
	outputdir += "/"

	if not os.path.exists(outputdir):
		os.makedirs(outputdir)

	megatable = pd.read_table(inputtable)
	reference = pysam.Fastafile(referencefile)

	table = open(outputdir + outputtable, 'w')
	table.write('\t'.join(list(megatable.columns.values)) + '\t' + 'P_MASK' + '\t' + 'P_SHIFT' + '\t' + 'P_HAMMING' + '\n')

	k = 0
	for index, row in log_progress(megatable.iterrows(), name = inputtable, every = 1, size = len(megatable)):
		k += 1
		pos = int(row['POS'])
		if row['STRAND'] == '+':
			start = pos + 1
			end = pos + window
			seq = Seq(reference.fetch(row['CHR'], start, end))
			seq = str(seq).upper()
		else:
			start = pos - window
			end = pos - 1
			seq = reference.fetch(row['CHR'], start, end)
			seq = seq.upper()
		# match = 2, mismatch = -0.5, open_gap = -1, close_gap = abs(-0.5)
		'''
		aligments = pairwise2.align.globalms(seq, primer, 2, 0.5, -100, -100)
		if len(aligments) != 0:
			best_score = max(aligments)
			seq_new = list(best_score[0])
			primer_align = list(best_score[1])
			print(k)
			print(''.join(x for x in primer_align))
			print(''.join(x for x in seq_new))
			p_start = best_score[3]
			p_end = best_score[4]
			mask = ["'"]
			mismatch = 0
			insertion = 0
			deletion = 0
			for i in range(p_start, p_end):
				if seq_new[i] == '-':
					mask.append('^')
					insertion += 1
				elif primer_align[i] == '-':
					mask.append('*')
					deletion += 1
				elif seq_new[i] != primer_align[i]:
					mask.append(0)
					mismatch += 1
				else:
					mask.append(1)
				ham = mismatch + insertion + deletion
			shift = p_start
		else:
			mask = np.repeat('0', len(primer) + 1, axis=0)
			mask[0] = "'"
			mismatch = len(primer)
			insertion = 0
			deletion = 0
			ham = mismatch + insertion + deletion
			shift = 22
		table.write('\t'.join(str(x) for x in list(row)) + '\t' + ''.join(str(x) for x in mask) + '\t' + str(shift) + '\t' + 
	str(mismatch) + '\t' + str(insertion) + '\t' + str(deletion) + '\t' + str(ham) + '\n')
		'''
		score = []
		for i in range(len(seq) - len(primer)):
			seq_new = seq[i : len(primer) + i]
			score.append(hamming(seq_new, primer))
		ham = min(score)
		shift = score.index(min(score))
		mask = []
		mask.append("'")
		seq_new = seq[shift : len(primer) + shift]
		for i in range(len(primer)):
			if seq_new[i] == primer[i]:
				mask.append('1')
			else:
				mask.append('0')
		table.write('\t'.join(str(x) for x in list(row)) + '\t' + ''.join(x for x in mask) + '\t' + 
			str(shift) + '\t' + str(ham) + '\n')
	table.close()