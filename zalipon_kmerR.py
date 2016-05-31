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
from collections import Counter

def hamming (x1, x2):
    j = 0
    for i in range(len(x1)):
        if x1[i] != x2[i]:
            j += 1
    return (j)

def main(inputtable, referencefile, window, k_min, k_max, outputdir, outputtable):
	outputdir += "/"

	if not os.path.exists(outputdir):
		os.makedirs(outputdir)

	megatable = pd.read_table(inputtable)
	reference = pysam.Fastafile(referencefile)

	table = open(outputdir + outputtable, 'w')
	table.write('KMER' + '\t' + 'AMOUNT' + '\n')

	kmer = []
	for index, row in log_progress(megatable.iterrows(), name = inputtable, every = 250, size = len(megatable)):
		if (str(row['Alu_hg38']) == 'Unknown') and (str(row['Alu_dbRIP_hg38']) == 'Unknown'):
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
			for i in range(k_min, k_max + 1):
				for j in range(len(seq) - i):
						kmer.append(seq[j : j + i])
	kmer_count = dict(Counter(kmer))
	for key, value in kmer_count.items():
		table.write(str(key) + '\t' + str(value) + '\n')
	table.close()