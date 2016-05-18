import pandas as pd
import re
from collections import defaultdict
from collections import Counter
import sys, os, re
from utils import *
from os import listdir
from os.path import isfile, join
from scipy.stats import gaussian_kde
import numpy as np
import networkx as nx

def hamming (x1, x2):
	j = 0
	for i in range(len(x1)):
		if x1[i] != x2[i]:
			j += 1
			if j > 1:
				return (0)
	if j == 1:
		return (1)
	return (0)

def main(inputtable_pcread, inputtable_for_output, outputdir, outputtable):
	outputdir += '/'

	megatable = pd.read_table(inputtable_pcread)
	megatable_colnames = list(megatable.columns.values)
	barcode_col = [x.split('_BARCODE_LIST')[0] for x in megatable_colnames if re.search('BARCODE_LIST', x)]
	scc_column = {}
	mask_column = {}
	megatable_output = pd.read_table(inputtable_for_output)
	megatable_output_index = megatable_output.index
	megatable_output = megatable_output.set_index(megatable_output.MEGACLUSTER_ID.astype('int'))
	for col in barcode_col:
		megatable_output_colnames = list(megatable_output.columns.values)
		nmbr_scc = megatable_output_colnames.index(col + '_NUM_BARCODES') + 1
		nmbr_mask = megatable_output_colnames.index(col + '_NUM_BARCODES') + 2
		megatable_output.insert(nmbr_scc, col + '_REAL_NUM_BARCODE', ['*' for x in range(len(megatable))])
		megatable_output.insert(nmbr_mask, col + '_BARCODE_MASK', ['*' for x in range(len(megatable))])
	for index, row in log_progress(megatable.iterrows(), name = inputtable_pcread, every = 250, size = len(megatable)):
		if (megatable_output.Alu_hg38[int(row['MEGACLUSTER_ID'])] == 'Unknown') and (megatable_output.Alu_dbRIP_hg38[int(row['MEGACLUSTER_ID'])] == 'Unknown'):
			scc = []
			mask = []
			barcode_cluster = []
			for col in barcode_col:
				barcode_list = re.findall(r'.........', row[col + '_BARCODE_LIST'])
				barcode_list.remove('RRRRRRRRR')
				barcode = list(set(barcode_list))
				barcode_cluster.append(barcode)
				G = nx.Graph()
				for i in range(len(barcode) - 1):
					G.add_node(barcode[i])
					for j in range(i + 1, len(barcode)):
						G.add_node(barcode[j])
						if hamming(barcode[i],  barcode[j]) == 1:
							G.add_edge(barcode[i], barcode[j])
				scc_len = [len(x) for x in list(nx.connected_components(G))]
				scc.append(len(scc_len))
			scc = [str(x) for x in scc]
			for i in range(len(barcode_col)):
				mask_col = []
				for j in range(len(barcode_col)):
					mask_col.append(len(list(set(barcode_cluster[i]).intersection(barcode_cluster[j]))))
				mask_col = ','.join(str(x) for x in mask_col)
				mask.append(mask_col)
			for i in range(len(barcode_col)):
				megatable_output.set_value(int(row['MEGACLUSTER_ID']), barcode_col[i] + '_REAL_NUM_BARCODE', scc[i])
				megatable_output.set_value(int(row['MEGACLUSTER_ID']), barcode_col[i] + '_BARCODE_MASK', mask[i])
	table = open(outputdir + outputtable, 'w')
	table.close()
	megatable_output.to_csv(outputdir + outputtable, index=None, sep='\t', mode='a')