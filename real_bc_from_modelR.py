import pandas as pd
import re
import sys, os, re
from utils import *
from os import listdir
from os.path import isfile, join
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

def main(inputtable_pcread, inputtable_for_output, model_bc_fit, outputdir, outputtable):
	outputdir += '/'

	df = pd.read_table(model_bc_fit)
	x = np.array(list(df['N_OBS']))
	y = np.array(list(df['FREE_POINTS']))
	min_y = [min([int(y) for y in x.split(',')]) for x in list(df['FREE_POINTS_LIST'])]
	max_y = [max([int(y) for y in x.split(',')]) for x in list(df['FREE_POINTS_LIST'])]

	p3_mean = np.poly1d(np.polyfit(x, y, 3))
	p3_min = np.poly1d(np.polyfit(x, min_y, 3))	
	p3_max = np.poly1d(np.polyfit(x, max_y, 3))

	megatable_pcread = pd.read_table(inputtable_pcread)
	megatable_output = pd.read_table(inputtable_for_output)
	megatable_output = megatable_output.set_index(megatable_output.MEGACLUSTER_ID.astype('int'))
	megatable_output_colnames = list(megatable_output.columns.values)
	remove_list = [x for x in megatable_output_colnames if re.search('_NUM_READS', x) or 
	 re.search('_NUM_BARCODES', x) or re.search('_REAL_NUM_BARCODE', x) or not re.search('P125_5_ACAGTC')]
	megatable_AluYa5_colnames = [x for x in megatable_output_colnames if not x in remove_list]
	megatable_AluYa5 = megatable_output.loc[megatable_output['Alu_hg38'] == 'Alu_Ya5' and megatable_output['P125_5_ACAGTC_L007__NUM_BARCODES'] > 49 and 
	 megatable_output['P125_5_ACAGTC_L007__NUM_BARCODES'] < 10001, megatable_AluYa5_colnames]

	megatable_AluYa5_output = megatable_AluYa5

	nmbr_scc = megatable_output_colnames.index('P125_5_ACAGTC_L007__NUM_BARCODES') + 1
	newcol = ['_REAL_BARCODES_HAMMING', '_REAL_BARCODES_MODEL_MEAN', '_REAL_BARCODES_MODEL_MIN', '_REAL_BARCODES_MODEL_MAX']
	for i in range(len(newcol)):
		megatable_AluYa5_output.insert(nmbr_scc + i, 'P125_5_ACAGTC_L007_' + newcol[i], [0.0 for x in range(len(megatable_AluYa5_output))])
	
	for index, row in log_progress(megatable_AluYa5.iterrows(), name = 'Alu_Ya5', every = 1, size = len(megatable_AluYa5)):
		barcode_list = re.findall(r'.........', str(megatable_pcread[megatable_pcread['MEGACLUSTER_ID'] == index]['P125_5_ACAGTC_L007__BARCODE_LIST']))
		barcode_list.remove('RRRRRRRRR')
		barcode = list(set(barcode_list))
		G = nx.Graph()
		for i in range(len(barcode) - 1):
			G.add_node(barcode[i])
			for j in range(i + 1, len(barcode)):
				G.add_node(barcode[j])
				if hamming(barcode[i],  barcode[j]) == 1:
					G.add_edge(barcode[i], barcode[j])
		scc = len(list(nx.connected_components(G)))
		points = len([0 for x in list(nx.degree(G).values()) if x == 0])
		output_value = [scc, p3_mean(points), p3_min(points), p3_max(points)]
		for i in range(len(newcol)):
			megatable_AluYa5_output.set_value(index, 'P125_5_ACAGTC_L007_' + newcol[i], output_value[i])
	table = open(outputdir + outputtable, 'w')
	table.close()
	megatable_output.to_csv(outputdir + outputtable, index=None, sep='\t', mode='a')