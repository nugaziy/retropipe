import pandas as pd
import intervaltree as it
import sys, os, re
from os import listdir
from os.path import isfile, join
from utils import *

def main(inputtable, inputlibrary, outputdir, outputtable, inswindow, readwindow, is_dinamic_window):
	inputlibrary += "/"
	outputdir += "/"

	if not os.path.exists(outputdir):
		os.makedirs(outputdir)

	onlyfiles = [f for f in listdir(inputlibrary) if isfile(join(inputlibrary, f))]
	megatable = pd.read_table(inputtable)
	megatable_id = list(megatable['MEGACLUSTER_ID'])

	megatable_group = megatable.groupby(['CHR', 'STRAND'])
	tree_dict = {}
	for name, group in  log_progress(megatable_group, name = 'Create tree_list with', every = 1, who = 'classes: chrom & strand'):
		if is_dinamic_window:
			readwindow_list = [int(0.9 * x) for x in list(group['TLEN'])]
		else:
			readwindow_list = [readwindow for x in list(group['TLEN'])]
		inswindow_list = [inswindow for x in list(group['TLEN'])]
		if name[1] == '+':
			start_group = [pos - read_w for read_w, pos in zip(readwindow_list, list(group['POS']))]
			end_group = [pos + ins_w for ins_w, pos in zip(inswindow_list, list(group['POS']))]
		else:
			start_group = [pos - ins_w for ins_w, pos in zip(inswindow_list, list(group['POS']))]
			end_group = [pos + read_w for read_w, pos in zip(readwindow_list, list(group['POS']))]
		tree_dict[name[0] + name[1]] = it.IntervalTree(it.Interval(start, end, data)
		 for start, end, data in zip(start_group, end_group, list(group['MEGACLUSTER_ID'])))

	for filename in onlyfiles:
		repcolumn = {x : 'Unknown' for x in megatable_id}
		filename = filename.rstrip()
		inputfile, ext = os.path.splitext(filename)
		replib = pd.read_table(inputlibrary + filename)
		replib_group = replib.groupby(['CHR', 'STRAND'])
		for name, group in log_progress(replib_group, name = 'For lib - ' + inputfile, every = 1, who = 'classes: chrom & strand'):
			if name[0] + name[1] in tree_dict:
				if str(name[1]) == '+':
					point = zip(list(group['START']), list(group['NAME']))
				else:
					point = zip(list(group['END']), list(group['NAME']))
				for pos, repname in point:
					finter = tree_dict[name[0] + name[1]][pos]
					if len(finter) > 0:
						for x in finter:
							repcolumn[x.data] = repname
		repcolumn = {inputfile : list(repcolumn.values())}
		repcolumn_df = pd.DataFrame(repcolumn)
		megatable = megatable.join(repcolumn_df)
	table = open(outputdir + outputtable, 'w')
	table.close()
	megatable.to_csv(outputdir + outputtable, index=None, sep='\t', mode='a')