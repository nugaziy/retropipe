import pandas as pd
import intervaltree as it
import sys, os, re
from os import listdir
from os.path import isfile, join
from utils import *

def getOverlap(a, b):
	return (max(0, min(a[1], b[1]) - max(a[0], b[0])))

def main(inputtable, inputlibraryfile, outputdir, outputtable):
	outputdir += "/"

	if not os.path.exists(outputdir):
		os.makedirs(outputdir)

	megatable = pd.read_table(inputtable)

	megatable_group = megatable.groupby(['CHR', 'STRAND'])
	tree_dict = {}
	for name, group in  log_progress(megatable_group, name = 'Create tree_list with', every = 1, who = 'classes: chrom & strand'):
		readwindow_list = [int(x) for x in list(group['TLEN'])]
		if name[1] == '+':
			start_group = [pos - read_w for read_w, pos in zip(readwindow_list, list(group['POS']))]
			end_group = [pos + 1 for ins_w, pos in zip(readwindow_list, list(group['POS']))]
		else:
			start_group = [pos for ins_w, pos in zip(readwindow_list, list(group['POS']))]
			end_group = [pos + read_w + 1 for read_w, pos in zip(readwindow_list, list(group['POS']))]
		tree_dict[name[0] + name[1]] = it.IntervalTree(it.Interval(start, end, data)
		 for start, end, data in zip(start_group, end_group, list(group['MEGACLUSTER_ID'])))

	megatable_id = megatable['MEGACLUSTER_ID']
	megatable = megatable.set_index(megatable_id)

	inputfile, ext = os.path.splitext(inputlibraryfile)
	megatable[inputfile] = pd.Series(['Unknown' for x in range(len(megatable))], index = megatable_id)
	megatable[inputfile + '_DIV'] = pd.Series([0 for x in range(len(megatable))], index = megatable_id)
	megatable[inputfile + '_OVERLAP'] = pd.Series([0.0 for x in range(len(megatable))], index = megatable_id)
	megatable[inputfile + '_STRAND'] = pd.Series(['*' for x in range(len(megatable))], index = megatable_id)

	replib = pd.read_table(inputlibraryfile)
	replib_group = replib.groupby(['CHR'])
	for name, group in log_progress(replib_group, name = 'For lib - ' + inputfile, every = 1, who = 'classes: chrom & strand'):
		for strand_i in ['+', '-']:
			if name + strand_i in tree_dict:
				rep_int = zip(list(group['START']), list(group['END']), zip(list(group['NAME']), list(group['DIV']), list(group['STRAND'])))
				for begin, end, rep_data in rep_int:
					finter = tree_dict[name + strand_i][begin:end]
					if len(finter) > 0:
						for x in finter:
							overlap = getOverlap([x.begin, x.end], [begin, end]) / (x.end - x.begin)
							if megatable[inputfile][x.data] == 'Unknown':
								megatable.set_value(x.data, inputfile, rep_data[0])
								megatable.set_value(x.data, inputfile + '_DIV', rep_data[1])
								megatable.set_value(x.data, inputfile + '_OVERLAP', overlap)
								megatable.set_value(x.data, inputfile + '_STRAND', rep_data[2])
							else:
								if overlap > megatable[inputfile + '_OVERLAP'][x.data]:
									megatable.set_value(x.data, inputfile, rep_data[0])
									megatable.set_value(x.data, inputfile + '_DIV', rep_data[1])
									megatable.set_value(x.data, inputfile + '_OVERLAP', overlap)
									megatable.set_value(x.data, inputfile + '_STRAND', rep_data[2])
	table = open(outputdir + outputtable, 'w')
	table.close()
	megatable.to_csv(outputdir + outputtable, index=None, sep='\t', mode='a')