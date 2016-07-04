import pandas as pd
import numpy as np
import intervaltree as it
import sys, os, re
from os import listdir
from os.path import isfile, join
from utils import *

def main(inputtable, inputalutable, outputdir, outputtable):
	outputdir += "/"

	if not os.path.exists(outputdir):
		os.makedirs(outputdir)

	alutable = pd.read_table(inputalutable)

	alutable_group = alutable.groupby(['CHR', 'STRAND'])
	tree_dict = {}
	for name, group in  log_progress(alutable_group, name = 'Create tree_list with', every = 1, who = 'classes: chrom & strand'):
		start_group = list(group['START'])
		end_group = list(group['END'])
		tree_dict[name[0] + name[1]] = it.IntervalTree(it.Interval(start, end, data)
		 for start, end, data in zip(start_group, end_group, list(group['NAME'])))

	megatable = pd.read_table(inputtable)

	megatable['RangeAluAGCT_NameAlu'] = pd.Series(['*' for x in range(len(megatable))], index = megatable.index)
	megatable['RangeAluAGCT_DistAGCT'] = pd.Series(['*' for x in range(len(megatable))], index = megatable.index)

	megatable_group = megatable.groupby(['CHR', 'STRAND'])
	for name, group in log_progress(replib_group, name = 'For table - ' + inputtable, every = 1, who = 'classes: chrom & strand'):
		if name[0] + name[1] in tree_dict:
			point = zip(list(group['POS']), list(megatable.index))
			for pos, idx in point:
				finter = tree_dict[name[0] + name[1]][pos]
				if len(finter) > 0:
					for x in finter:
						megatable.set_value(idx, inputfile, x.data)

	megatable.to_csv(outputdir + outputtable, index=None, sep='\t')