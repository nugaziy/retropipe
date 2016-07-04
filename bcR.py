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
		return(0)
	else:
		return (1)

def which (rowsum):
	for x in rowsum:
		if x > 2:
			return(1)
	return(0)

def main(inputtable):

	megatable = pd.read_table(inputtable)

	k = 0
	k_1 = 0
	for index, row in log_progress(megatable.iterrows(), name = inputtable, every = 250):
		k += 1
		barcode = row['BARCODE_LIST'].split(',')
		hamming_table = pd.DataFrame(columns = barcode, index = barcode)
		for i in range(len(barcode)):
			for j in  range(len(barcode)):
				hamming_table.set_value(barcode[i], barcode[j], hamming(barcode[i], barcode[j]))
		rowsum = hamming_table.sum(axis = 1)
		k_1 += which(rowsum)
	print(k_1 / k)
