from Bio import SeqIO
import sys, os, re
import numpy as np
from os import listdir
from os.path import isfile, join
from utils import *
import multiprocessing as mp
import pandas as pd

def recovery(df, inputdir, filename1, filename2, name):
	df_index = df.index
	df = df.set_index(df['READNAME'])
	R1_reads = SeqIO.parse(inputdir + filename1, "fastq")
	R2_reads = SeqIO.parse(inputdir + filename2, "fastq")
	for r1, r2 in log_progress(zip(R1_reads, R2_reads),
     name = name, size = count_fastq_records(inputdir + filename1), every = 250):
		if r1.id in df['READNAME']:
			df.set_value(r1.id, 'READ1', r1.seq)
			df.set_value(r1.id, 'READ2', r2.seq)
	df = df.set_index(df_index)
	return (df)

def main(inputtable, inputdir, outputdir):
	inputdir += "/"
	outputdir += "/"

	# Read files in folder
	onlyfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]

	r1_files = {}
	r2_files = {}

	for filename in onlyfiles:
		filename = filename.rstrip()
		if re.search('good', filename):
			if re.search('R1', filename):
				key_filename = filename.split('R1')[0]
				r1_files[key_filename] = filename
			elif re.search('R2', filename):
				key_filename = filename.split('R2')[0]
				r2_files[key_filename] = filename

	conform_files = []
	nonconform_files = []

	for key in r1_files:
		if key in r2_files:
			conform_files.append([r1_files[key], r2_files[key]])
			del r2_files[key]
		else: nonconform_files.append(r1_files[key])

	nonconform_files = nonconform_files + list(r2_files.values())

	if not os.path.exists(outputdir):
		os.makedirs(outputdir)

	if len(nonconform_files) != 0:
		print ('I can\'t read this files' + str(nonconform_files))

	filedict = {}

	for filename1, filename2 in conform_files:
		readsname = filename1.split('R1')[0]
		readsname = readsname.rsplit('.', 1)[0]
		filedict[readsname] = [filename1, filename2]

	data = pd.read_table(inputtable)
	data_new = pd.DataFrame(columns = list(data.columns.values))
	data_group = data.groupby(['FILENAME'])
	for name, group in data_group:
		if name in list(filedict.keys()):
			file1 = list(filedict[name])[0]
			file2 = list(filedict[name])[1]
			df_new = recovery(group, inputdir, file1, file2, name)
			print(list(df_new.columns.values))
			data_new = data_new.append(df_new)


	table = open(outputdir + 'megatable_fq.txt', 'w')
	table.close()
	data_new.to_csv(outputdir + 'megatable_fq.txt', index=None, columns = list(data.columns.values), sep='\t', mode='a')
