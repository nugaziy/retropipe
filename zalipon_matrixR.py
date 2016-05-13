from Bio import SeqIO
import sys, os, re
import numpy as np
from os import listdir
from os.path import isfile, join
from utils import *
import multiprocessing as mp
import pandas as pd
import subprocess
from simplesam import Reader, Writer

def main(inputtable, primer, main_flank_len, refway, outputdir, outputtable, bwaway):
	outputdir += "/"

	if not os.path.exists(outputdir):
		os.makedirs(outputdir)

	megatable = pd.read_table(inputtable)

	inputtablename, ext = os.path.splitext(os.path.basename(inputtable))
	outputtablename, ext = os.path.splitext(outputtable)

	tablefq = open(outputdir + outputtablename + '_zalipon.fq', 'w')

	for index, row in log_progress(megatable.iterrows(), name = 'Create fastq from megatable: ' + inputtablename,
	 every = 1, size = len(megatable)):
		seq = primer + str(row['ALU_BEST']) + str(row['READ1_BEST'])[0 : main_flank_len]
		qual = list(np.repeat('O', len(seq), axis=0))
		tablefq.write('@' + str(row['READNAME']) + '__' + str(int(row['MEGACLUSTER_ID'])) + '\n' + 
			seq + '\n' + 
			'+' + '\n' + 
			''.join(x for x in qual) + '\n')

	tablefq.close()

	logfile = open(outputdir + 'bwa_zalipon_output.log', 'w')
	logfile.write("##### BWA ALN #####")

	bwaaln = bwaway + ' aln ' + refway + ' ' + outputdir + outputtablename + '_zalipon.fq' + ' > ' + outputdir + outputtablename + '.sai'
	print (bwaaln)
	p_aln = subprocess.Popen(bwaaln, stderr=subprocess.PIPE, shell = True)
	logline = p_aln.stderr.read().decode()
	print (logline)
	logfile.write(logline)
	logfile.write("##### BWA SAMSE #####")
	bwasamse = bwaway + ' samse ' + refway + ' ' + outputdir + outputtablename + '.sai' + ' ' + outputdir + outputtablename + '_zalipon.fq' + ' > ' + outputdir + outputtablename + '.sam'
	print(bwasamse)
	p_samse = subprocess.Popen(bwasamse, stderr=subprocess.PIPE, shell = True)
	logline = p_samse.stderr.read().decode()
	print(logline)
	logfile.write(logline)
	logfile.close()

	in_sam_filename = outputdir + outputtablename + '.sam'
	in_sam_file = open(in_sam_filename, 'r')
	tablefile = open(outputdir + outputtablename + '_zalipontable.txt', 'w')
	in_sam = Reader(in_sam_file)
	tablefile.write('MEGACLUSTER_ID\tMDFLAG\tMATCH\n')

	for read in log_progress(in_sam, name = 'Create ZaliponTable from samfile: ' + inputtablename,
	 size = count_fastq_records(inputtable) * 4, every = 1):
		megacluster_id = read.qname.split('__')[4]
		if read.mapped:
			for x in read._tags:
				if re.search('MD', x):
					mdflag = x
				else:
					mdflag = '*'
			match = re.findall(r'\d+', mdflag)
			match = sum([ int(x) for x in match])
		else:
			mdflag = '*'
			match = 0
		tablefile.write(str(megacluster_id) + '\t' + mdflag + '\t' + str(match) + '\n')
	in_sam_file.close()
	tablefile.close()

	megatable_add = pd.read_table(outputdir + outputtablename + '_zalipontable.txt')
	mdcolumn = ['*' for x in list(megatable['MEGACLUSTER_ID'])]
	for index, row in log_progress(megatable.iterrows(), name = 'Get output: ' + inputtablename + ' to ' + outputtablename,
	 every = 1, size = len(megatable)):
		megacluster_id = row['MEGACLUSTER_ID']
		minitable = megatable_add[megatable_add['MEGACLUSTER_ID'] == megacluster_id]
		if len(minitable) == 1:
			mdcolumn[int(index)] = str(minitable['MDFLAG'].values[0])
		else:
			mdmatch = list(minitable['MATCH'])
			mdmatch_index = mdmatch.index(max(mdmatch))
			mdcolumn[int(index)] = str(minitable['MDFLAG'].values[mdmatch_index])

	mdcolumn = {'MDFLAG_MATRIX' : mdcolumn}
	mdcolumn_df = pd.DataFrame(mdcolumn)
	megatable = megatable.join(mdcolumn_df)

	endtable = open(outputdir + outputtable, 'w')
	endtable.close()
	megatable.to_csv(outputdir + outputtable, index=None, sep='\t', mode='a')