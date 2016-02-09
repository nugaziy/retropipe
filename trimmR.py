from Bio import SeqIO
import sys, os, re
import numpy as np
from os import listdir
from os.path import isfile, join

# Add functions:


def log_progress(sequence, name, every=None, size=None):
	from ipywidgets import IntProgress, HTML, VBox
	from IPython.display import display
	
	is_iterator = False
	if size is None:
		try:
			size = len(sequence)
		except TypeError:
			is_iterator = True
	if size is not None:
		if every is None:
			if size <= 200:
				every = 1
			else:
				every = size / 200     # every 0.5%
	else:
		assert every is not None, 'sequence is iterator, set every'
	
	if is_iterator:
		progress = IntProgress(min=0, max=1, value=1)
		progress.bar_style = 'info'
	else:
		progress = IntProgress(min=0, max=size, value=0)
	label = HTML()
	box = VBox(children=[label, progress])
	display(box)
	
	index = 0
	try:
		for index, record in enumerate(sequence, 1):
			if index == 1 or index % every == 0:
				if is_iterator:
					label.value = name + ': ' + '{index} reads / ?'.format(index=index)
				else:
					progress.value = index
					label.value = u'{index} / {size}'.format(
					index=index,
						size=size
					)
			yield record
	except:
		progress.bar_style = 'danger'
		raise
	else:
		progress.bar_style = 'success'
		progress.value = index
		label.value = str(name + ': ' + str(index) + ' reads' or '?')



def hamming (x1, x2, m):
	j = 0
	for i in range(len(x1)):
		if x1[i] != x2[i]:
			j += 1
			if j > m:
				return (False)
	return (True)



def trim_primers(record, primer, m):
	record = str(record.format('fastq'))
	record = record.split('\n')
	record = record[1]
	badgood = {'good':False, 'bad':np.array([0, 0, 0])}
	len_primer = len(primer)
	ham_prim = (hamming(primer, record[0:len_primer], m))or(hamming(primer, record[1:(len_primer+1)], m))
	if ham_prim:
		badgood['good'] = True
	else:
		badgood['bad'] = np.array([1, 0, 0])
	return badgood



def trim_ads(record, ad1, ad2, barlen, m):
	record = str(record.format('fastq'))
	record = record.split('\n')
	record = record[1]
	badgood = {'good':False, 'bad':np.array([0, 0, 0]), 'barcode':None}
	len_ad1 = len(ad1)
	len_ad2 = len(ad2)
	seq1 = record[0:len_ad1]
	seq1_shift = record[1:(len_ad1+1)]
	seq2 = record[(len_ad1+barlen):(len_ad1+barlen+len_ad2)]
	seq2_shift = record[(len_ad1+barlen+1):(len_ad1+barlen+len_ad2+1)]
	ham_ad1 = (hamming(ad1, seq1, m))or(hamming(ad1, seq1_shift, m))
	ham_ad2 = (hamming(ad2, seq2, m))or(hamming(ad1, seq2_shift, m))
	if (ham_ad1)and(ham_ad2):
		badgood['good'] = True
		if hamming(ad2, seq2, m) > hamming(ad1, seq2_shift, m):
			badgood['barcode'] = record[len_ad1:(len_ad1+barlen)]
		else: badgood['barcode'] = record[(len_ad1+1):(len_ad1+barlen+1)]
	else:
		if not((ham_ad1)or(ham_ad2)):
			badgood['bad'] = np.array([0, 1, 1])
		elif not(ham_ad1):
			badgood['bad'] = np.array([0, 1, 0])
		else:
			badgood['bad'] = np.array([0, 0, 1])
	return badgood



def concate(x1, x2):
	result = ''
	for i in range(len(x1)):
		result = result + str(x1[i]) + '-' + str(x2[i])
		if (i != (len(x1) - 1)): result = result + ','
	return (result)



def trim_reads(filename1, filename2, inputdir, outputdir, mist, primer, ad1, ad2, barlen):
	readsname = filename1.split('R1')[0]
	readsname = readsname.rsplit('.', 1)[0]
	
	outputfile1, ext = os.path.splitext(filename1)
	outputfile2, ext = os.path.splitext(filename2)
	
	goodr1 = open(outputdir + outputfile1 + '_good.fastq', 'w')
	goodr2 = open(outputdir + outputfile2 + '_good.fastq', 'w')
	badr1 = open(outputdir + outputfile1 + '_bad.fastq', 'w')
	badr2 = open(outputdir + outputfile2 + '_bad.fastq', 'w')
	
	original_R1_reads = SeqIO.parse(inputdir + filename1, "fastq")
	original_R2_reads = SeqIO.parse(inputdir + filename2, "fastq")
	
	count = np.array([0, 0, 0])
	elem = ('primer', 'ad', 'green')
	
	for zipi in log_progress(zip(original_R1_reads, original_R2_reads), name = readsname, every = 1):
		r1,r2 = zipi
		fr1 = trim_primers(r1, primer, mist)
		if fr1['good']:
			fr2 = trim_ads(r2, ad1, ad2, barlen, mist)
			if fr2['good']:
				goodread = str(r2.format('fastq'))
				goodread = goodread.split('\n')
				goodread[0] = goodread[0] + ' barcode:' + str(fr2['barcode'])
				goodread = '\n'.join(goodread)
				goodr1.write(str(r1.format('fastq')))
				goodr2.write(goodread)
			else:
				badread = str(r1.format('fastq'))
				badread = badread.split('\n')
				badread[0] = badread[0] + ' reason:' + concate(elem, (np.char.mod('%d', fr2['bad'])))
				badread = '\n'.join(badread)
				badr1.write(badread)
				badr2.write(str(r2.format('fastq')))
				count = np.sum([count, fr2['bad']], axis=0)
		else:
			badread = str(r1.format('fastq'))
			badread = badread.split('\n')
			badread[0] = badread[0] + ' reason:' + concate(elem, (np.char.mod('%d', fr1['bad'])))
			badread = '\n'.join(badread)
			badr1.write(badread)
			badr2.write(str(r2.format('fastq')))
			count = np.sum([count, fr1['bad']], axis=0)
	
	goodr1.close()
	goodr2.close()
	badr1.close()
	badr2.close()
	
	count_elem = concate(elem, (np.char.mod('%d', count)))
	print ('Mistake: place-amount: ' + count_elem + ' for ' + readsname + '\n')



def main(inputdir, outputdir, mist, primer, ad1, ad2, barlen):
	# Read files in folder
	onlyfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]
	
	r1_files = {}
	r2_files = {}
	
	for filename in onlyfiles:
		filename = filename.rstrip()
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
	
	for pair in conform_files:
		filename1 = str(pair[0])
		filename2 = str(pair[1])
		responses = trim_reads(filename1, filename2,
                                            inputdir, outputdir, mist,
                                            primer, ad1, ad2, barlen)
	
	if len(nonconform_files) != 0:
		print ('I can\'t read this files' + str(nonconform_files))


