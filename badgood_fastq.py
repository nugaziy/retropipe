#!usr/bin/python3

from Bio import SeqIO
import sys, argparse, os
import numpy as np

parser = argparse.ArgumentParser(description='Trim fastq reads (r1, r2) by primer and adapters')
parser.add_argument('-i', '--input', help='input files', nargs='+')
parser.add_argument('-o', '--output', help='output directory')
parser.add_argument('-m','--mistake', type = int, help='amount of mistake')
parser.add_argument('-p','--primer', help='primer')
parser.add_argument('-a','--ads', help='adapter and ad.green', nargs='+')
parser.add_argument('-b','--barlen', type = int, help='barcode length')
args = parser.parse_args()

dir1, filename1 = os.path.split(args.input[0])
dir2, filename2 = os.path.split(args.input[1])

outputfile1, outputfile1_ext = os.path.splitext(filename1)
outputfile2, outputfile2_ext = os.path.splitext(filename2)

mist = args.mistake
primer = args.primer
ad1 = args.ads[0]
ad2 = args.ads[1]
barlen = args.barlen


goodr1 = open(args.output + outputfile1 + '_good.fastq', 'w')
goodr2 = open(args.output + outputfile2 + '_good.fastq', 'w')
badr1 = open(args.output + outputfile1 + '_bad.fastq', 'w')
badr2 = open(args.output + outputfile2 + '_bad.fastq', 'w')




def hamming (x1, x2, m):
	j = 0
	for i in range(len(x1)):
		if x1[i] != x2[i]:
			j += 1
			if j > m:
				return (False)
	return (True)

def trim_primers(record, primer, m):
	badgood = {"good":False, "bad":np.array([0, 0, 0])}
	len_primer = len(primer)
	ham_prim = (hamming(primer, record[0:len_primer], m))or(hamming(primer, record[1:(len_primer+1)], m))
	if ham_prim:
		badgood['good'] = True
	else:
		badgood['bad'] = np.array([1, 0, 0])
	return badgood

def trim_ads(record, ad1, ad2, barlen, m):
	badgood = {"good":False, "bad":np.array([0, 0, 0])}
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
	else:
		if not((ham_ad1)or(ham_ad2)):
			badgood['bad'] = np.array([0, 1, 1])
		elif not(ham_ad1):
			badgood['bad'] = np.array([0, 1, 0])
		else:
			badgood['bad'] = np.array([0, 0, 1])
	return badgood



original_R1_reads = SeqIO.parse(args.input[0], "fastq")
original_R2_reads = SeqIO.parse(args.input[1], "fastq")

count = np.array([0, 0, 0])

for zipi in zip(original_R1_reads, original_R2_reads):
	r1,r2 = zipi
	fr1 = trim_primers(r1, primer, mist)
	if fr1['good']:
		fr2 = trim_ads(r2, ad1, ad2, barlen, mist)
		if fr2['good']:
			goodr1.write(str(r1.format('fastq')))
			goodr2.write(str(r2.format('fastq')))
		else:
			badread = str(r1.format('fastq'))
			badread = badread.split('\n')
			badread[0] = badread[0] + ' m:' + ''.join(np.char.mod('%d', fr2['bad']))
			badread = '\n'.join(badread)
			badr1.write(badread)
			badr2.write(str(r2.format('fastq')))
			count = np.sum([count, fr2['bad']], axis=0)
	else:
		badread = str(r1.format('fastq'))
		badread = badread.split('\n')
		badread[0] = badread[0] + ' m:' + ''.join(np.char.mod('%d', fr1['bad']))
		badread = '\n'.join(badread)
		badr1.write(badread)
		badr2.write(str(r2.format('fastq')))
		count = np.sum([count, fr1['bad']], axis=0)


goodr1.close()
goodr2.close()
badr1.close()
badr2.close()

print (count)
