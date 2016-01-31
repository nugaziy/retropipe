#!usr/bin/python3

from Bio import SeqIO
import sys, argparse, os

parser = argparse.ArgumentParser(description='inputfiles and outputdir')
parser.add_argument('-i','--input', help='inputfiles', nargs='+')
parser.add_argument('-o','--output', help='outputdir', required=True)
parser.add_argument('-m','--mistake', type = int, help='amount of mistake', required=True)
parser.add_argument('-p','--primer', help='primer', required=True)
parser.add_argument('-a','--ads', help='adapters', nargs='+')
parser.add_argument('-b','--barlen', type = int, help='barlen', required=True)
args = parser.parse_args()

dir1, filename1 = os.path.split(args.input[0])
dir2, filename2 = os.path.split(args.input[1])

outputfile1, outputfile1_ext = os.path.splitext(filename1)
outputfile2, outputfile2_ext = os.path.splitext(filename2)


goodr1 = open(args.output + outputfile1 + '_good.fastq', 'w')
goodr2 = open(args.output + outputfile2 + '_good.fastq', 'w')
badr1 = open(args.output + outputfile1 + '_bad.fastq', 'w')
badr2 = open(args.output + outputfile2 + '_bad.fastq', 'w')


mist = args.mistake
primer = args.primer
ad1 = args.ads[0]
ad2 = args.ads[1]
barlen = args.barlen

def hamming (x1, x2, m):
	j = 0
	for i in range(len(x1)):
		if x1[i] != x2[i]:
			j += 1
			if j > (m + 1):
				return (j)
	return (j)

def trim_primers(record, primer, m):
	badgood = {"good":None, "bad":None, "reason": None}
	len_primer = len(primer)
	min_ham = min (hamming (primer, record[0:len_primer], m), hamming (primer, record[1:(len_primer+1)], m))
	if min_ham <= m:
		badgood['good'] = record
	else:
		badgood['bad'] = record
		badgood['reason'] = min_ham
	return badgood

def trim_ads(record, ad1, ad2, barlen, m):
	badgood = {"good":None, "bad":None, "reason": None}
	len_ad1 = len(ad1)
	len_ad2 = len(ad2)
	seq1 = record[0:len_ad1]
	seq1_shift = record[1:(len_ad1+1)]
	seq2 = record[(len_ad1+barlen):(len_ad1+barlen+len_ad2)]
	seq2_shift = record[(len_ad1+barlen+1):(len_ad1+barlen+len_ad2+1)]
	min_ham = min (hamming (ad1+ad2, seq1+seq2, m), hamming (ad1+ad2, seq1_shift+seq2_shift, m))
	if min_ham <= m:
		badgood['good'] = record
	else:
		badgood['bad'] = record
		badgood['reason'] = min_ham
	return badgood


def reverse(record): 
	record['bad'], record['good'] = record['good'], record['bad']
	return record


original_R1_reads = SeqIO.parse(args.input[0], "fastq")
original_R2_reads = SeqIO.parse(args.input[1], "fastq")

for zipi in zip(original_R1_reads, original_R2_reads):
	r1,r2 = zipi
	fr1 = trim_primers(r1, primer, mist)
	fr2 = trim_ads(r2, ad1, ad2, barlen, mist)
	if not((fr1['good'])and(fr2['good'])):
		if fr1['good']:
			fr1 = reverse(fr1)
			fr1['reason'] = str(fr2['reason']) + ':d'
			fr2['reason'] = str(fr2['reason']) + ':d'
		else:
			fr2 = reverse(fr2)
			fr2['reason'] = str(fr1['reason']) + ':d'
			fr1['reason'] = str(fr1['reason']) + ':d'
	if fr1['good']: goodr1.write(str(fr1['good'].format('fastq')))
	else:
		badread = str(fr1['bad'].format('fastq'))
		badread = badread.split('\n')
		badread[0] = badread[0] + ' m:' + str(fr1['reason'])
		badread = '\n'.join(badread)
		badr1.write(badread)
	if fr2['good']: goodr2.write(str(fr2['good'].format('fastq')))
	else:
		badread = str(fr2['bad'].format('fastq'))
		badread = badread.split('\n')
		badread[0] = badread[0] + ' m:' + str(fr2['reason'])
		badread = '\n'.join(badread)
		badr2.write(badread)


goodr1.close()
goodr2.close()
badr1.close()
badr2.close()
