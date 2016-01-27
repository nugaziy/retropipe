#!usr/bin/python3

from Bio import SeqIO

def hamming (x1, x2):
	j = 0
	for i in range(len(x1)):
		if x1[i] != x2[i]:
			j += 1
	return (j)

def trim_primers(records, primer):
	len_primer = len(primer)
	for record in records:
		min_ham = min (hamming (primer, record[0:len_primer]), hamming (primer, record[1:(len_primer+1)]))
		if min_ham <= 1:
			yield record


def trim_ads(records, ad1, ad2, barlen):
	len_ad1 = len(ad1)
	len_ad2 = len(ad2)
	for record in records:
		seq1 = record[0:len_ad1]
		seq1_shift = record[1:(len_ad1+1)]
		seq2 = record[(len_ad1+barlen):(len_ad1+barlen+len_ad2)]
		seq2_shift = record[(len_ad1+barlen+1):(len_ad1+barlen+len_ad2+1)]
		min_ham = min (hamming (ad1+ad2, seq1+seq2), hamming (ad1+ad2, seq1_shift+seq2_shift))
		if min_ham <= 1:
			yield record


primer = 'GAGCCACCGCGC'
ad1 = 'GCGTGCTGCGG'
ad2 = 'AGGGCGGT'
barlen = 9

original_R1_reads = SeqIO.parse("index7.R1.fastq", "fastq")
original_R2_reads = SeqIO.parse("index7.R2.fastq", "fastq")
trimmed_R1_reads = trim_primers(original_R1_reads, primer)
trimmed_R2_reads = trim_ads(original_R2_reads, ad1, ad2, barlen)
count_R1 = SeqIO.write(trimmed_R1_reads, "trimmed_R1.fastq", "fastq")
count_R2 = SeqIO.write(trimmed_R2_reads, "trimmed_R2.fastq", "fastq")
