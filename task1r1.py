inputfile = open('index7.R1.fastq', 'r')
goodoutput = open('index7.good.R1.fastq', 'w')
badoutput = open('index7.bad.R1.fastq', 'w')

primer = 'GAGCCACCGCGC'

def hamming (x1, x2):
	j = 0
	for i in range(len(x1)):
		if x1[i] != x2[i]:
			j += 1
	return (j)

seqobj = []
notshift = []
shift = []
for i, line in enumerate(inputfile):
	line = line.rstrip('\n')
	seqobj.append(line)
	if (i+1) % 4 == 0:
		seq1 = seqobj[1][0:len(primer)]
		mist1 = hamming(primer, seq1)
		seq2 = seqobj[1][1:(len(primer)+1)]
		mist2 = hamming(primer, seq2)
		mist = min(mist1, mist2)
		if mist <= 1:
			for item in seqobj:
				goodoutput.write("%s\n" % item)
		else:
			reason = "h" + str(mist)
			seqobj[0] = seqobj[0] + " " + reason
			for item in seqobj:
				badoutput.write("%s\n" % item)
		seqobj = []

inputfile.close()
goodoutput.close()
badoutput.close()
