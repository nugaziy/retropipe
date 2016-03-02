from Bio import SeqIO
import sys, os, re
import numpy as np
from os import listdir
from os.path import isfile, join
from utils import *
from collections import namedtuple
info = namedtuple('info', 'good read errors')

def hamming (x1, x2, m):
    j = 0
    for i in range(len(x1)):
        if x1[i] != x2[i]:
            j += 1
            if j > m:
                return (False)
    return (True)



def trim_primers(record, primer, m, elem_remove):
    record_seq = record.seq
    len_primer = len(primer)
    ham_prim = hamming(primer, record_seq[0:len_primer], m)
    ham_prim_shift = hamming(primer, record_seq[1:(len_primer+1)], m)
    if ham_prim or ham_prim_shift:
        if ham_prim: 
            alu = record_seq[len_primer:(len_primer+6)]
            for elem in elem_remove:
                if record_seq[(len_primer+6):].find(elem, 0) != -1:
                    return info(good = False, read = None,
                            errors = np.array([0, 0, 0, 1]))
            record = record[(len_primer+6):]
            record.description += ' alu.seq:' + str(alu)
            return info(good = True, read = record,
                    errors = np.array([0, 0, 0, 0]))
        else:
            alu = record_seq[(len_primer+1):(len_primer+7)]
            for elem in elem_remove:
                if record_seq[(len_primer+7):].find(elem, 0) != -1:
                    return info(good = False, read = None,
                            errors = np.array([0, 0, 0, 1]))
            record = record[(len_primer+7):]
            record.description += ' alu.seq:' + str(alu)
            return info(good = True, read = record,
                    errors = np.array([0, 0, 0, 0]))
    else:
        return info(good = False, read = None,
                errors = np.array([1, 0, 0, 0]))



def trim_ads(record, ad1, ad2, barlen, m, elem_remove):
    record_seq = record.seq
    len_ad1 = len(ad1)
    len_ad2 = len(ad2)
    seq1 = record_seq[0:len_ad1]
    seq1_shift = record_seq[1:(len_ad1+1)]
    seq2 = record_seq[(len_ad1+barlen):(len_ad1+barlen+len_ad2)]
    seq2_shift = record_seq[(len_ad1+barlen+1):(len_ad1+barlen+len_ad2+1)]
    ham_ad1 = hamming(ad1, seq1, m)
    ham_ad1_shift = hamming(ad1, seq1_shift, m)
    ham_ad2 = hamming(ad2, seq2, m)
    ham_ad2_shift = hamming(ad1, seq2_shift, m)
    if (ham_ad1)and(ham_ad2):
        barcode = record_seq[(len_ad1+1):(len_ad1+barlen+1)]
        for elem in elem_remove:
                if record_seq[(len_ad1+barlen+len_ad2):].find(elem, 0) != -1:
                    return info(good = False, read = None,
                            errors = np.array([0, 0, 0, 1]))
        record = record[(len_ad1+barlen+len_ad2):]
        record.description += ' barcode:' + str(barcode)
        return info(good = True, read = record,
                errors = np.array([0, 0, 0, 0]))
    elif (ham_ad1_shift)and(ham_ad2_shift):
        barcode = record_seq[(len_ad1+1):(len_ad1+barlen+1)]
        for elem in elem_remove:
            if record_seq[(len_ad1+barlen+len_ad2+1):].find(elem, 0) != -1:
                return info(good = False, read = None,
                        errors = np.array([0, 0, 0, 1]))
        record = record[(len_ad1+barlen+len_ad2+1):]
        record.description += ' barcode:' + str(barcode)
        return info(good = True, read = record,
                errors = np.array([0, 0, 0, 0]))
    else:
        if not((ham_ad1)or(ham_ad2)):
            return info(good = False, read = None,
                    errors = np.array([0, 1, 1, 0]))
        elif not(ham_ad1):
            return info(good = False, read = None,
                    errors = np.array([0, 1, 0, 0]))
        else:
            return info(good = False, read = None,
                    errors = np.array([0, 0, 1, 0]))



def concate(x1, x2):
    result = ''
    for i in range(len(x1)):
        result = result + str(x1[i]) + '-' + str(x2[i])
        if (i != (len(x1) - 1)): result = result + ','
    return (result)



def trim_reads(filename1, filename2, inputdir, outputdir, mist, primer, ad1, ad2, barlen, elem_remove):
    readsname = filename1.split('R1')[0]
    readsname = readsname.rsplit('.', 1)[0]
    
    outputfile1, ext = os.path.splitext(filename1)
    outputfile2, ext = os.path.splitext(filename2)
    
    #goodr1 = outputdir + outputfile1 + '_good.fastq'
    #goodr2 = outputdir + outputfile2 + '_good.fastq'
    #badr1 = outputdir + outputfile1 + '_bad.fastq'
    #badr2 = outputdir + outputfile2 + '_bad.fastq'
    
    goodr1 = open(outputdir + outputfile1 + '_good.fastq', 'w')
    goodr2 = open(outputdir + outputfile2 + '_good.fastq', 'w')
    badr1 = open(outputdir + outputfile1 + '_bad.fastq', 'w')
    badr2 = open(outputdir + outputfile2 + '_bad.fastq', 'w')
    
    original_R1_reads = SeqIO.parse(inputdir + filename1, "fastq")
    original_R2_reads = SeqIO.parse(inputdir + filename2, "fastq")
    
    count = np.array([0, 0, 0, 0])
    elem = ('primer', 'ad', 'green', 'flank')
    count_reads = {'readname':readsname, 'all':0, 'good':0, 'bad':0, 'primer':0, 'ad':0, 'green':0}
    for r1, r2 in log_progress(zip(original_R1_reads, original_R2_reads), name = readsname, size = count_fastq_records(inputdir + filename1), every = 250):
    #for r1, r2 in zip(original_R1_reads, original_R2_reads):
        count_reads['all'] += 1
        fr1 = trim_primers(r1, primer, mist, elem_remove)
        if fr1.good:
            fr2 = trim_ads(r2, ad1, ad2, barlen, mist, elem_remove)
            if fr2.good:
                count_reads['good'] += 1
                goodr1.write(fr1.read.format('fastq'))
                goodr2.write(fr1.read.format('fastq'))
            else:
                r2.description += ' reason:' + concate(elem, (np.char.mod('%d', fr2.errors)))
                badr1.write(r1.format('fastq'))
                badr2.write(r2.format('fastq'))
                count = np.sum([count, fr2.errors], axis=0)
        else:
            r2.description += ' reason:' + concate(elem, (np.char.mod('%d', fr1.errors)))
            badr1.write(r1.format('fastq'))
            badr2.write(r2.format('fastq'))
            count = np.sum([count, fr1.errors], axis=0)
    
    goodr1.close()
    goodr2.close()
    badr1.close()
    badr2.close()
    
    count_reads['primer'] = count[0]
    count_reads['ad'] = count[1]
    count_reads['green'] = count[2]
    count_reads['flank'] = count[3]
    
    count_reads['good'] = round((count_reads['good'] / count_reads['all']), 2)
    count_reads['bad'] = round((1 - count_reads['good']), 2)
    count_elem = concate(elem, (np.char.mod('%d', count)))
    print ('For ' + readsname + ': mistake(place-amount) = ' + count_elem + ';  ')
    print ('reads: {0:d}, good: {1:.2f}, bad: {2:.2f}\n'.format(count_reads['all'],
                                                                    count_reads['good'],
                                                                    count_reads['bad']
                                                                   )
          )
    return (count_reads)



def main(inputdir, outputdir, mist, primer, ad1, ad2, barlen, elem_remove):
    inputdir += "/"
    outputdir += "/"

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

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    statistics = open(outputdir + 'statistics.txt', 'w')
    statistics.write('readname\t' + 'reads\t' + 'good.pt\t' + 'bad.pt\t' + 'primer\t' + 'ad\t' + 'green\t' + 'flank\n')
    for filename1, filename2 in conform_files:
        stat_out = trim_reads(filename1, filename2,
                              inputdir, outputdir, mist,
                              primer, ad1, ad2, barlen, elem_remove)
        statistics.write("\t".join([stat_out['readname'], 
                                   str(stat_out['all']), 
                                   str(stat_out['good']),
                                   str(stat_out['bad']),
                                   str(stat_out['primer']),
                                   str(stat_out['ad']),
                                   str(stat_out['green']),
                                   str(stat_out['flank'])]) +
                                   '\n')
    
    statistics.close()
    

    if len(nonconform_files) != 0:
        print ('I can\'t read this files' + str(nonconform_files))


