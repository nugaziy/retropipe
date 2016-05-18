from Bio import SeqIO
import sys, os, re
import numpy as np
from os import listdir
from os.path import isfile, join
from utils import *
import subprocess
from collections import namedtuple
import multiprocessing as mp

info = namedtuple('info', 'good read alu_barcode errors')

def hamming (x1, x2, m):
    j = 0
    for i in range(len(x1)):
        if x1[i] != x2[i]:
            j += 1
            if j > m:
                return (False)
    return (True)

def trim_primers (record, primer, shift, m, elem_remove, search_win):
    #test
    record_seq = record.seq
    len_primer = len(primer)
    pr_elem_remove = primer
    for i in range(shift):
        ham_prim = hamming(primer, record_seq[i : len_primer + i], m)
        if ham_prim:
            alu = record_seq[len_primer + i : len_primer + 6 + i]
            for elem in elem_remove:
                if record_seq[len_primer + 6 + i :].find(elem, 0) != -1:
                    return (info(good = False, read = None, alu_barcode = None,
                            errors = np.array([0, 0, 0, 1, 0, 0])))
            if record_seq[len_primer + 6 + i : len_primer + 6 + i + search_win].find(pr_elem_remove, 0) != -1:
                    return (info(good = False, read = None, alu_barcode = None,
                            errors = np.array([0, 0, 0, 0, 1, 0])))
            record = record[len_primer + 6 + i :]
            alu_bar = '__abq:' + str(alu)
            record.description = ''
            record.name = ''
            return (info(good = True, read = record, alu_barcode = alu_bar,
                    errors = np.array([0, 0, 0, 0, 0, 0])))
    return (info(good = False, read = None, alu_barcode = None,
                errors = np.array([1, 0, 0, 0, 0, 0])))

def trim_ads (record, ad1, ad2, barlen, shift, m, elem_remove, search_win):
    record_seq = record.seq
    len_ad1 = len(ad1)
    len_ad2 = len(ad2)
    ad_elem_remove = [ad1, ad2]
    for i in range(shift):
        seq1 = record_seq[i : len_ad1 + i]
        seq2 = record_seq[len_ad1 + barlen + i : len_ad1 + barlen + len_ad2 + i]
        ham_ad1 = hamming(ad1, seq1, m)
        ham_ad2 = hamming(ad2, seq2, m)
        if (ham_ad1)and(ham_ad2):
            barcode = record_seq[len_ad1 + i : len_ad1 + barlen + i]
            barcode_q = [chr(x + 33) for x in record.letter_annotations['phred_quality']]
            barcode_q = ''.join(barcode_q)
            barcode_q = barcode_q[len_ad1 + i : len_ad1 + barlen + i]
            for elem in elem_remove:
                if record_seq[len_ad1 + barlen + len_ad2 + i :].find(elem, 0) != -1:
                    return (info(good = False, read = None, alu_barcode = None,
                            errors = np.array([0, 0, 0, 1, 0, 0])))
            for ad_elem in ad_elem_remove:
                if record_seq[len_ad1 + barlen + len_ad2 + i : len_ad1 + barlen + len_ad2 + i + search_win].find(ad_elem, 0) != -1:
                    return (info(good = False, read = None, alu_barcode = None,
                        errors = np.array([0, 0, 0, 0, 1, 0])))
            record = record[len_ad1 + barlen + len_ad2 + i :]
            alu_bar = '__abq:' + str(barcode) + '__abq:' + str(barcode_q)
            record.description = ''
            record.name = ''
            if record_seq[len_ad1 + barlen + len_ad2 + i : len_ad1 + barlen + len_ad2 + i + 2] == 'CT':
                return (info(good = True, read = record, alu_barcode = alu_bar,
                        errors = np.array([0, 0, 0, 0, 0, 0])))
            else:
                mb_return = info(good = False, read = None, alu_barcode = None,
                    errors = np.array([0, 0, 0, 0, 0, 1]))
        elif not((ham_ad1)or(ham_ad2)):
            mb_return = info(good = False, read = None, alu_barcode = None,
                    errors = np.array([0, 1, 1, 0, 0, 0]))
        elif not(ham_ad1):
            mb_return = info(good = False, read = None, alu_barcode = None,
                    errors = np.array([0, 1, 0, 0, 0, 0]))
        else:
            mb_return = info(good = False, read = None, alu_barcode = None,
                    errors = np.array([0, 0, 1, 0, 0, 0]))
    return (mb_return)



def concate(x1, x2):
    result = ''
    for i in range(len(x1)):
        result = result + str(x1[i]) + '-' + str(x2[i])
        if (i != (len(x1) - 1)): result = result + ','
    return (result)



def trim_reads(filename1, filename2, inputdir, outputdir, shift,
 mist1, mist2, primer, ad1, ad2, barlen, elem_remove, search_win):
    readsname = filename1.split('R1')[0]
    readsname = readsname.rsplit('.', 1)[0]
    
    outputfile1, ext = os.path.splitext(filename1)
    outputfile2, ext = os.path.splitext(filename2)
    
    goodr1 = open(outputdir + outputfile1 + '_good.fq', 'w')
    goodr2 = open(outputdir + outputfile2 + '_good.fq', 'w')
    badr1 = open(outputdir + outputfile1 + '_bad.fq', 'w')
    badr2 = open(outputdir + outputfile2 + '_bad.fq', 'w')
    
    original_R1_reads = SeqIO.parse(inputdir + filename1, "fastq")
    original_R2_reads = SeqIO.parse(inputdir + filename2, "fastq")
    
    count = np.array([0, 0, 0, 0, 0, 0])
    elem = ('primer', 'ad', 'green', 'flank_simple', 'flank_strange', 'non_ct')
    count_reads = {'readname':readsname, 'all':0, 'good':0, 'bad':0, 'primer':0, 'ad':0, 'green':0, 'flank_simple':0,
    'flank_strange':0, 'non_ct':0}
    for r1, r2 in log_progress(zip(original_R1_reads, original_R2_reads), name = readsname, size = count_fastq_records(inputdir + filename1), every = 250):
        count_reads['all'] += 1
        fr1 = trim_primers(r1, primer, shift, mist1, elem_remove, search_win)
        if fr1.good:
            fr2 = trim_ads(r2, ad1, ad2, barlen, shift, mist2, elem_remove, search_win)
            if fr2.good:
                count_reads['good'] += 1
                fr1.read.id += fr1.alu_barcode + fr2.alu_barcode
                fr2.read.id += fr1.alu_barcode + fr2.alu_barcode
                goodr1.write(fr1.read.format('fastq'))
                goodr2.write(fr2.read.format('fastq'))
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
    count_reads['flank_simple'] = count[3]
    count_reads['flank_strange'] = count[4]
    count_reads['non_ct'] = count[5]
    
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



def main(inputdir, outputdir, shift,
 mist1, mist2, primer, ad1, ad2, barlen, elem_remove, search_win, statfilename):
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
    
    statistics = open(outputdir + statfilename, 'w')
    statistics.write('readname\t' + 'reads\t' + 'good.pt\t' + 'bad.pt\t' + 'primer\t' + 'ad\t' + 'green\t' + 'flank_simple\t' + 
        'flank_strange' + '\t' + 'non_ct\n')

    '''
    # Setup a list of processes that we want to run
    processes = [mp.Process(target=trim_reads, args=(filename1, filename2,
                              inputdir, outputdir, mist1, mist2,
                              primer, ad1, ad2, barlen, elem_remove)) for filename1, filename2 in conform_files]

    # Run processes
    for p in processes:
        p.start()
 
    # Exit the completed processes
    for p in processes:
        p.join()

    for p in processes:
        stat_out = output.get()
        statistics.write("\t".join([stat_out['readname'], 
                                   str(stat_out['all']), 
                                   str(stat_out['good']),
                                   str(stat_out['bad']),
                                   str(stat_out['primer']),
                                   str(stat_out['ad']),
                                   str(stat_out['green']),
                                   str(stat_out['flank_simple'])]),
                                   str(stat_out['flank_strange']) +
                                   '\n')

    '''

    for filename1, filename2 in conform_files:

        ext1, inputfile1 = os.path.splitext(filename1)
        ext2, inputfile2 = os.path.splitext(filename2)
        '''
        is_ext_gz = False
        if ext1 == 'gz':
            print ('unpack ' + inputfile1 + '\t')
            if not os.path.exists(outputdir + 'forunpackingfastq'):
                os.makedirs(outputdir + 'forunpackingfastq')
            unpack_line = 'tar -zxvf ' + inputdir + filename1 + '> ' + outputdir + forunpackingfastq + '/r1_unpack.fastq'
            new_input1 = outputdir + forunpackingfastq + '/'
            is_ext_gz = True
            p = subprocess.Popen (unpack_line, stdout=subprocess.PIPE, shell = True)
            p.stdout
        else: new_input = ''
        if ext2 == 'gz':
            if is_ext_gz:
                print ('unpack ' + inputfile2 + '\t')
                if not os.path.exists(outputdir + 'forunpackingfastq'):
                    os.makedirs(outputdir + 'forunpackingfastq')
                unpack_line = 'tar -zxvf ' + inputdir + filename2 + '> ' + outputdir + forunpackingfastq + '/r2_unpack.fastq'
                is_ext_gz = True
                p = subprocess.Popen (unpack_line, stdout=subprocess.PIPE, shell = True)
                p.stdout
            else:
                cp_line = 'cp ' + inputdir + filename1 + '> ' + outputdir + forunpackingfastq + '/r1_unpack.fastq'
        else:
            if is_ext_gz:
                cp_line = 'cp ' + inputdir + filename2 + '> ' + outputdir + forunpackingfastq + '/r2_unpack.fastq'
                p = subprocess.Popen (cp_line, stdout=subprocess.PIPE, shell = True)
                p.stdout

        '''
        stat_out = trim_reads(filename1, filename2, inputdir, outputdir, shift,
 mist1, mist2, primer, ad1, ad2, barlen, elem_remove, search_win)
        statistics.write("\t".join([stat_out['readname'], 
                                    str(stat_out['all']), 
                                    str(stat_out['good']),
                                    str(stat_out['bad']),
                                    str(stat_out['primer']),
                                    str(stat_out['ad']),
                                    str(stat_out['green']),
                                    str(stat_out['flank_simple']),
                                    str(stat_out['flank_strange']),
                                    str(stat_out['non_ct'])]) + 
                                    '\n')

    statistics.close()
    

    if len(nonconform_files) != 0:
        print ('I can\'t read this files' + str(nonconform_files))


