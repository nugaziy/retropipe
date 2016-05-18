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

def hamming (x1, x2):
    j = 0
    for i in range(len(x1)):
        if x1[i] != x2[i]:
            j += 1
    return (j)

def get_best_alu (x, standart_alu):
    x_counter = dict(Counter(x))
    x_vmax = max(x_counter.values())
    x_max = []
    for key, value in x_counter.items():
        if value == x_vmax:
            x_max.append(key)
    if len(x_max) > 1:
        ham = []
        for i in x_max:
            ham.append(hamming(i, standart_alu))
        ham = np.array(ham)
        nmbr = ham.argmin()
        if len(nmbr) > 1:
            nmbr = random.choice(nmbr)
        x_max = x_max[nmbr]
    x_hmin = hamming(str(x_max[0]), standart_alu)
    result = {'seq' : str(x_max[0]), 'amount' : x_vmax, 'hamming' : x_hmin}
    return(result)

def megaclustering(df, window, megacluster_id, chrom, strand, standart_alu, table1, table2, files_name):
    is_cluster_open = False
    for index, row in df.iterrows():
        if not is_cluster_open:
            is_cluster_open = True
            num_reads_by_files = defaultdict()
            num_barcodes_by_files = defaultdict()
            filecluster_list = defaultdict()
            barcode = defaultdict()
            barcode_q = defaultdict()
            alu = []
            for i in files_name:
                num_reads_by_files[i] = 0
                num_barcodes_by_files[i] = 0
                filecluster_list[i] = '0'
                barcode[i] = ['RRRRRRRRR']
                barcode_q[i] = ['RRRRRRRRR']
            # transform alu_list and extend
            alu_list = np.asarray(row['ALU_LIST'].split(','))
            alu_amount = row['ALU_AMOUNT'].split(',')
            alu_amount = [int(x) for x in alu_amount]
            alu_current = list(np.repeat(alu_list, alu_amount, axis = 0))
            alu.extend(alu_current)
            # transform barcode_list and extend
            barcode_list = re.findall(r'.........', row['BARCODE_LIST'])
            barcode_list_q = re.findall(r'.........', row['BARCODE_Q'])
            barcode[index[0]].extend(barcode_list)
            barcode_q[index[0]].extend(barcode_list_q)
            ###
            megacluster_id += 1
            pos = int(row['POS'])
            pos_list = [int(row['POS'])]
            info = defaultdict(list)
            best_read1 = row['READ1_BEST']
            best_read2 = row['READ2_BEST']
            best_fname = row['FILENAME']
            best_rname = row['READNAME']
            best_tlen = row['TLEN']
            current_mdflag = row['MDFLAG_BEST']
            best_mdflag_int = re.findall(r'\d+', current_mdflag)
            best_mdflag_int = sum([ int(x) for x in best_mdflag_int ])
            best_mdflag = current_mdflag
            best_cigar = row['CIGAR_BEST']
            num_reads_by_files[index[0]] += int(row['NUM_READS'])
            if filecluster_list[index[0]] == '0':
                filecluster_list[index[0]] = str(int(row['CLUSTER_ID']))
            else:
                filecluster_list[index[0]] += ',' + str(int(row['CLUSTER_ID']))
            
        else:
            if abs(int(row['POS']) - pos) <= window:
                current_mdflag = row['MDFLAG_BEST']
                current_mdflag_int = re.findall(r'\d+', current_mdflag)
                current_mdflag_int = sum([ int(x) for x in current_mdflag_int ])
                num_reads_by_files[index[0]] += int(row['NUM_READS'])
                # transform alu_list and extend
                alu_list = np.asarray(row['ALU_LIST'].split(','))
                alu_amount = row['ALU_AMOUNT'].split(',')
                alu_amount = [int(x) for x in alu_amount]
                alu_current = list(np.repeat(alu_list, alu_amount, axis = 0))
                alu.extend(alu_current)
                # transform barcode_list and extend
                barcode_list = re.findall(r'.........', row['BARCODE_LIST'])
                barcode_list_q = re.findall(r'.........', row['BARCODE_Q'])
                barcode[index[0]].extend(barcode_list)
                barcode_q[index[0]].extend(barcode_list_q)
                ###
                if filecluster_list[index[0]] == '0':
                    filecluster_list[index[0]] = str(int(row['CLUSTER_ID']))
                else:
                    filecluster_list[index[0]] += ',' + str(int(row['CLUSTER_ID']))
                if current_mdflag_int > best_mdflag_int:
                    best_read1 = row['READ1_BEST']
                    best_read2 = row['READ2_BEST']
                    best_fname = row['FILENAME']
                    best_rname = row['READNAME']
                    best_tlen = row['TLEN']
                    best_mdflag_int = current_mdflag_int
                    best_mdflag = current_mdflag
                    best_cigar = row['CIGAR_BEST']
            else:
                is_cluster_open = False
                alu_x = dict(Counter(alu).most_common(1))
                alu_best = {'seq' : str(list(alu_x.keys())[0]), 'amount' : list(alu_x.values())[0],
                 'hamming' : hamming(str(list(alu_x.keys())[0]), standart_alu)}
                '''
                pos_density = gaussian_kde(pos_list)
                xs = np.linspace(min(pos_list) - 1, max(pos_list) + 1, len(pos_list) * 100)
                pos_arr = np.column_stack((np.array(xs), np.array(pos_density(xs))))
                best_pos = int(round(pos_arr[pos_arr[:, 1].argmax(), 0], 0))
                '''
                best_pos = list(dict(Counter(pos_list).most_common(1)).keys())[0]
                for i in files_name:
                    num_barcodes_by_files[i] = len(set(barcode[i])) - 1
                    barcode[i] = ''.join(str(x) for x in barcode[i])
                    barcode_q[i] = ''.join(str(x) for x in barcode_q[i])
                cigar_del = sum([ int(x) for x in re.findall(r'(\d+)+D', best_cigar)])
                cigar_ins = sum([ int(x) for x in re.findall(r'(\d+)+I', best_cigar)])
                md_mismatch = len(re.findall('[A-Z]', best_mdflag.split('MD:Z:')[1])) - cigar_del
                md_sum = best_mdflag_int + md_mismatch
                table1.write(str(megacluster_id) + '\t' + best_fname + '\t' + best_rname + '\t' + 
                    'chr' + chrom + '\t' + str(best_pos) + '\t' + 
strand + '\t' + alu_best['seq'] + '\t' + str(alu_best['amount']) + '\t' + str(alu_best['hamming']) + '\t' + best_read1 + '\t' +
best_read2 + '\t' + str(best_tlen) + '\t' + 
best_cigar + '\t' + best_mdflag + '\t' + str(md_sum) + '\t' + str(md_mismatch) + '\t' + str(cigar_ins) + '\t' + str(cigar_del) + '\t' + 
'\t'.join(str(num_reads_by_files[x]) + '\t' + str(num_barcodes_by_files[x]) for x in files_name) + '\n')
                table2.write(str(megacluster_id) + '\t' + ';'.join(str(key) for key, value in filecluster_list.items()) + '\t' + 
                    ';'.join(','.join([value]) for key, value in filecluster_list.items()) + '\t' + best_fname + '\t' +
                    best_rname + '\t' + 'chr' + chrom + '\t' + str(best_pos) + '\t' + 
strand + '\t' + alu_best['seq'] + '\t' + str(alu_best['amount']) + '\t' + str(alu_best['hamming']) + '\t' + best_read1 + '\t' +
best_read2 + '\t' + str(best_tlen) + '\t' + 
best_cigar + '\t' + best_mdflag + '\t' + str(md_sum) + '\t' + str(md_mismatch) + '\t' + str(cigar_ins) + '\t' + str(cigar_del) + '\t' + 
'\t'.join(str(num_reads_by_files[x]) + '\t' + str(num_barcodes_by_files[x]) for x in files_name) + '\t' + 
'\t'.join(str(barcode[x]) + '\t' + str(barcode_q[x]) for x in files_name) + '\n')
                is_cluster_open = True
                num_reads_by_files = defaultdict()
                num_barcodes_by_files = defaultdict()
                filecluster_list = defaultdict()
                barcode = defaultdict()
                barcode_q = defaultdict()
                alu = []
                for i in files_name:
                    num_reads_by_files[i] = 0
                    num_barcodes_by_files[i] = 0
                    filecluster_list[i] = '0'
                    barcode[i] = ['RRRRRRRRR']
                    barcode_q[i] = ['RRRRRRRRR']
                # transform alu_list and extend
                alu_list = np.asarray(row['ALU_LIST'].split(','))
                alu_amount = row['ALU_AMOUNT'].split(',')
                alu_amount = [int(x) for x in alu_amount]
                alu_current = list(np.repeat(alu_list, alu_amount, axis = 0))
                alu.extend(alu_current)
                # transform barcode_list and extend
                barcode_list = re.findall(r'.........', row['BARCODE_LIST'])
                barcode_list_q = re.findall(r'.........', row['BARCODE_Q'])
                barcode[index[0]].extend(barcode_list)
                barcode_q[index[0]].extend(barcode_list_q)
                ###
                megacluster_id += 1
                pos = int(row['POS'])
                pos_list = [int(row['POS'])]
                info = defaultdict(list)
                best_read1 = row['READ1_BEST']
                best_read2 = row['READ2_BEST']
                best_fname = row['FILENAME']
                best_rname = row['READNAME']
                best_tlen = row['TLEN']
                current_mdflag = row['MDFLAG_BEST']
                best_mdflag_int = re.findall(r'\d+', current_mdflag)
                best_mdflag_int = sum([ int(x) for x in best_mdflag_int ])
                best_mdflag = current_mdflag
                best_cigar = row['CIGAR_BEST']
                num_reads_by_files[index[0]] += int(row['NUM_READS'])
                if filecluster_list[index[0]] == '0':
                    filecluster_list[index[0]] = str(int(row['CLUSTER_ID']))
                else:
                    filecluster_list[index[0]] += ',' + str(int(row['CLUSTER_ID']))
    if is_cluster_open:
        alu_x = dict(Counter(alu).most_common(1))
        alu_best = {'seq' : str(list(alu_x.keys())[0]), 'amount' : list(alu_x.values())[0],
         'hamming' : hamming(str(list(alu_x.keys())[0]), standart_alu)}
        '''
        pos_density = gaussian_kde(pos_list)
        xs = np.linspace(min(pos_list) - 1, max(pos_list) + 1, len(pos_list) * 100)
        pos_arr = np.column_stack((np.array(xs), np.array(pos_density(xs))))
        best_pos = int(round(pos_arr[pos_arr[:, 1].argmax(), 0], 0))
        '''
        best_pos = list(dict(Counter(pos_list).most_common(1)).keys())[0]
        for i in files_name:
            num_barcodes_by_files[i] = len(set(barcode[i])) - 1
            barcode[i] = ''.join(str(x) for x in barcode[i])
            barcode_q[i] = ''.join(str(x) for x in barcode_q[i])
        cigar_del = sum([ int(x) for x in re.findall(r'(\d+)+D', best_cigar)])
        cigar_ins = sum([ int(x) for x in re.findall(r'(\d+)+I', best_cigar)])
        md_mismatch = len(re.findall('[A-Z]', best_mdflag.split('MD:Z:')[1])) - cigar_del
        md_sum = best_mdflag_int + md_mismatch
        table1.write(str(megacluster_id) + '\t' + best_fname + '\t' + best_rname + '\t' + 
'chr' + chrom + '\t' + str(best_pos) + '\t' + 
strand + '\t' + alu_best['seq'] + '\t' + str(alu_best['amount']) + '\t' + str(alu_best['hamming']) + '\t' + best_read1 + '\t' +
best_read2 + '\t' + str(best_tlen) + '\t' + 
best_cigar + '\t' + best_mdflag + '\t' + str(md_sum) + '\t' + str(md_mismatch) + '\t' + str(cigar_ins) + '\t' + str(cigar_del) + '\t' + 
'\t'.join(str(num_reads_by_files[x]) + '\t' + str(num_barcodes_by_files[x]) for x in files_name) + '\n')
        table2.write(str(megacluster_id) + '\t' + ';'.join(str(key) for key, value in filecluster_list.items()) + '\t' + 
';'.join(','.join([value]) for key, value in filecluster_list.items()) + '\t' + best_fname + '\t' +
best_rname + '\t' + 'chr' + chrom + '\t' + str(best_pos) + '\t' + 
strand + '\t' + alu_best['seq'] + '\t' + str(alu_best['amount']) + '\t' + str(alu_best['hamming']) + '\t' + best_read1 + '\t' +
best_read2 + '\t' + str(best_tlen) + '\t' + 
best_cigar + '\t' + best_mdflag + '\t' + str(md_sum) + '\t' + str(md_mismatch) + '\t' + str(cigar_ins) + '\t' + str(cigar_del) + '\t' + 
'\t'.join(str(num_reads_by_files[x]) + '\t' + str(num_barcodes_by_files[x]) for x in files_name) + '\t' + 
'\t'.join(str(barcode[x]) + '\t' + str(barcode_q[x]) for x in files_name) + '\n')
    return (megacluster_id)


def main(inputdir, outputdir, outputtable, window, standart_alu):
    #test
    inputdir += "/"
    outputdir += "/"
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    # Read files in folder
    onlyfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]
    
    colnames = ['CLUSTER_ID', 'ID_LIST', 'FILENAME', 'READNAME', 'CHR', 'POS', 'STRAND',
     'READ1_BEST', 'READ2_BEST', 'TLEN', 'CIGAR_BEST', 'MDFLAG_BEST', 'NUM_READS', 'NUM_BARCODES',
      'ALU_LIST', 'ALU_AMOUNT', 'BARCODE_LIST', 'BARCODE_Q']
    outputtable, ext = os.path.splitext(outputtable)
    table1 = open(outputdir + outputtable + '_humanread.txt', 'w')
    table2 = open(outputdir + outputtable + '_pcread.txt', 'w')
    start_point = defaultdict()
    for filename in onlyfiles:
        filename = filename.rstrip()
        inputfile, ext = os.path.splitext(filename)
        if ext == '.txt' and re.search('_bigtable_pcread', inputfile):
            data_index = inputfile.split('_bigtable_pcread')[0]
            start_point[data_index] = 0
    table1.write('MEGACLUSTER_ID\tFILENAME\tREADNAME\tCHR\tPOS\tSTRAND\tALU_BEST\tALU_AMOUNT\tALU_HAMMING\tREAD1_BEST\tREAD2_BEST\t' + 
        'TLEN\tCIGAR_BEST\tMDFLAG_BEST\tMDFLAG_SUM\tMISMATCH\tINSERTION\tDELETION\t' + 
        '\t'.join(x + '_NUM_READS'  + '\t' + x + '_NUM_BARCODES' for x in list(start_point.keys())) + '\n')
    table2.write('MEGACLUSTER_ID\tFILE_IN_MEGACLUSTER\tFILEID_IN_MEGACLUSTER\t'
        'FILENAME\tREADNAME\tCHR\tPOS\tSTRAND\tALU_BEST\tALU_AMOUNT\tALU_HAMMING\tREAD1_BEST\tREAD2_BEST\t' + 
        'TLEN\tCIGAR_BEST\tMDFLAG_BEST\tMDFLAG_SUM\tMISMATCH\tINSERTION\tDELETION\t' + 
        '\t'.join(x + '_NUM_READS'  + '\t' + x + '_NUM_BARCODES' for x in list(start_point.keys())) + '\t' +
        '\t'.join(x + '_BARCODE_LIST' + '\t' + x + '_BARCODE_Q' for x in list(start_point.keys())) + '\t' + '\n')
    megacluster_id = 0
    chromline = list(range(1, 23))
    chromline.append('X')
    chromline.append('Y')
    for i_chrom in log_progress(chromline, name = 'megatable', every = 1, who = 'classes: chrom'):
        bigdata_chrom = pd.DataFrame(columns = colnames)
        for filename in onlyfiles:
            filename = filename.rstrip()
            inputfile, ext = os.path.splitext(filename)
            if ext == '.txt' and re.search('_bigtable_pcread', inputfile):
                data_index = inputfile.split('_bigtable_pcread')[0]
                data = pd.read_table(inputdir + filename, skiprows = 1, header = None)
                data.columns = colnames
                data_chrom = data[data['CHR'] == 'chr' + str(i_chrom)]
                start_point[data_index] += len(data_chrom)
                bigdata_chrom = pd.concat([bigdata_chrom, pd.concat([data_chrom], keys=[data_index])])
        bigdata_chrom = bigdata_chrom.groupby(['STRAND'])
        for name, group in bigdata_chrom:
            group = group.sort_values(['POS'])
            megacluster_id = megaclustering(group, window, megacluster_id, str(i_chrom), name, standart_alu,
             table1, table2, list(start_point.keys()))
    table1.close()
    table2.close()
