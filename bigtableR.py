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


def clustering(df, window, cluster_id, chrom, strand, table1, table2):
    is_cluster_open = False
    for index, row in df.iterrows():
        if not is_cluster_open:
            is_cluster_open = True
            cluster_id += 1
            pos = row['POS']
            pos_list = [row['POS']]
            info = defaultdict(list)
            best_read1 = row['READ1']
            best_read2 = row['READ2']
            best_tlen = row['TLEN']
            current_mdflag = row['MDFLAG_R1']
            best_mdflag_int = re.findall(r'\d+', current_mdflag)
            best_mdflag_int = sum([ int(x) for x in best_mdflag_int ])
            best_mdflag = current_mdflag
            best_cigar = row['CIGAR_R1']
            info['id'].append(row['ID'])
            info['barcode'].append(row['BARCODE'])
            info['alu'].append(row['ALU'])
        else:
            if abs(row['POS'] - pos) <= window:
                pos_list.append(row['POS'])
                info['id'].append(row['ID'])
                info['barcode'].append(row['BARCODE'])
                info['alu'].append(row['ALU'])
                current_mdflag = row['MDFLAG_R1']
                current_mdflag_int = re.findall(r'\d+', current_mdflag)
                current_mdflag_int = sum([ int(x) for x in current_mdflag_int ])
                if current_mdflag_int > best_mdflag_int:
                    best_read1 = row['READ1']
                    best_read2 = row['READ2']
                    best_tlen = row['TLEN']
                    best_mdflag_int = current_mdflag_int
                    best_mdflag = current_mdflag
                    best_cigar = row['CIGAR_R1']
            else:
                is_cluster_open = False
                if len(set(pos_list))> 1 and len(set(info['alu'])) > 1:
                    best_pos = list(dict(Counter(pos_list).most_common(1)).keys())[0]
                    '''
                    pos_density = gaussian_kde(pos_list)
                    xs = np.linspace(min(pos_list) - 1, max(pos_list) + 1, len(pos_list) * 100)
                    pos_arr = np.column_stack((np.array(xs), np.array(pos_density(xs))))
                    best_pos = int(round(pos_arr[pos_arr[:, 1].argmax(), 0], 0))
                    '''
                else:
                    best_pos = pos_list[0]
                info['barcode'] = list(set(info['barcode']))
                alu_dict = dict(Counter(info['alu']))
                table1.write(str(cluster_id) + '\t' + chrom + '\t' + str(best_pos) + '\t' + strand + '\t' + 
','.join(str(key) for key, value in alu_dict.items()) + '\t' + ','.join(str(value) for key, value in alu_dict.items()) + '\t' + 
best_read1 + '\t' + best_read2 + '\t' + str(best_tlen) + '\t' + best_cigar + '\t' + best_mdflag + '\t' + str(len(info['id'])) + '\t' +
str(len(info['barcode'])) + '\n')
                table2.write(str(cluster_id) + '\t' + ','.join(str(x) for x in info['id']) + '\t' + 
','.join(str(x) for x in info['barcode']) + '\n')
                is_cluster_open = True
                cluster_id += 1
                pos = row['POS']
                pos_list = [row['POS']]
                info = defaultdict(list)
                best_read1 = row['READ1']
                best_read2 = row['READ2']
                best_tlen = row['TLEN']
                current_mdflag = row['MDFLAG_R1']
                best_mdflag_int = re.findall(r'\d+', current_mdflag)
                best_mdflag_int = sum([ int(x) for x in best_mdflag_int ])
                best_mdflag = current_mdflag
                best_cigar = row['CIGAR_R1']
                info['id'].append(row['ID'])
                info['barcode'].append(row['BARCODE'])
                info['alu'].append(row['ALU'])
    if is_cluster_open:
        if len(set(pos_list))> 1 and len(set(info['alu'])) > 1:
            best_pos = list(dict(Counter(pos_list).most_common(1)).keys())[0]
            '''
            pos_density = gaussian_kde(pos_list)
            xs = np.linspace(min(pos_list) - 1, max(pos_list) + 1, len(pos_list) * 100)
            pos_arr = np.column_stack((np.array(xs), np.array(pos_density(xs))))
            best_pos = int(round(pos_arr[pos_arr[:, 1].argmax(), 0], 0))
            '''
        else:
            best_pos = pos_list[0]
        info['barcode'] = list(set(info['barcode']))
        alu_dict = dict(Counter(info['alu']))
        table1.write(str(cluster_id) + '\t' + chrom + '\t' + str(best_pos) + '\t' + strand + '\t' + 
','.join(str(key) for key, value in alu_dict.items()) + '\t' + ','.join(str(value) for key, value in alu_dict.items()) + '\t' + 
best_read1 + '\t' + best_read2 + '\t' + str(best_tlen) + '\t' + best_cigar + '\t' + best_mdflag + '\t' + str(len(info['id'])) + '\t' +
str(len(info['barcode'])) + '\n')
        table2.write(str(cluster_id) + '\t' + ','.join(str(x) for x in info['id']) + '\t' + 
','.join(str(x) for x in info['barcode']) + '\n')
    return (cluster_id)


def main(inputdir, outputdir, window):
    inputdir += "/"
    outputdir += "/"
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    # Read files in folder
    onlyfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]
    
    for filename in onlyfiles:
        filename = filename.rstrip()
        inputfile, ext = os.path.splitext(filename)
        if ext == '.txt' and re.search('_table', inputfile):
            table1_name = filename.split('_table')[0] + '_bigtable_all.txt'
            table2_name = filename.split('_table')[0] + '_bigtable_baralu.txt'
            table1 = open(outputdir + table1_name, 'w')
            table1.write('CLUSTER_ID\tCHR\tPOS\tSTRAND\tALU_LIST\tALU_AMOUNT\tREAD1_BEST\tREAD2_BEST\tTLEN\t' + 
                'CIGAR_BEST\tMDFLAG_BEST\tNUM_BARCODES\tNUM_READS\n')
            table2 = open(outputdir + table2_name, 'w')
            table2.write('CLUSTER_ID\tID_LIST\tBARCODE_LIST\n')
            data = pd.read_table(inputdir + filename)
            data_group = data.groupby(['CHR', 'STRAND'])
            cluster_id = 0
            for name, group in log_progress(data_group, name = inputfile, every = 1, who = 'classes: chrom & strand'):
                group = group.sort_values(['POS'])
                cluster_id = clustering(group, window, cluster_id, name[0], name[1], table1, table2)
            print ('Done ' + inputfile + '\n')
            table1.close()
            table2.close()

