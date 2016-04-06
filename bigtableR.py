import pandas as pd
import re
from collections import defaultdict
import sys, os, re
from utils import *
from os import listdir
from os.path import isfile, join

def clustering(df, perm_gap, cluster_id, chrom, strand, table1, table2):
    is_cluster_open = False
    for index, row in df.iterrows():
        if not is_cluster_open:
            is_cluster_open = True
            cluster_id += 1
            start = row['START']
            end = row['END']
            info = defaultdict(list)
            best_read = row['READ1']
            current_mdflag = row['MDFLAG_R1']
            best_mdflag_int = re.findall(r'\d+', current_mdflag)
            best_mdflag_int = sum([ int(x) for x in best_mdflag_int ])
            best_mdflag = current_mdflag
            best_cigar = row['CIGAR_R1']
            info['id'].append(row['ID'])
            info['barcode'].append(row['BARCODE'])
            info['alu'].append(row['ALU'])
        else:
            if abs(row['START'] - start) <= perm_gap:
                info['id'].append(row['ID'])
                info['barcode'].append(row['BARCODE'])
                info['alu'].append(row['ALU'])
                current_mdflag_int = re.findall(r'\d+', current_mdflag)
                current_mdflag_int = sum([ int(x) for x in current_mdflag_int ])
                if current_mdflag_int > best_mdflag_int:
                    best_read = row['READ1']
                    best_mdflag_int = current_mdflag_int
                    best_mdflag = current_mdflag
                    best_cigar = row['CIGAR_R1']
            else:
                is_cluster_open = False
                info['barcode'] = list(set(info['barcode']))
                info['alu'] = list(set(info['alu']))
                table1.write(str(cluster_id) + '\t' + chrom + '\t' + strand + '\t' + 
str(start) + '\t' + str(end) + '\t' + best_read + '\t' + best_cigar + '\t' + best_mdflag + '\t' +
str(len(info['barcode'])) + '\t' + str(len(info['id'])) + '\n')
                table2.write(str(cluster_id) + '\t' + ','.join(str(x) for x in info['id']) + '\t' + 
','.join(str(x) for x in info['barcode']) + '\t' + ','.join(str(x) for x in info['alu']) + '\n')
                is_cluster_open = True
                cluster_id += 1
                start = row['START']
                end = row['END']
                info = defaultdict(list)
                best_read = row['READ1']
                current_mdflag = row['MDFLAG_R1']
                best_mdflag_int = re.findall(r'\d+', current_mdflag)
                best_mdflag_int = sum([ int(x) for x in best_mdflag_int ])
                best_mdflag = current_mdflag
                best_cigar = row['CIGAR_R1']
                info['id'].append(row['ID'])
                info['barcode'].append(row['BARCODE'])
                info['alu'].append(row['ALU'])
    if is_cluster_open:
        info['barcode'] = list(set(info['barcode']))
        info['alu'] = list(set(info['alu']))
        table1.write(str(cluster_id) + '\t' + chrom + '\t' + strand + '\t' + 
str(start) + '\t' + str(end) + '\t' + best_read + '\t' + best_cigar + '\t' + best_mdflag + '\t' +
str(len(info['barcode'])) + '\t' + str(len(info['id'])) + '\n')
        table2.write(str(cluster_id) + '\t' + ','.join(str(x) for x in info['id']) + '\t' + 
','.join(str(x) for x in info['barcode']) + '\t' + ','.join(str(x) for x in info['alu']) + '\n')
    return (cluster_id)


def main(inputdir, outputdir, perm_gap):
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
            table1.write('CLUSTER_ID\tCHR\tSTRAND\tSTART\tEND\tREAD1_BEST\tCIGAR_BEST\tMDFLAG_BEST\tNUM_BARCODES\tNUM_READS\n')
            table2 = open(outputdir + table2_name, 'w')
            table2.write('CLUSTER_ID\tID_LIST\tBARCODE_LIST\tALU_LIST\n')
            data = pd.read_table(inputdir + filename)
            data_group = data.groupby(['CHR', 'STRAND'])
            cluster_id = 0
            for name, group in data_group:
                group = group.sort_values(['START'])
                cluster_id = clustering(group, perm_gap, cluster_id, name[0], name[1], table1, table2)
            print ('Done ' + inputfile + '\n')
            table1.close()
            table2.close()
