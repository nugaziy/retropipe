import pandas as pd
import re
from collections import defaultdict
import sys, os, re
from utils import *
from os import listdir
from os.path import isfile, join

def megaclustering(df, window, megacluster_id, chrom, strand, table, files_name):
    is_cluster_open = False
    for index, row in df.iterrows():
        if not is_cluster_open:
            is_cluster_open = True
            num_reads_by_files = defaultdict()
            num_barcodes_by_files = defaultdict()
            for i in files_name:
                num_reads_by_files[i] = 0
                num_barcodes_by_files[i] = 0
            megacluster_id += 1
            pos = int(row['POS'])
            best_read = row['READ1_BEST']
            current_mdflag = row['MDFLAG_BEST']
            best_mdflag_int = re.findall(r'\d+', current_mdflag)
            best_mdflag_int = sum([ int(x) for x in best_mdflag_int ])
            best_mdflag = current_mdflag
            best_cigar = row['CIGAR_BEST']
            num_reads_by_files[index[0]] += int(row['NUM_READS'])
            num_barcodes_by_files[index[0]] += int(row['NUM_BARCODES'])
        else:
            if abs(int(row['POS']) - pos) <= window:
                current_mdflag_int = re.findall(r'\d+', current_mdflag)
                current_mdflag_int = sum([ int(x) for x in current_mdflag_int ])
                num_reads_by_files[index[0]] += int(row['NUM_READS'])
                num_barcodes_by_files[index[0]] += int(row['NUM_BARCODES'])
                if current_mdflag_int > best_mdflag_int:
                    best_read = row['READ1_BEST']
                    best_mdflag_int = current_mdflag_int
                    best_mdflag = current_mdflag
                    best_cigar = row['CIGAR_BEST']
            else:
                is_cluster_open = False
                table.write(str(megacluster_id) + '\t' + 'chr' + chrom + '\t' + str(pos) + '\t' + 
strand + '\t' + best_read + '\t' + best_cigar + '\t' + best_mdflag + '\t' + 
'\t'.join(str(num_reads_by_files[x]) for x in files_name) + '\t' + '\t'.join(str(num_barcodes_by_files[x]) for x in files_name) + '\n')
                is_cluster_open = True
                num_reads_by_files = defaultdict()
                num_barcodes_by_files = defaultdict()
                for i in files_name:
                    num_reads_by_files[i] = 0
                    num_barcodes_by_files[i] = 0
                megacluster_id += 1
                pos = int(row['POS'])
                best_read = row['READ1_BEST']
                current_mdflag = row['MDFLAG_BEST']
                best_mdflag_int = re.findall(r'\d+', current_mdflag)
                best_mdflag_int = sum([ int(x) for x in best_mdflag_int ])
                best_mdflag = current_mdflag
                best_cigar = row['CIGAR_BEST']
                num_reads_by_files[index[0]] += int(row['NUM_READS'])
                num_barcodes_by_files[index[0]] += int(row['NUM_BARCODES'])
    if is_cluster_open:
        table.write(str(megacluster_id) + '\t' + 'chr' + chrom + '\t' + str(pos) + '\t' + 
strand + '\t' + best_read + '\t' + best_cigar + '\t' + best_mdflag + '\t' + 
'\t'.join(str(num_reads_by_files[x]) for x in files_name) + '\t' + '\t'.join(str(num_barcodes_by_files[x]) for x in files_name) + '\n')
    return (megacluster_id)


def main(inputdir, outputdir, window):
	#test
    inputdir += "/"
    outputdir += "/"
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    # Read files in folder
    onlyfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]
    
    colnames = ['CLUSTER_ID', 'CHR', 'POS', 
'STRAND', 'READ1_BEST', 'CIGAR_BEST', 'MDFLAG_BEST', 'NUM_BARCODES', 'NUM_READS']
    table = open(outputdir + 'megatable.txt', 'w')
    start_point = defaultdict()
    for filename in onlyfiles:
        filename = filename.rstrip()
        inputfile, ext = os.path.splitext(filename)
        if ext == '.txt' and re.search('_bigtable_all', inputfile):
            data_index = inputfile.split('_bigtable_all')[0]
            start_point[data_index] = 0
    table.write('MEGACLUSTER_ID\tCHR\tPOS\tSTRAND\tREAD1_BEST\tCIGAR_BEST\tMDFLAG_BEST\t' + 
'\t'.join(x + 'NUM_READS' for x in list(start_point.keys())) + '\t' + 
'\t'.join(x + 'NUM_BARCODES' for x in list(start_point.keys())) + '\n')
    megacluster_id = 0
    for i_chrom in log_progress(range(22), name = 'megatable', every = 1, who = 'classes: chrom'):
        bigdata_chrom = pd.DataFrame(columns = colnames)
        for filename in onlyfiles:
            filename = filename.rstrip()
            inputfile, ext = os.path.splitext(filename)
            if ext == '.txt' and re.search('_bigtable_all', inputfile):
                data_index = inputfile.split('_bigtable_all')[0]
                data = pd.read_table(inputdir + filename, skiprows = start_point[data_index], header = None)
                data.columns = colnames
                data_chrom = data[data['CHR'] == 'chr' + str(i_chrom + 1)]
                start_point[data_index] += len(data_chrom)
                bigdata_chrom = pd.concat([bigdata_chrom, pd.concat([data_chrom], keys=[data_index])])
        bigdata_group = bigdata_chrom.groupby(['STRAND'])
        for name, group in bigdata_group:
            group = group.sort_values(['POS'])
            megacluster_id = megaclustering(group, window, megacluster_id, str(i_chrom + 1), name, table, list(start_point.keys()))
    table.close()
