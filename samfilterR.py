from simplesam import Reader, Writer
import sys, os, re
from utils import *
from os import listdir
from os.path import isfile, join

def sam2table(inputdir, in_sam_filename, outputdir, tablefilename, out_sam_errorfilename):
    in_sam_file = open(inputdir + in_sam_filename, 'r')
    tablefile = open(outputdir + tablefilename, 'w')
    out_sam_errorfile = open(outputdir + out_sam_errorfilename, 'w')
    in_sam = Reader(in_sam_file)
    out_sam_error = Writer(out_sam_errorfile)
    tablefile.write('ID\tCHR\tPOS\tSTRAND\tREAD1\tREAD2\tBARCODE\tALU\tCIGAR_R1\tCIGAR_R2\tMDFLAG_R1\tMDFLAG_R2\n')
    id_count = 0
    for read in log_progress(in_sam, name = in_sam_filename, size = count_fastq_records(inputdir + in_sam_filename), every = 250):
        if read.flag in [99, 83]:
            read1 = read
            read1_bool = True
        elif read.flag in [147, 163]:
            read2 = read
            if read1_bool:
                read1_bool = False
                barcode = read2.qname.split('__ab:')[2]
                alu_seq = read1.qname.split('__ab:')[1]
                if read1.reverse:
                    # it's strand of insertion (not read)
                    strand = '+'
                    pos = read2.pos + abs(read1.tlen) - 1
                else:
                    strand = '-'
                    pos = read1.pos
                id_count += 1
                mdflag_r1 = read1._tags[1]
                mdflag_r2 = read2._tags[1]
                tablefile.write(str(int(id_count)) + '\t' + read1.rname + '\t' + str(pos) + '\t' + strand + '\t' + 
read1.seq + '\t' +read2.seq + '\t' + barcode + '\t' + alu_seq + '\t' + 
read1.cigar + '\t' + read2.cigar + '\t' + mdflag_r1 + '\t' + mdflag_r2 + '\n')
        else:
            out_sam_error.write(read)
    in_sam_file.close()
    tablefile.close()
    out_sam_errorfile.close()

def main(inputdir, outputdir):
    inputdir += "/"
    outputdir += "/"
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    # Read files in folder
    onlyfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]
    for filename in onlyfiles:
        filename = filename.rstrip()
        inputfile, ext = os.path.splitext(filename)
        if ext == '.sam':
            in_sam_filename = filename
            tablefilename = inputfile + '_table.txt'
            out_sam_errorfilename = inputfile + '_error.sam'
            sam2table(inputdir, in_sam_filename, outputdir, tablefilename, out_sam_errorfilename)
