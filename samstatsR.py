import sys, os, re
import subprocess
from os import listdir
from os.path import isfile, join

def main(inputdir, outputdir, filename, statsline):
    inputdir += "/"
    outputdir += "/"

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    logfile = open(outputdir + filename, 'w')

    # Read files in folder
    onlyfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]
    for filename in onlyfiles:
        filename = filename.rstrip()
        inputfile, ext = os.path.splitext(filename)
        if ext == '.sam':
            print('STATISTICS FOR ' + inputfile)
            logfile.write('STATISTICS FOR ' + inputfile + '\n')
            sam_filename = filename
            samline = statsline + ' ' + inputdir + sam_filename + ' | ' + 'grep ^SN | cut -f 2-'
            print (samline)
            p = subprocess.Popen (samline, stdout=subprocess.PIPE, shell = True)
            logline = p.stdout.read().decode()
            print (logline)
            logfile.write(logline)
            logfile.write('\n#################################\n')

    logfile.close()