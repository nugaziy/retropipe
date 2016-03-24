import sys, os, re
import subprocess
from os import listdir
from os.path import isfile, join

def main(inputdir, refway, outputdir, memline):
    inputdir += "/"
    outputdir += "/"

    # Read files in folder
    onlyfiles = [f for f in listdir(inputdir) if isfile(join(inputdir, f))]
    
    r1_files = {}
    r2_files = {}
    
    for filename in onlyfiles:
        filename = filename.rstrip()
        if re.search('good', filename):
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

    logfile = open(outputdir + 'bwamem_output.log', 'w')

    for filename1, filename2 in conform_files:
        readsname = filename1.split('R1')[0]
        readsname = readsname.rsplit('.', 1)[0]
        memline = memline + ' ' + refway + ' ' + inputdir + filename1 + ' ' + inputdir + filename2 + ' > ' + outputdir + readsname + '.sam'
        print (memline)
        p = subprocess.Popen (memline, stderr=subprocess.PIPE, shell = True)
        logline = p.stderr.read().decode()
        print (logline)
        logfile.write(logline)
        logfile.write('\n#################################\n')

    logfile.close()

    if len(nonconform_files) != 0:
        print ('I can\'t read this files' + str(nonconform_files))
