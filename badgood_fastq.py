# Import module (main code)
import trimmR
#import imp
#imp.reload(trimmR)

# Variable parameters
#################################
# The number of permissible error:
mist = 1
# Primer, ad1 = Adapter 1, ad2 = Adapter 3 (aka Green)
primer = 'GAGCCACCGCGC'
ad1 = 'GCGTGCTGCGG'
ad2 = 'AGGGCGGT'
# Length of barcode
barlen = 9
# List of wrong elements in flank
elem_remove = ['ACGT']

# Input FASTQ files folder path.
inputdir = 'input/'
# Output folder for processed FASTQ files.
outputdir = 'output/'
#################################

# Main function

trimmR.main(inputdir, outputdir, mist, primer, ad1, ad2, barlen, elem_remove)
