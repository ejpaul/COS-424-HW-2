'''
Created on Mar 23, 2015

@author: jonathanshor
'''
import nltk, numpy
import sys, time
from optparse import OptionParser

TRAIN_FNAME = "intersected_final_chr1_cutoff_20_train_revised.bed"
SAMPLE_FNAME = "intersected_final_chr1_cutoff_20_sample.bed"
TEST_FNAME = "intersected_final_chr1_cutoff_20_test.bed"

def read_bed_dat(myfile):
    return numpy.loadtxt(myfile, dtype="string", delimiter="\t")

def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--path", dest="path", help='read bed data fom PATH', metavar='PATH')
	(options, args) = parser.parse_args()     
	path = options.path
	print "PATH = " + path
	trains = read_bed_dat(path + TRAIN_FNAME)
	sample = read_bed_dat(path + SAMPLE_FNAME)
    
	print "Sample[0] =%s" % sample[0]
	print "Trains[0] =%s" % trains[0]
	print "len(Trains) =%s" % len(trains)

if __name__ == '__main__':
    main(sys.argv[1:])
