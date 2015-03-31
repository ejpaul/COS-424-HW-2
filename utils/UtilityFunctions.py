'''
Created on Mar 23, 2015

@author: jonathanshor
@author: epaul626
'''
import numpy as np
import sys
from optparse import OptionParser

TRAIN_FNAME = "intersected_final_chr1_cutoff_20_train_revised.bed"
SAMPLE_FNAME = "intersected_final_chr1_cutoff_20_sample.bed"
TEST_FNAME = "intersected_final_chr1_cutoff_20_test.bed"

# Expects two 'Beta' both of size nX1
# Accepts an intercept and calc's accordingly
# Returns (Coeff of Determination (r^2), RMSE)
def calc_r2_RMSE(preds, gTruth, intercept = 0):
	samps = len(preds)
	if samps != len(gTruth):
		print "Incompatible 'Beta's' passed to calc_r2"
		return (np.nan, np.nan)
	else:
		rss = 0.
		tss = 0.
		emean = np.mean(gTruth)
		for i in range(0, samps):
			rss += (gTruth[i] - preds[i] - intercept) ** 2
			tss += (gTruth[i] - emean) ** 2
		mse = rss / samps
		rmse = mse ** 0.5
		r2 = 1 - (rss / tss)
		return (r2, rmse)

# Accepts path to location of SAMPLE_FNAME
def read_bed_dat_sample(mypath):
	return np.loadtxt(mypath + SAMPLE_FNAME, dtype=[('Chrom', np.str_, 4), ('Start', np.int32), \
									('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32), ('450k', np.int8)])

# Accepts path to location of TEST_FNAME
def read_bed_dat_test(mypath):
	return np.loadtxt(mypath + TEST_FNAME, dtype=[('Chrom', np.str_, 4), ('Start', np.int32), \
									('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32), ('450k', np.int8)])

# Accepts path to location of TRAIN_FNAME
def read_bed_dat_train(mypath):
	return np.loadtxt(mypath + TRAIN_FNAME, dtype=[('Chrom', np.str_, 4), ('Start', np.int32), \
									('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32, (33)), ('450k', np.int8)])

def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--path", dest="path", help='read bed data fom PATH', metavar='PATH')
	(options, args) = parser.parse_args()     
	path = options.path
	print "PATH = " + path
	train = read_bed_dat_train(path)
	sample = read_bed_dat_sample(path)

	print "sample[0] = %s" % sample[0]
	print "train[0] = %s" % train[0]
	print "len(train) = %s" % len(train)
	print "train['Beta'][0] = %s" % train['Beta'][0]

if __name__ == '__main__':
	main(sys.argv[1:])
