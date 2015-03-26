'''
Created on Mar 24, 2015

@author: epaul626

Utility function for filling in nans for training
'''
import numpy as np
import sys
from scipy.spatial import cKDTree
from UtilityFunctions import * 

TRAIN_FNAME = "intersected_final_chr1_cutoff_20_train_revised.bed"
SAMPLE_FNAME = "intersected_final_chr1_cutoff_20_sample.bed"
TEST_FNAME = "intersected_final_chr1_cutoff_20_test.bed"

# Fills in nans by averaging over 32 samples in training bed
# Pass in beta array for training bed
# Returns complete array
def fill_mean(train_beta):
	train_beta_filled = train_beta
	if (train_beta.shape[1] != 33):
		print "Array must have 33 samples"
		sys.exit(2)
	for i in range(0, len(train_beta)):
		# Calculate mean of non-nan samples
		mean = np.mean(train_beta[i][~np.isnan(train_beta[i])])
		# Replace nan samples with mean
		for j in range(0, 33):
			if (np.isnan(train_beta[i,j])):
				train_beta_filled[i,j] = mean
	return train_beta_filled

# Fills in nans with uniformly distributed random number
# Pass in beta array for training bed
# Returns complete array
def fill_rand(train_beta):
        train_beta_filled = train_beta
        if (train_beta.shape[1] != 33):
                print "Array must have 33 samples"
                sys.exit(2)
        for i in range(0, len(train_beta)):
                for j in range(0,33):
                        if (np.isnan(train_beta[i,j])):
                                train_beta_filled[i,j] = np.random.uniform(0,1)
        return train_beta_filled

# Fills in nans with mean over k nearest CpG sites
# Pass in beta array for training bed
# Returns complete array
def fill_neighbors(sites, train_beta, k):
	kdtree = cKDTree(sites, leafsize=100)
	train_beta_filled = train_beta
        if (train_beta.shape[1] != 33):
                print "Array must have 33 samples"
                sys.exit(2)
	if (sites.shape[1] != 2):
                print "Array must have 2 columns corresponding to start and end sites"
                sys.exit(2)
	for i in range(0,33):
		sample = train_beta[:, i]
        	for j in range(0, len(train_beta)):
                        if (np.isnan(sample[j])):
				(dist, j) = kdtree.query(sites, k)
				train_beta_filled[j, i] = np.mean(sample[j])	
        return train_beta_filled

def main(argv):
        parser = OptionParser()
        parser.add_option("-p", "--path", dest="path", help='read bed data fom PATH', metavar='PATH')
        (options, args) = parser.parse_args()
        path = options.path
        print "PATH = " + path
        train = read_bed_dat_train(path + TRAIN_FNAME)
	sites = np.transpose(np.array((train['Start'], train['End'])))
	beta = fill_mean(train['Beta'])
	print "Sample mean replacement: %s" % beta[0]
	beta = fill_rand(train['Beta'])
	print "Random replacement: %s" % beta[0]
	beta = fill_neighbors(sites, train['Beta'], 10)
	print "k=10 neighbors replacement: %s" % beta[0]

if __name__ == '__main__':
    main(sys.argv[1:])
