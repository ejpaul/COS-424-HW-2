'''
Created on Mar 24, 2015

@author: epaul626

Utility function for filling in nans for training
'''
import numpy as np
import sys
from optparse import OptionParser
from scipy.spatial import cKDTree
from UtilityFunctions import read_bed_dat_train 

def fill_mean(train_beta):
# Fills in nans by averaging over 32 samples in training bed
# Pass in beta array for training bed
# Returns complete array
	train_beta_filled = train_beta
	for i in range(0, len(train_beta)):
		# Replace nan samples with mean of non-nan samples
		train_beta_filled[i][np.isnan(train_beta[i])] = np.mean(train_beta[i][~np.isnan(train_beta[i])])
	return train_beta_filled

def fill_rand(train_beta):
# Fills in nans with uniformly distributed random number
# Pass in beta array for training bed
# Returns complete array
	train_beta_filled = train_beta
	train_beta_filled[np.isnan(train_beta)]=np.random.uniform(0,1,train_beta[np.isnan(train_beta)].shape)
	return train_beta_filled

def fill_neighbors(site, train_beta, k):
# Fills in nans with mean over k nearest CpG sites
# Pass in beta array for training bed
# Returns complete array
	kdtree = cKDTree(site, leafsize=100)
	train_beta_filled = train_beta
	if (train_beta.shape[1] != 33):
		print "Array must have 33 samples"
		sys.exit(2)
	if (site.shape[1] != 1):
		print "Array must have 1 columns corresponding to start site"
		sys.exit(2)
	for i in range(0,33):
		sample = train_beta[:, i]
        	for j in range(0, len(train_beta)):
			# If NaN found, query the tree for the k nearest points, with an upper bound of 
			# 100 bp. m is an array of ints corresponding to location of the neighbors within
			# position[]
			if (np.isnan(sample[j])):
				(dist, j) = kdtree.query(site, k)
				train_beta_filled[j, i] = np.mean(sample[j])	
	return train_beta_filled

def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--path", dest="path", help='read bed data fom PATH', metavar='PATH')
	(options, args) = parser.parse_args()
	path = options.path
	print "PATH = " + path
	train = read_bed_dat_train(path)
	site = train['Start']

	print "Pre-sample mean replacement: %s" % train['Beta'][0]
	beta = fill_mean(train['Beta'].copy())
	print "Sample mean replacement: %s" % beta[0]
	print "Pre-random replacement: %s" % train['Beta'][0]
	beta = fill_rand(train['Beta'].copy())
	print "Random replacement: %s" % beta[0]
	print "Pre-neighbors replacement: %s" % train['Beta'][0]
	beta = fill_neighbors(site, train['Beta'].copy(), 10)
	print "k=10 neighbors replacement: %s" % beta[0]

if __name__ == '__main__':
    main(sys.argv[1:])
