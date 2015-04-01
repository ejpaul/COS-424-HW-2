'''
created on Mar 29, 2015

@author: epaul626
'''
import numpy as np
import sys
from optparse import OptionParser
from sklearn.ensemble import RandomForestRegressor
sys.path.append('../utils/')
from UtilityFunctions import *
from FillData import *

# Produces 'X' feature array of nearest neighbor beta values that will be fed
# into regressor
# Beta is an array of shape (# of CpG sites, nsamples)
# Sites is an array of shape (# of CpG sites, 2) corresponding to 
# start and ending bp of site
def feat_neighbors(sites, train_beta, sample):
	X = np.zeros((len(train_beta), 38))
	Y = np.zeros((len(train_beta)))
	# Feature 1: CpG start site
	X[:,0] = sites[:,0]
	for i in range(0, len(train_beta)):
                # Find 2 nearest neighbors to site i
                (indices,distance)  = find_neighbors(sites, i)
                index1 = indices[0]
                index2 = indices[1]
		for j in range(0, 33):
			# Feature 2: beta at neighbor 1
			X[i, 1] = train_beta[i+index1, j]
			# Feature 3: distance to neighbor 1
			X[i, 2] = distance[0]
			# Feature 4: beta at neighbor 2
			X[i, 3] = train_beta[i+index2, j]
			# Feature 5: distance to neighbor 2
			X[i, 4] = distance[1]
			# Features 6-38: 33 sample beta values at CpG site
			X[i, 5+j] = train_beta[i, j]
	return (X, Y)

# Find nearest 2 sites to current CpG site
# Returns list of indices of 2 nearest neighbors
# Where the index in 'sites' is i + (index1, index2)
def find_neighbors(sites, i):
	start_curr = sites[i, 0]
	end_curr = sites[i, 1]
	d_l1 = start_curr - sites[i-1, 1]
	d_l2 = start_curr - sites[i-2, 1]
	d_r1 = sites[i-1, 0] - end_curr
	d_r2 = sites[i-2, 0] - end_curr
	# Look to left
	if (i == 0):
		return ([1, 2], [d_r1, d_r2])
	if (i == 1):
		d_l2 = -1
	# Look to right
	if (i == len(sites)-1):
		return ([-1, -2], [d_l1, d_l2])
	if (i == len(sites)-2):
		d_r2 = -1
	# Compare distances
	if (d_l1 < d_r1 and d_l1 != -1):
		if (d_l2 < d_r1):
			return ([-1, -2], [d_l1, d_l2])
		else:
			return ([-1, 1], [d_l1, d_r1])
	if (d_r2 < d_l1 and d_r2 != -1):
		return ([1, 2], [d_r1, d_r2])
	else:
		return ([1, -1], [d_r1, d_l1])
 
def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--path", dest="path", help='read bed data fom PATH', metavar='PATH')
	(options, args) = parser.parse_args()
	path = options.path
        train = read_bed_dat_train(path)
	sample = read_bed_dat_sample(path)
        sites = np.transpose(np.array((train['Start'], train['End'])))
# Fill in NaNs in training with mean over 33 samples in training bed
        train_beta = fill_mean(train['Beta'])
# Produce feature array 'X' and vector of beta values 'Y'
	(X, Y) = feat_neighbors(sites.copy(), train_beta.copy(), sample.copy())
# Initialize regressor with default parameters
	model = RandomForestRegressor()
# Fit regressor using training data
	model.fit(X, Y)

if __name__ == '__main__':
    main(sys.argv[1:])
