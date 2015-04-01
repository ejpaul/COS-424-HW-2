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
	X = np.zeros((len(train_beta), 37))
	Y = np.zeros((len(train_beta)))
	# Feature 1: CpG start site
	X[:,0] = sites[:,0]
	# Feature 2: CpG end site
	X[:,1] = sites[:,1]
	for i in range(0, len(train_beta)):
		# Find 2 nearest neighbors to site i
		indices = find_neighbors(sites, i)
		index1 = indices[0]
		index2 = indices[1]
		# Feature 3: beta at neighbor 1
		X[i, 2] = train_beta[i+index1, j]
		# Feature 4: distance to neighbor 1
		# Feature 5: beta at neighbor 2
		# Feature 6: distance to neighbor 2
		X[i, 3] = train_beta[i+index2, j]
		# Features 5-37: 33 sample beta values at CpG site
		for j in range(0, 33):
			X[i, 3+j] = train_beta[i, j]


	return (X, Y)

# Find nearest 2 sites to current CpG site
# Returns list of indices of 2 nearest neighbors
# Where the index in 'sites' is i + (index1, index2)
def find_neighbors(sites, i):
	start_curr = sites[i, 0]
	end_curr = sites[i, 1]
	# Look to left
	if (i > 0):
		d_l1 = start_curr - sites[i-1, 1]
	else: 
		return ([1, 2], [d_
	if (i > 1):
		d_l2 = start_curr - sites[i-2, 1]
	else:
		d_l2 = -1
	# Look to right
	if (i < len(sites)-2):
		d_r1 = sites[i-1, 0] - end_curr
	else:
		return [-1, -2]
	if (i < len(sites)-3):
		d_r2 = sites[i-2, 0] - end_curr
	else:
		d_r2 = -1
	# Compare distances
	if (d_l1 < d_r1 and d_l1 != -1):
		site_1 = -1
		if (d_l2 < d_r1):
			site_2 = -2
		else:
			 site_2 = 1
		return [site_1, site_2] 
	site_1 = 1
	if (d_r2 < d_l1 and d_r2 != -1):
		site_2 = 2
	else:
		site_2 = -1	
	return [site_1, site_2]
 
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
