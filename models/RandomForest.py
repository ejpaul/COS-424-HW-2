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

# Produces 'X' feature array of nearest neighbor beta values that will be fed
# into regressor
# Beta is an array of shape (# of CpG sites, nsamples)
# Sites is an array of shape (# of CpG sites, 2) corresponding to 
# start and ending bp of site
def feat_neighbors(sites, train_beta):
	X = np.zeros(len(train_beta)*33, 5)
	Y = np.zeros(len(train_beta)*33, 1)
	# Iterate over CpG sites
	for i in range(0, len(train_beta)):
		indices = find_neighbors(sites, i)[0]
		index1 = indices[0]
		index2 = indices[1]
		for j in range(0, 33):
			k = i*33 + j
			# Feature 1: nearest neighbor
			X[k, 0] = train_beta[j, i+index1]
			# Feature 2: second nearest neighbor
			X[k, 1] = train_beta[j, i+index2]
			# Feature 3: CpG start site
			X[k, 2] = sites[i, 0]
			# Feature 4: CpG end site
			X[k, 3] = sites[i, 1]
			# Feature 5: Sample number
			X[k, 4] = j
			# Beta value at sample site
			Y[k] = train_beta[j, i]
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
		return [1, 2]
	if (i > 1):
		d_l2 = start_curr - sites[i-2, 1]
	else:
		d_l2 = -1
	# Look to right
	if (i < len(train_beta)-2):
		d_r1 = sites[train_beta-1, 0] - end_curr
	else:
		return [-1, -2]
	if (i < len(train_beta)-3):
		d_r2 = sites[train_beta-2, 0] - end_curr
	else:
		d_r2 = -1
	# Compare
	if (d_l1 < d_r1):
		site_1 = -1
		if (d_l2 < d_r1):
			site_2 = -2
		else:
			 site_2 = 1
		return [site_1, site_2] 
	site_1 = 1
	if (d_r2 < d_l1):
		site_2 = 2
	else:
		site_2 = -1	
	return [site_1, site_2]
 
def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--path", dest="path", help='read bed data fom PATH', metavar='PATH')
	(options, args) = parser.parse_args()
	path = options.path
        train = read_bed_dat_train(path + TRAIN_PATH)
        sites = np.transpose(np.array((train['Start'], train['End'])))
# Fill in NaNs in training with mean over 33 samples in training bed
        train_beta = fill_mean(train['Beta'])
# Produce feature array 'X' and vector of beta values 'Y'
	(X, Y) = feat_neighbors(sites, train_beta)
# Initialize regressor with default parameters
	model = RandomForestRegressor()
# Fit regressor using training data
	model.fit(X,Y)

if __name__ == '__main__':
    main(sys.argv[1:])
