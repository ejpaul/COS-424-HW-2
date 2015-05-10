'''
created on Mar 29, 2015

@author: epaul626
'''
import numpy as np
import sys, time
from optparse import OptionParser
from sklearn.ensemble import RandomForestRegressor
sys.path.append('../utils/')
import UtilityFunctions as uf
from FillData import fill_mean, fill_neighbors

def feat_neighbors(sites, train_beta, train_extras, sample, test):
# Produces 'X' feature array of nearest neighbor beta values that will be fed
# into regressor
# Beta is an array of shape (# of CpG sites, nsamples)
# Sites is an array of shape (# of CpG sites, ) corresponding to 
# the start of each bp site
# train_extras: array of shape (# of CpG sites, 3) which
	# Adds 6 extra features: Exon, DHS, CGI, GC_100, GC_400, GC_1000
	num_extras = train_extras.shape[1]
	X = np.zeros((len(train_beta), 38+num_extras))
	for i in range(0, len(train_beta)):
		# Feature 1: CpG start site
		X[i,0] = sites[i]
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
	# Features 39+: Extra features
	X[:,-num_extras:] = train_extras
	# Predict on feature set not on 450k chip
	Xstar = X[sample['450k']==0]
	gTruth = test['Beta'][sample['450k']==0]
	Xstar = Xstar[~np.isnan(gTruth)]
	gTruth = gTruth[~np.isnan(gTruth)]

	# Only train on non-NaN values of Y
	X = X[~np.isnan(sample['Beta'])]
	Y = sample['Beta'][~np.isnan(sample['Beta'])]
	return (X, Y, Xstar, gTruth)

def feat_neighbors_island(sites, train_beta, train_extras, sample, test, island):
# Produces 'X' feature array of nearest neighbor beta values that will be fed
# into regressor
# Beta is an array of shape (# of CpG sites, nsamples)
# Sites is an array of shape (# of CpG sites, ) corresponding to
# the start of each bp site
# train_extras: array of shape (# of CpG sites, 3) which
        # Adds 6 extra features: Exon, DHS, CGI, GC_100, GC_400, GC_1000
        num_extras = train_extras.shape[1]
        X = np.zeros((len(train_beta), 38+num_extras))
        for i in range(0, len(train_beta)):
                # Feature 1: CpG start site
                X[i,0] = sites[i]
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
        # Features 39+: Extra features
        X[:,-num_extras:] = train_extras
	# Select on islands
	X = X[island == 1]
        sample = sample[island == 1]
        test = test[island == 1]
        # Predict on feature set not on 450k chip
        Xstar = X[sample['450k']==0]
        gTruth = test['Beta'][sample['450k']==0]
        Xstar = Xstar[~np.isnan(gTruth)]
        gTruth = gTruth[~np.isnan(gTruth)]

        # Only train on non-NaN values of Y
        X = X[~np.isnan(sample['Beta'])]
        Y = sample['Beta'][~np.isnan(sample['Beta'])]
        return (X, Y, Xstar, gTruth)


def feat_neighbors_thresh(sites, train_beta, train_extras, sample, test, corr, thresh):
# Produces 'X' feature array of nearest neighbor beta values that will be fed
# into regressor
# Beta is an array of shape (# of CpG sites, nsamples)
# Sites is an array of shape (# of CpG sites, ) corresponding to
# the start of each bp site
# train_extras: array of shape (# of CpG sites, 3) which
        # Adds 6 extra features: Exon, DHS, CGI, GC_100, GC_400, GC_1000
        num_extras = train_extras.shape[1]
        X = np.zeros((len(train_beta), 38+num_extras))
        for i in range(0, len(train_beta)):
                # Feature 1: CpG start site
                X[i,0] = sites[i]
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
# Features 39+: Extra features
        X[:,-num_extras:] = train_extras
	# Train on sites with > thresh correlation value
	# Trim off first and last row
	corr = corr[1:-1]
	X = X[corr[:,0] > thresh]
	sample = sample[corr[:,0] > thresh]
	test = test[corr[:,0] > thresh]
	corr = corr[corr[:,0] > thresh]
	X = X[corr[:,1] > thresh]
	sample = sample[corr[:,1] > thresh]
	test = test[corr[:,1] > thresh]
	corr = corr[corr[:,1] > thresh]

        # Predict on feature set not on 450k chip
        Xstar = X[sample['450k']==0]
        gTruth = test['Beta'][sample['450k']==0]
        Xstar = Xstar[~np.isnan(gTruth)]
        gTruth = gTruth[~np.isnan(gTruth)]

        # Only train on non-NaN values of Y
        X = X[~np.isnan(sample['Beta'])]
        Y = sample['Beta'][~np.isnan(sample['Beta'])]
        return (X, Y, Xstar, gTruth)

# Find nearest 2 sites to current CpG site
# Returns list of indices of 2 nearest neighbors
# Where the index in 'sites' is i + (index1, index2)
def find_neighbors(sites, i):
	curr = sites[i]
	if (i > 0):
		d_l1 = curr - sites[i-1]
	if (i > 1):
		d_l2 = curr - sites[i-2]
	if (i < len(sites)-1):
		d_r1 = sites[i+1] - curr
	if (i < len(sites)-2):
		d_r2 = sites[i+2] - curr
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
	parser.add_option("-i", dest="island", action="store_true", default=False)
	parser.add_option("-n", dest="fill_neighb", action="store_true", default=False)
	parser.add_option("-m", dest="fill_mean", action="store_true", default=False)
	parser.add_option("-c", dest="chroms", type='int', default=1)
	parser.add_option("-g", dest="GCwindow", type='int', default=0, help='100,400, or 1000. Else all will be included.')
	parser.add_option("-t", dest="thresh", type='float', default=0.0, help='Correlation value threshold.')
	parser.add_option("-e", dest="estimators", type='int', default=20, help='Number of trees in the forest.')
	(options, _args) = parser.parse_args()
	path = options.path

# 	fill_n = options.fill_neighb
	fill_m = options.fill_mean
	chroms = options.chroms
	trees = options.estimators
	thresh = options.thresh
	island = options.island
	start_time = time.time()
	if fill_m:
		fill_str = 'mean'
	else:
		fill_str = 'neigh'
	paras = "chr=%s_fill=%s_GC=%s_thrsh=%s_trees=%s" % (chroms,fill_str, options.GCwindow, thresh, trees)

	# Read in full feature data
	train = uf.read_bed_dat_feat(path, chrom=chroms, ftype='train')
	sample = uf.read_bed_dat_feat(path, chrom=chroms, ftype='sample')
	test = uf.read_bed_dat_feat(path, chrom=chroms, ftype='test')
	print "Beds read %s" % (time.time() - start_time)
	sites = train['Start']
# Fill in NaNs in training with mean over 33 samples in training bed
	train_extras = np.c_[train['Exon'], train['DHS'], train['CGI']]
	if options.GCwindow:
		try:
			train_extras = np.c_[train_extras, train['GC_%s' % str(options.GCwindow)]]
		except ValueError:
			train_extras = np.c_[train_extras, train['GC_100'], train['GC_400'], train['GC_1000']]
	else:
		train_extras = np.c_[train_extras, train['GC_100'], train['GC_400'], train['GC_1000']]
	train_extras = np.c_[train_extras, uf.read_corrs(path, chroms)]
	print "train_extras[0]: %s, shape: %s" % (train_extras[1:-1][0], train_extras[1:-1].shape)
	if fill_m:
		train_beta = fill_mean(train['Beta'])
	else:
		train_beta = fill_neighbors(sites, train['Beta'], 10)
# Trim first and last row off to eliminate nan values in corrs
# Produce feature array 'X' and vector of beta values 'Y'
# Produce feature array 'X*' to predict on and ground truth beta values 'Y*'
	if thresh:
		(X, Y, Xstar, Ystar) = feat_neighbors_thresh(sites[1:-1].copy(), train_beta[1:-1].copy(), \
							train_extras[1:-1].copy(), sample[1:-1].copy(), test[1:-1].copy(), uf.read_corrs(path, chroms), thresh)
	elif island:
		(X, Y, Xstar, Ystar) = feat_neighbors_island(sites[1:-1].copy(), train_beta[1:-1].copy(), \
                                                        train_extras[1:-1].copy(), sample[1:-1].copy(), test[1:-1].copy(), train['CGI'])
	else:
		(X, Y, Xstar, Ystar) = feat_neighbors(sites[1:-1].copy(), train_beta[1:-1].copy(), \
							train_extras[1:-1].copy(), sample[1:-1].copy(), test[1:-1].copy())	
	print "X[0]: %s" % X[0]
	full_feat_start = time.time()
	# Initialize regressor with default parameters
	model = RandomForestRegressor(oob_score=True, n_estimators=trees)
# Fit regressor using training data
	model.fit(X, Y)
# Predict on Xstar values 
	Yhat = model.predict(Xstar)
# Calculate r2 and RMSE
	(r2, RMSE) = uf.calc_r2_RMSE(Yhat, Ystar)
	print "Runtime: %f" % (time.time()-full_feat_start)
	print "RandomForest Runtime: %f" % (time.time()-start_time)
	print str(model.feature_importances_)
	print "oob: " + str(model.oob_score_)
	print "r2 : %f" % (r2)
	print "RMSE: %f" % (RMSE)
	print "Total Runtime: %f" % (time.time()-start_time)
	uf.storePreds(path, Yhat, paras, start_time)

if __name__ == '__main__':
	main(sys.argv[1:])
