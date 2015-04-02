'''
Created on Mar 25, 2015

@author: epaul626
'''
import numpy as np
import sys, time
from optparse import OptionParser
sys.path.append('../utils/')
from UtilityFunctions import *

# Uses empirical mean of training data at CpG site to predict NaN data in the 
# sample file. 
# Pass in training beta array (width = 33) and sample beta vector of same length
# Returns sample vector with filled in NaN values. 
def predict_simple_mean(train, sample):
	sample_predicted = sample
	if (train.shape[1] != 33):
		print "Array must have 33 samples"
		sys.exit(2)
	if (train.shape[0] != sample.shape[0]):
		print "Train and sample must have same length"
		sys.exit(2)
	# Replace sample NaNs with mean of training samples at CpG site
	for i in range(0, len(sample)):
		if (np.isnan(sample[i])):
			sample_predicted[i] = np.mean(train[i][~np.isnan(train[i])])
	return sample_predicted

def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--path", dest="path", help='read bed data fom PATH', metavar='PATH')
	(options, args) = parser.parse_args()
	path = options.path

	time_start = time.time()
	train = read_bed_dat_train(path)
	sample = read_bed_dat_sample(path)
	# Pull Testing data
	test = read_bed_dat_test(path)
	gTruth = test['Beta'][sample['450k']==0]
	gTruth = gTruth[~np.isnan(gTruth)]
	# Predict using simple mean 
	pred_time = time.time()
	predict = predict_simple_mean(train['Beta'], sample['Beta'])
	predict = predict[sample['450k']==0]
	predict = predict[~np.isnan(gTruth)]
	print "Empirical Mean Prediction"
	(r2, RMSE) = calc_r2_RMSE(predict, gTruth)
	print "r2: %s" % r2
	print "RMSE: %s" % RMSE
	print "Prediction time: %s" % (time.time() - pred_time)
	print "Elapsed time: %s" % (time.time() - time_start)

if __name__ == '__main__':
    main(sys.argv[1:])

