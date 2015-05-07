'''
Created on May 6, 2015

@author: epaul626
'''
import numpy as np
import sys, time
from scipy.stats import pearsonr
from optparse import OptionParser
from UtilityFunctions import read_bed_dat_train
import pylab as P

PRE_FNAME = "intersected_final_chr"
TRAIN_FNAME = "_cutoff_20_train_revised.bed"
TRAIN_INAME = "_cutoff_20_train_revised_island.bed"
SAMPLE_FNAME = "_cutoff_20_sample.bed"
SAMPLE_INAME = "_cutoff_20_sample_island.bed"
TEST_FNAME = "_cutoff_20_test.bed"
TEST_INAME = "_cutoff_20_test_island.bed"

def get_upstream(train):
# Expects array of 'Beta'
# Returns array of upstream neighbor values for each site
# of same shape
	upstream = np.roll(train, 1)
	upstream = upstream[0:-2]
	# Insert 0s for last row
	np.insert(upstream, -1, 0)
	return upstream

def get_downstream(train):
# Expects array of 'Beta'
# Returns array of downstream neighbor values for each site
# of same shape
	downstream = np.roll(train, len(train)-1)
	downstream = downstream[1:-1]
	# Insert 0s for first row
	np.insert(downstream, 0, 0)
	return downstream

def corr_calc(train, neighb):
# Expects array of 'Beta' and array of neighbors of same shape
# for correlation calculation
	# For each position, calculate correlation
	# of beta with neighboring beta
	
	# Vector of correlation coefficients for each site	
	corr = np.zeros((len(train),))	
	for i in range(0, len(corr)):
		corr[i] = pearsonr(train[i], neighb[i])[0]
	return corr

def corr_neighb(train, neighb):
# For each position, calculate correlation of beta with neighbor 
# Uses 0.5 cutoff for correlation value of feature
	corr = corr_calc(train, neighb)
	feat = corr*neighb
	return feat	
	
def dist_neighb(train, neighb):
# Weights beta value at neighbor by 1/dist for input as feature vector
# Train and neighb are entire csv file, not just betas
	start_train = train['Start']
	start_neighb = neighb['Start']
	dist = np.abs(start_train - start_neighb)
	weight = 1./dist
	feat = neighb['Beta']*weight
	return feat	

def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
	(options, _args) = parser.parse_args()     
	path = options.path
	train = read_bed_dat_train(path)
	beta = train['Beta']
	
	neighb = get_downstream(beta)
	corr = corr_calc(beta, neighb)

	P.hist(corr, range=(-1,1),bins=30)

	P.figure()
	start = train['Start']
	P.plot(start, corr, 'r.')
	P.show()

if __name__ == '__main__':
	main(sys.argv[1:])
