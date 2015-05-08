'''
Created on May 6, 2015

@author: epaul626
'''
import numpy as np
import sys, time
from scipy.stats import pearsonr
from optparse import OptionParser
from UtilityFunctions import read_bed_dat_feat
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
	upstream = np.empty_like(train)
	upstream[1:] = train[:-1]
	upstream[0] = np.zeros(upstream[0].shape)
	return upstream

def get_downstream(train):
# Expects array of 'Beta'
# Returns array of downstream neighbor values for each site
# of same shape
	downstream = np.empty_like(train)
	downstream[:-1] = train[1:]
	downstream[-1] = np.zeros(downstream[-1].shape)
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

# Stores correlation values for upstream and downstream neighbor
# for all training data in 'path'
# stores in 'train_corr.txt' in 'path' 
# Format: corr_upstream	corr_downstream 
def store_corr(path, chrm):
	train = read_bed_dat_train(path, chrm)
        beta = train['Beta']
        neighb = get_downstream(beta)
        corr_down = corr_calc(beta, neighb)
	neighb = get_upstream(beta)
	corr_up = corr_calc(beta, neighb)
	dat = np.column_stack((corr_up, corr_down))
	np.savetxt('chr' + str(chrm) + '_train_corr.txt', dat)

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
