'''
Created on Mar 23, 2015

@author: jonathanshor
@author: epaul626
'''
import numpy as np
import sys, time
from optparse import OptionParser

PRE_FNAME = "intersected_final_"
TRAIN_FNAME = "_cutoff_20_train_revised.bed"
TRAIN_INAME = "_cutoff_20_train_revised_island.bed"
SAMPLE_FNAME = "_cutoff_20_sample.bed"
SAMPLE_INAME = "_cutoff_20_sample_island.bed"
TEST_FNAME = "_cutoff_20_test.bed"
TEST_INAME = "_cutoff_20_test_island.bed"
INT_TR_FNAME = "_train_full.bed"
INT_S_FNAME = "_sample_full.bed"
All_TE_FNAME = "_test.bed"
PREDSUM_FNAME = 'PredictionSummary.txt'

def calc_r2_RMSE(preds, gTruth, intercept = 0):
# Expects two 'Beta' both of size nX1
# Accepts an intercept and calc's accordingly
# Returns (Coeff of Determination (r^2), RMSE)
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

def read_bed_dat_train(mypath, chrom=1, addIsland=False, Island=False):
# Accepts path to location of PRE_FNAME + chrom + TRAIN_FNAME
# chrom defaults to 1 
	if chrom == 0:
		chrom_s = 'all'
	else:
		chrom_s = 'chr' + str(chrom)

	if Island and chrom==1:
		return np.loadtxt(mypath + PRE_FNAME + chrom_s + TRAIN_INAME, \
						dtype=[('Chrom', np.str_, 5), ('Start', np.int32), ('End', np.int32), \
							('Strand', np.str_, 1), ('Beta', np.float32, (33)), ('450k', np.int8)])

	bed = np.loadtxt(mypath + PRE_FNAME + chrom_s + TRAIN_FNAME, \
					dtype=[('Chrom', np.str_, 5), ('Start', np.int32), ('End', np.int32), \
						('Strand', np.str_, 1), ('Beta', np.float32, (33)), ('450k', np.int8)])
	return bed

def read_bed_dat_sample(mypath, chrom=1, addIsland=False, Island=False):
# Accepts path to location of PRE_FNAME + chrom + SAMPLE_FNAME
# chrom defaults to 1 
	if chrom == 0:
		chrom_s = 'all'
	else:
		chrom_s = 'chr' + str(chrom)

	if Island and chrom==1:
		return np.loadtxt(mypath + PRE_FNAME + chrom_s + SAMPLE_INAME, \
						dtype=[('Chrom', np.str_, 5), ('Start', np.int32),
							('End', np.int32), ('Strand', np.str_, 1), \
							('Beta', np.float32), ('450k', np.int8)])

	bed = np.loadtxt(mypath + PRE_FNAME + chrom_s + SAMPLE_FNAME, dtype=[('Chrom', np.str_, 5), ('Start', np.int32), \
									('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32), ('450k', np.int8)])
	return bed

def read_bed_dat_test(mypath, chrom=1, addIsland=False, Island=False):
# Accepts path to location of PRE_FNAME + chrom + TEST_FNAME
# chrom defaults to 1 
	if chrom == 0:
		prefname = ''
		chrom_s = 'all'
		fname = All_TE_FNAME
	else:
		prefname = PRE_FNAME
		chrom_s = 'chr' + str(chrom)
		fname = TEST_FNAME

	if Island and chrom==1:
		return np.loadtxt(mypath + PRE_FNAME + chrom_s + TEST_INAME, \
						dtype=[('Chrom', np.str_, 5), ('Start', np.int32), ('End', np.int32), \
							('Strand', np.str_, 1), ('Beta', np.float32), ('450k', np.int8)])

	bed = np.loadtxt(mypath + prefname + chrom_s + fname, dtype=[('Chrom', np.str_, 5), ('Start', np.int32), \
									('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32), ('450k', np.int8)])
	return bed

def read_bed_dat_feat(mypath, chrom=1, ftype='train'):
# Read featured bed file from mypath of type given by ftype
# If chrom == 0, read file containing full genome data
# Return record dtype nparray
# Note - removed 'test' option since this does not have full features
	if ftype == 'test':
		return read_bed_dat_test(mypath, chrom=chrom)

	if chrom == 0:
		chrom_s = 'all'
	else:
		chrom_s = 'chr' + str(chrom)

	beta_len = 1
	if ftype == 'train':
		fname = INT_TR_FNAME
		beta_len = 33
	elif ftype == 'sample':
		fname = INT_S_FNAME
	else:
		raise RuntimeError('Bad bed type name')
	
	dt = [('Chrom', np.str_, 5), ('Start', np.int32), ('End', np.int32), \
		('Strand', np.str_, 1), ('Beta', np.float32, (beta_len)), ('450k', np.byte)]
	if ftype == 'train' or ftype == 'sample':
		dt += [('Exon', np.byte), ('DHS', np.int8), ('CGI', np.byte), ('GC_100', np.float32), ('GC_400', np.float32), ('GC_1000', np.float32)]

	if ftype == 'train':
		cols = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,43,45,49,54,59)
	if ftype == 'sample':
		cols = (0,1,2,3,4,5,6,11,13,17,22,27)
	# Read in bed file
	return np.loadtxt(mypath + chrom_s + fname, dtype = dt, usecols=cols)

def read_bed_dat_train_feat(mypath, chrom=1):
# read_bed_dat_feat wrapper for train beds
	return read_bed_dat_feat(mypath, chrom, 'train')
def read_bed_dat_test_feat(mypath, chrom=1):
# read_bed_dat_feat wrapper for test beds
	return read_bed_dat_feat(mypath, chrom, 'test')
def read_bed_dat_sample_feat(mypath, chrom=1):
# read_bed_dat_feat wrapper for sample beds
	return read_bed_dat_feat(mypath, chrom, 'sample')

def read_corrs(mypath, chrom=1):
# Return neighbor correlation values in (# of CpG sites, 2) array
# Col 0: upstream neighbor, Col 1: downstream
	if chrom == 0:
		chrom_s = 'all'
	else:
		chrom_s = 'chr' + str(chrom)
	return np.loadtxt(mypath + chrom_s + "_train_corr.txt")


def storePreds(path, yHats, paras, start_time, strs = False, summary = True):
# Stores ndarray to txt in path+"predictions_"+paras+"_csec="+time.time()-start_time+".txt"
# paras also written as a comment to first line
# yHats will be formatted as floats 
	if strs:
		np.savetxt(path+"predictions_"+str(paras)+"_csec="+str(int(time.time()-start_time))+".txt", yHats, fmt="%s", header=paras)
	else:
		np.savetxt(path+"predictions_"+str(paras)+"_csec="+str(int(time.time()-start_time))+".txt", yHats, fmt="%f", header=paras)
# 	if summary:
# 		with open(path + PREDSUM_FNAME, 'a') as predf:
# 			pass
# 		predf.close()
	

def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
	(options, _args) = parser.parse_args()
	path = options.path
	
	train_feat_all = read_bed_dat_feat(path, chrom=1)
	sample_feat_all = read_bed_dat_feat(path, chrom=1, ftype='sample')
	test_feat_all = read_bed_dat_feat(path, chrom=1, ftype='test')

	corrs = read_corrs(path, chrom=1)
	print "corrs[0] = %s" % corrs[0]

	print "train_feat_all[0] = %s" % train_feat_all[0]
	print "sample_feat_all[0] = %s" % sample_feat_all[0]
	print "test_feat_all[0] = %s" % test_feat_all[0]

if __name__ == '__main__':
	main(sys.argv[1:])
