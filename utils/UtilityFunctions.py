'''
Created on Mar 23, 2015

@author: jonathanshor
@author: epaul626
'''
import numpy as np
import sys, time
from optparse import OptionParser

PRE_FNAME = "intersected_final_chr"
TRAIN_FNAME = "_cutoff_20_train_revised.bed"
TRAIN_INAME = "_cutoff_20_train_revised_island.bed"
SAMPLE_FNAME = "_cutoff_20_sample.bed"
SAMPLE_INAME = "_cutoff_20_sample_island.bed"
TEST_FNAME = "_cutoff_20_test.bed"
TEST_INAME = "_cutoff_20_test_island.bed"
INT_TR_FNAME = "_train_full_feats.bed"
INT_S_FNAME = "_sample_full_feats.bed"
All_TE_FNAME = "all_test.bed"

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

def insert_island_data(bed, islandFile, isTrain=False):
# Accepts nadrray and a full path to the island bed to read
# Returns new ndarray with new 'Island' boolean field
	new = np.empty(bed.shape, dtype = bed.dtype.descr + [('Island', np.bool)])
	for col in bed.dtype.names:
		new[col] = bed[col]
	new['Island'] = False
	dtype=[('Chrom', np.str_, 4), ('Start', np.int32), ('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32), ('450k', np.int8)]
	if isTrain:
		dtype=[('Chrom', np.str_, 4), ('Start', np.int32), ('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32, (33)), ('450k', np.int8)]
	inFile = np.loadtxt(islandFile, dtype)
	for i in range(len(inFile)):
		new['Island'][np.where(new['Start'] == inFile['Start'][i])] = True
	return new

def read_bed_dat_sample(mypath, chrom=1, addIsland=False, Island=False):
# Accepts path to location of PRE_FNAME + chrom + SAMPLE_FNAME
# chrom defaults to 1 
	if Island and chrom==1:
		return np.loadtxt(mypath + PRE_FNAME + str(chrom) + SAMPLE_INAME, \
						dtype=[('Chrom', np.str_, 4), ('Start', np.int32),
							('End', np.int32), ('Strand', np.str_, 1), \
							('Beta', np.float32), ('450k', np.int8)])

	bed = np.loadtxt(mypath + PRE_FNAME + str(chrom) + SAMPLE_FNAME, dtype=[('Chrom', np.str_, 4), ('Start', np.int32), \
									('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32), ('450k', np.int8)])
	if addIsland and chrom==1:
		return insert_island_data(bed, mypath + PRE_FNAME + str(chrom) + TEST_INAME)
	else:
		return bed

def read_bed_dat_test(mypath, chrom=1, addIsland=False, Island=False):
# Accepts path to location of PRE_FNAME + chrom + TEST_FNAME
# chrom defaults to 1 
	if Island and chrom==1:
		return np.loadtxt(mypath + PRE_FNAME + str(chrom) + TEST_INAME, \
						dtype=[('Chrom', np.str_, 4), ('Start', np.int32), ('End', np.int32), \
							('Strand', np.str_, 1), ('Beta', np.float32), ('450k', np.int8)])

	bed = np.loadtxt(mypath + PRE_FNAME + str(chrom) + TEST_FNAME, dtype=[('Chrom', np.str_, 4), ('Start', np.int32), \
									('End', np.int32), ('Strand', np.str_, 1), \
									('Beta', np.float32), ('450k', np.int8)])
	if addIsland and chrom==1:
		return insert_island_data(bed, mypath + PRE_FNAME + str(chrom) + TEST_INAME)
	else:
		return bed

def read_bed_dat_feat(mypath, chrom=1, ftype='train'):
# Read featured bed file from mypath of type given by ftype
# If chrom == 0, read file containing full genome data
# Return record dtype nparray
# Note - removed 'test' option since this does not have full features
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
	
	dt = [('Chrom', np.str_, 4), ('Start', np.int32), ('End', np.int32), \
		('Strand', np.str_, 1), ('Beta', np.float32, (beta_len)), ('450k', np.byte)]
	if ftype == 'train' or ftype == 'sample':
		dt += [('Exon', np.byte), ('DHS', np.int8), ('CGI', np.byte)]

 	if ftype == 'train':
		cols = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,-3,-1)
	if ftype == 'sample':
		cols = (0,1,2,3,4,5,6,-3,-1)

	return np.loadtxt(mypath + chrom_s + fname, dtype = dt, usecols=cols)
# 					dtype=[('Chrom', np.str_, 4), ('Start', np.int32), ('End', np.int32), \
# 						('Strand', np.str_, 1), ('Beta', np.float32, (beta_len)), ('450k', np.int8), \
# 						('Exon', np.byte), ('DHS', np.int8), ('CGI', np.byte)])


def read_bed_dat_train_feat(mypath, chrom=1):
# read_bed_dat_feat wrapper for train beds
	return read_bed_dat_feat(mypath, chrom, 'train')
def read_bed_dat_test_feat(mypath, chrom=1):
# read_bed_dat_feat wrapper for test beds
	return read_bed_dat_feat(mypath, chrom, 'test')
def read_bed_dat_sample_feat(mypath, chrom=1):
# read_bed_dat_feat wrapper for sample beds
	return read_bed_dat_feat(mypath, chrom, 'sample')

def read_bed_dat_train(mypath, chrom=1, addIsland=False, Island=False):
# Accepts path to location of PRE_FNAME + chrom + TRAIN_FNAME
# chrom defaults to 1 
	if Island and chrom==1:
		return np.loadtxt(mypath + PRE_FNAME + str(chrom) + TRAIN_INAME, \
						dtype=[('Chrom', np.str_, 4), ('Start', np.int32), ('End', np.int32), \
							('Strand', np.str_, 1), ('Beta', np.float32, (33)), ('450k', np.int8)])

	bed = np.loadtxt(mypath + PRE_FNAME + str(chrom) + TRAIN_FNAME, \
					dtype=[('Chrom', np.str_, 4), ('Start', np.int32), ('End', np.int32), \
						('Strand', np.str_, 1), ('Beta', np.float32, (33)), ('450k', np.int8)])
	if addIsland and chrom==1:
		return insert_island_data(bed, mypath + PRE_FNAME + str(chrom) + TEST_INAME)
	else:
		return bed


def storePreds(path, yHats, paras, start_time):
# Stores ndarray to txt in path+"predictions_"+paras+"_csec="+time.time()-start_time+".txt"
# paras also written as a comment to first line
# yHats will be formatted as floats 
	np.savetxt(path+"predictions_"+str(paras)+"_csec="+str(int(time.time()-start_time))+".txt", yHats, fmt="%f", header=paras)

def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--path", dest="path", help='read bed data from PATH', metavar='PATH')
	(options, _args) = parser.parse_args()     
	path = options.path

	train_feat_all = read_bed_dat_feat(path, chrom=1)
	sample_feat_all = read_bed_dat_feat(path, chrom=1, ftype='sample')
	print "train_feat_all[0] = %s" % train_feat_all[0]
        print "sample_feat_all[0] = %s" % sample_feat_all[0]

	train = read_bed_dat_train(path)
	sample = read_bed_dat_sample(path)
	train_feat = read_bed_dat_train_feat(path)
	print "train_feat_all[0] = %s" % train_feat_all[0]
	print "sample_feat_all[0] = %s" % sample_feat_all[0]
	print "test_feat_all[0] = %s" % test_feat_all[0]
	print "train_feat[0] = %s" % train_feat[0]
	print "sample[0] = %s" % sample[0]
	print "train[0] = %s" % train[0]
	print "len(train) = %s" % len(train)
	print "train['Beta'][0] = %s" % train['Beta'][0]

if __name__ == '__main__':
	main(sys.argv[1:])
